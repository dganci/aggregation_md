#include "SamplingMonitor.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <deque>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace cg {
namespace {

struct ColvarTable {
    std::vector<std::string> fields;
    std::unordered_map<std::string, std::size_t> col;
    std::vector<std::vector<double>> rows;
};

struct TransitionCounts { int transitions = 0, up = 0, down = 0; };

std::string js(bool x) { return x ? "true" : "false"; }
std::string pair_contact_name(int i, int j) { return "cn_" + std::to_string(i) + "_" + std::to_string(j); }

ColvarTable read_colvars(const std::vector<std::filesystem::path>& paths) {
    ColvarTable table;

    for (const auto& path : paths) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot read COLVAR file: " + path.string());

        std::vector<std::string> local_fields;
        bool have_fields = false;
        for (std::string line; std::getline(in, line);) {
            line = trim(line);
            if (line.empty()) continue;

            if (starts_with(line, "#!")) {
                const auto tokens = split_ws(line);
                if (tokens.size() >= 3 && tokens[1] == "FIELDS") {
                    local_fields.assign(tokens.begin() + 2, tokens.end());
                    have_fields = true;
                    if (table.fields.empty()) {
                        table.fields = local_fields;
                        for (std::size_t i = 0; i < table.fields.size(); ++i) table.col[table.fields[i]] = i;
                    } else if (local_fields != table.fields) {
                        throw std::runtime_error("COLVAR header mismatch in file: " + path.string());
                    }
                }
                continue;
            }

            if (!have_fields && table.fields.empty())
                throw std::runtime_error("COLVAR file has no '#! FIELDS' header: " + path.string());

            const auto tokens = split_ws(line);
            if (tokens.size() != table.fields.size())
                throw std::runtime_error("COLVAR row has " + std::to_string(tokens.size()) +
                                         " columns, expected " + std::to_string(table.fields.size()) +
                                         " in file: " + path.string());

            std::vector<double> row;
            row.reserve(tokens.size());
            for (const auto& t : tokens) row.push_back(std::stod(t));
            table.rows.push_back(std::move(row));
        }
    }

    if (table.rows.empty()) throw std::runtime_error("No COLVAR data rows found.");
    return table;
}

std::vector<double> column_values(const ColvarTable& table, const std::string& name) {
    const auto it = table.col.find(name);
    if (it == table.col.end()) throw std::runtime_error("Required COLVAR column not found: " + name);

    std::vector<double> out;
    out.reserve(table.rows.size());
    for (const auto& row : table.rows) out.push_back(row[it->second]);
    return out;
}

std::vector<double> moving_average(const std::vector<double>& x, int window) {
    if (window <= 1 || x.empty()) return x;

    std::vector<double> y;
    y.reserve(x.size());
    std::deque<double> q;
    double sum = 0.0;
    for (double v : x) {
        q.push_back(v);
        sum += v;
        if (static_cast<int>(q.size()) > window) {
            sum -= q.front();
            q.pop_front();
        }
        y.push_back(sum / static_cast<double>(q.size()));
    }
    return y;
}

std::vector<int> discretize_three_states(const std::vector<double>& x) {
    if (x.empty()) return {};
    const auto [mn_it, mx_it] = std::minmax_element(x.begin(), x.end());
    const auto range = *mx_it - *mn_it;
    if (range <= 1e-12) return std::vector<int>(x.size(), 1);

    const double low = *mn_it + range / 3.0;
    const double high = *mn_it + 2.0 * range / 3.0;
    std::vector<int> states;
    states.reserve(x.size());
    for (double v : x) states.push_back(v < low ? 0 : (v > high ? 2 : 1));
    return states;
}

std::vector<int> compress_states_with_residence(const std::vector<int>& states, int min_residence) {
    if (states.empty()) return {};
    min_residence = std::max(1, min_residence);

    std::vector<int> out;
    int current = states.front();
    int residence = 0;
    const auto flush = [&](int state, int count) {
        if (count >= min_residence && (out.empty() || out.back() != state)) out.push_back(state);
    };

    for (int s : states) {
        if (s == current) ++residence;
        else {
            flush(current, residence);
            current = s;
            residence = 1;
        }
    }
    flush(current, residence);
    return out;
}

TransitionCounts count_transitions(const std::vector<int>& states) {
    TransitionCounts c;
    for (std::size_t i = 1; i < states.size(); ++i) {
        const auto delta = states[i] - states[i - 1];
        if (!delta) continue;
        ++c.transitions;
        delta > 0 ? ++c.up : ++c.down;
    }
    return c;
}

int largest_connected_component(int n, const std::vector<std::vector<int>>& adj) {
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    int best = 0;

    for (int s = 0; s < n; ++s) {
        if (seen[static_cast<std::size_t>(s)]) continue;
        int size = 0;
        std::queue<int> q;
        q.push(s);
        seen[static_cast<std::size_t>(s)] = 1;

        while (!q.empty()) {
            const auto u = q.front();
            q.pop();
            ++size;
            for (int v : adj[static_cast<std::size_t>(u)]) {
                if (!seen[static_cast<std::size_t>(v)]) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    q.push(v);
                }
            }
        }
        best = std::max(best, size);
    }
    return best;
}

} // namespace

std::string SamplingMetrics::to_json(int last_chunk) const {
    std::ostringstream o;
    o << std::setprecision(10)
      << "{"
      << "\"last_chunk\":" << last_chunk << ','
      << "\"total_time_us\":" << total_time_us << ','
      << "\"n_frames\":" << n_frames << ','
      << "\"cn_min\":" << cn_min << ','
      << "\"cn_max\":" << cn_max << ','
      << "\"cn_range\":" << cn_range << ','
      << "\"rg_global_min\":" << rg_global_min << ','
      << "\"rg_global_max\":" << rg_global_max << ','
      << "\"rg_global_range\":" << rg_global_range << ','
      << "\"cn_state_transitions\":" << cn_state_transitions << ','
      << "\"cn_state_up\":" << cn_state_up << ','
      << "\"cn_state_down\":" << cn_state_down << ','
      << "\"unique_contact_patterns\":" << unique_contact_patterns << ','
      << "\"largest_cluster_unique\":" << largest_cluster_unique << ','
      << "\"largest_cluster_max\":" << largest_cluster_max << ','
      << "\"checks\":{"
      << "\"enough_time\":" << js(enough_time) << ','
      << "\"enough_cn_transitions\":" << js(enough_cn_transitions) << ','
      << "\"bidirectional_cn_motion\":" << js(bidirectional_cn_motion) << ','
      << "\"enough_cn_range\":" << js(enough_cn_range) << ','
      << "\"enough_rg_range\":" << js(enough_rg_range) << ','
      << "\"enough_contact_patterns\":" << js(enough_contact_patterns) << ','
      << "\"enough_cluster_diversity\":" << js(enough_cluster_diversity)
      << "},\"stop\":" << js(stop) << '}';
    return o.str();
}

SamplingMetrics evaluate_sampling(const std::vector<std::filesystem::path>& colvar_paths,
                                  int n_prot,
                                  double dt_colvar_ps,
                                  const SamplingRules& rules) {
    const auto table = read_colvars(colvar_paths);
    const auto cn = column_values(table, "cn_total");
    const auto rg = column_values(table, "rg_global");
    const auto [cn_min_it, cn_max_it] = std::minmax_element(cn.begin(), cn.end());
    const auto [rg_min_it, rg_max_it] = std::minmax_element(rg.begin(), rg.end());

    SamplingMetrics m;
    m.n_frames = static_cast<int>(table.rows.size());
    m.total_time_us = static_cast<double>(m.n_frames) * dt_colvar_ps / 1'000'000.0;
    m.cn_min = *cn_min_it;
    m.cn_max = *cn_max_it;
    m.cn_range = m.cn_max - m.cn_min;
    m.rg_global_min = *rg_min_it;
    m.rg_global_max = *rg_max_it;
    m.rg_global_range = m.rg_global_max - m.rg_global_min;

    const auto transitions = count_transitions(
        compress_states_with_residence(
            discretize_three_states(moving_average(cn, rules.smooth_window_frames)),
            rules.min_residence_frames));
    m.cn_state_transitions = transitions.transitions;
    m.cn_state_up = transitions.up;
    m.cn_state_down = transitions.down;

    std::vector<std::pair<int, int>> pairs;
    std::vector<std::size_t> cols;
    for (int i = 1; i <= n_prot; ++i) {
        for (int j = i + 1; j <= n_prot; ++j) {
            const auto it = table.col.find(pair_contact_name(i, j));
            if (it == table.col.end()) throw std::runtime_error("Missing expected pair contact column: " + pair_contact_name(i, j));
            pairs.emplace_back(i - 1, j - 1);
            cols.push_back(it->second);
        }
    }

    std::set<std::string> patterns;
    std::set<int> cluster_sizes;
    const auto step = static_cast<std::size_t>(std::max(1, rules.pattern_downsample));

    for (std::size_t frame = 0; frame < table.rows.size(); frame += step) {
        std::string pattern;
        pattern.reserve(cols.size());
        std::vector<std::vector<int>> adj(static_cast<std::size_t>(n_prot));

        for (std::size_t k = 0; k < cols.size(); ++k) {
            const bool contact = table.rows[frame][cols[k]] >= rules.pair_contact_threshold;
            pattern.push_back(contact ? '1' : '0');
            if (!contact) continue;
            const auto [a, b] = pairs[k];
            adj[static_cast<std::size_t>(a)].push_back(b);
            adj[static_cast<std::size_t>(b)].push_back(a);
        }

        patterns.insert(std::move(pattern));
        const auto lcc = largest_connected_component(n_prot, adj);
        cluster_sizes.insert(lcc);
        m.largest_cluster_max = std::max(m.largest_cluster_max, lcc);
    }

    m.unique_contact_patterns = static_cast<int>(patterns.size());
    m.largest_cluster_unique = static_cast<int>(cluster_sizes.size());
    m.enough_time = m.total_time_us >= rules.min_total_us;
    m.enough_cn_transitions = m.cn_state_transitions >= rules.min_cn_state_transitions;
    m.bidirectional_cn_motion = std::min(m.cn_state_up, m.cn_state_down) >= rules.min_bidirectional_events;
    m.enough_cn_range = m.cn_range >= rules.min_cn_range;
    m.enough_rg_range = m.rg_global_range >= rules.min_rg_global_range;
    m.enough_contact_patterns = m.unique_contact_patterns >= rules.min_unique_contact_patterns;
    m.enough_cluster_diversity = m.largest_cluster_unique >= rules.min_largest_cluster_unique;
    m.stop = m.enough_time && m.enough_cn_transitions && m.bidirectional_cn_motion &&
             m.enough_cn_range && m.enough_rg_range && m.enough_contact_patterns &&
             m.enough_cluster_diversity;
    return m;
}

void append_metrics_jsonl(const std::filesystem::path& path, const SamplingMetrics& metrics, int last_chunk) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::app);
    if (!out) throw std::runtime_error("Cannot write metrics JSONL: " + path.string());
    out << metrics.to_json(last_chunk) << '\n';
}

} // namespace cg
