#include "SamplingMonitor.hpp"
#include "ClusterAnalysis.hpp"
#include "ColvarTable.hpp"
#include "SamplingEstimators.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace cg {
namespace {

std::string js(bool x) { return x ? "true" : "false"; }
std::string pair_contact_name(int i, int j) { return "cn_" + std::to_string(i) + "_" + std::to_string(j); }

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
      << "\"patterns_first_half\":" << patterns_first_half << ','
      << "\"pattern_growth_ratio\":" << pattern_growth_ratio << ','
      << "\"pattern_jsd_halves_bits\":" << pattern_jsd_halves_bits << ','
      << "\"populated_contact_patterns\":" << populated_contact_patterns << ','
      << "\"largest_cluster_unique\":" << largest_cluster_unique << ','
      << "\"largest_cluster_max\":" << largest_cluster_max << ','
      << "\"largest_cluster_mean\":" << largest_cluster_mean << ','
      << "\"lcc_growth_events\":" << lcc_growth_events << ','
      << "\"lcc_shrink_events\":" << lcc_shrink_events << ','
      << "\"assembly_growth_events\":" << assembly_growth_events << ','
      << "\"assembly_shrink_events\":" << assembly_shrink_events << ','
      << "\"cn_effective_samples\":" << cn_effective_samples << ','
      << "\"cv_its1_ps\":" << cv_its1_ps << ','
      << "\"cv_plateau_found\":" << js(cv_plateau_found) << ','
      << "\"cv_time_over_its\":" << cv_time_over_its << ','
      << "\"duplicate_frames_dropped\":" << duplicate_frames_dropped << ','
      << "\"checks\":{"
      << "\"enough_time\":" << js(enough_time) << ','
      << "\"enough_cn_transitions\":" << js(enough_cn_transitions) << ','
      << "\"bidirectional_cn_motion\":" << js(bidirectional_cn_motion) << ','
      << "\"enough_cn_range\":" << js(enough_cn_range) << ','
      << "\"enough_rg_range\":" << js(enough_rg_range) << ','
      << "\"enough_contact_patterns\":" << js(enough_contact_patterns) << ','
      << "\"exploration_saturated\":" << js(exploration_saturated) << ','
      << "\"pattern_distribution_settled\":" << js(pattern_distribution_settled) << ','
      << "\"enough_independent_samples\":" << js(enough_independent_samples) << ','
      << "\"enough_assembly_events\":" << js(enough_assembly_events) << ','
      << "\"enough_for_cv_training\":" << js(enough_for_cv_training) << ','
      << "\"enough_cluster_diversity\":" << js(enough_cluster_diversity)
      << "},\"stop\":" << js(stop) << '}';
    return o.str();
}

std::vector<double> rg_com_from_distances(const ColvarTable& table, int n_prot) {
    std::vector<std::size_t> cols;
    for (int i = 1; i <= n_prot; ++i) {
        for (int j = i + 1; j <= n_prot; ++j) {
            const auto name = "d_" + std::to_string(i) + "_" + std::to_string(j);
            const auto it = table.col.find(name);
            if (it == table.col.end())
                throw std::runtime_error(
                    "COLVAR has neither an 'rg_com' column nor the '" + name + "' column needed to "
                    "derive it. Either it was written by a different --n-prot, or it predates the "
                    "descriptor set entirely.");
            cols.push_back(it->second);
        }
    }
    const double norm = 1.0 / (static_cast<double>(n_prot) * n_prot);
    std::vector<double> out;
    out.reserve(table.rows.size());
    for (const auto& row : table.rows) {
        double sq = 0.0;
        for (const auto c : cols) sq += row[c] * row[c];
        out.push_back(std::sqrt(sq * norm));
    }
    return out;
}

SamplingMetrics evaluate_sampling(const std::vector<std::filesystem::path>& colvar_paths,
                                  int n_prot,
                                  double dt_colvar_ps,
                                  const SamplingRules& rules) {
    const auto table = read_colvars(colvar_paths);
    const auto cn = column_values(table, "cn_total");
    const auto rg = table.col.count("rg_com")
                        ? column_values(table, "rg_com")
                        : rg_com_from_distances(table, n_prot);
    const auto [cn_min_it, cn_max_it] = std::minmax_element(cn.begin(), cn.end());
    const auto [rg_min_it, rg_max_it] = std::minmax_element(rg.begin(), rg.end());

    SamplingMetrics m;
    m.n_frames = static_cast<int>(table.rows.size());
    m.total_time_us = table.total_time_ps(dt_colvar_ps) / 1'000'000.0;
    m.duplicate_frames_dropped = static_cast<int>(table.dropped_non_monotonic_rows);
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

    std::unordered_map<std::string, int> patterns_first_half, patterns_second_half;
    const std::size_t half = table.rows.size() / 2;
    std::set<int> cluster_sizes;
    const auto step = static_cast<std::size_t>(std::max(1, rules.pattern_downsample));

    std::size_t sampled_frames = 0;
    double lcc_sum = 0.0;
    int previous_lcc = -1;
    std::vector<int> lcc_series;
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

        ++(frame < half ? patterns_first_half : patterns_second_half)[std::move(pattern)];
        const auto lcc = largest_connected_component(n_prot, adj);
        cluster_sizes.insert(lcc);
        lcc_series.push_back(lcc);
        m.largest_cluster_max = std::max(m.largest_cluster_max, lcc);

        if (previous_lcc >= 0) {
            if (lcc > previous_lcc) ++m.lcc_growth_events;
            else if (lcc < previous_lcc) ++m.lcc_shrink_events;
        }
        previous_lcc = lcc;
        lcc_sum += lcc;
        ++sampled_frames;
    }

    std::unordered_map<std::string, int> patterns = patterns_first_half;
    for (const auto& [pattern, count] : patterns_second_half) patterns[pattern] += count;
    m.unique_contact_patterns = static_cast<int>(patterns.size());
    m.patterns_first_half = static_cast<int>(patterns_first_half.size());
    m.pattern_growth_ratio =
        m.patterns_first_half > 0
            ? static_cast<double>(m.unique_contact_patterns) / m.patterns_first_half
            : std::numeric_limits<double>::infinity();
    m.pattern_jsd_halves_bits = jensen_shannon_bits(patterns_first_half, patterns_second_half);
    const auto populated_floor = static_cast<double>(sampled_frames) / 100.0;
    for (const auto& [pattern, count] : patterns)
        if (count >= populated_floor) ++m.populated_contact_patterns;
    m.largest_cluster_unique = static_cast<int>(cluster_sizes.size());
    m.largest_cluster_mean = sampled_frames ? lcc_sum / static_cast<double>(sampled_frames) : 0.0;

    const double sample_spacing_ps = (table.dt_ps > 0.0 ? table.dt_ps : dt_colvar_ps) * static_cast<double>(step);
    const int residence_samples = sample_spacing_ps > 0.0
        ? std::max(1, static_cast<int>(std::lround(rules.event_residence_ps / sample_spacing_ps)))
        : 1;
    const auto assembly = count_transitions(compress_states_with_residence(lcc_series, residence_samples));
    m.assembly_growth_events = assembly.up;
    m.assembly_shrink_events = assembly.down;

    m.cn_effective_samples = effective_sample_size(cn);
    decide_stop(m, rules);
    return m;
}

void decide_stop(SamplingMetrics& m, const SamplingRules& rules) {
    m.enough_time = m.total_time_us >= rules.min_total_us;
    m.enough_cn_transitions = m.cn_state_transitions >= rules.min_cn_state_transitions;
    m.bidirectional_cn_motion = std::min(m.cn_state_up, m.cn_state_down) >= rules.min_bidirectional_events;
    m.enough_cn_range = m.cn_range >= rules.min_cn_range;
    m.enough_rg_range = m.rg_global_range >= rules.min_rg_global_range;
    m.enough_contact_patterns = m.unique_contact_patterns >= rules.min_unique_contact_patterns;
    m.enough_cluster_diversity = m.largest_cluster_unique >= rules.min_largest_cluster_unique;
    m.exploration_saturated = m.pattern_growth_ratio <= rules.max_pattern_growth_ratio;
    m.pattern_distribution_settled = m.pattern_jsd_halves_bits <= rules.max_pattern_jsd_bits;
    m.enough_independent_samples = m.cn_effective_samples >= rules.min_effective_samples;
    m.enough_assembly_events =
        std::min(m.assembly_growth_events, m.assembly_shrink_events) >= rules.min_assembly_events;
    m.enough_for_cv_training = rules.min_time_over_its <= 0.0 ||
        (m.cv_plateau_found && m.cv_time_over_its >= rules.min_time_over_its);
    m.stop = m.enough_time && m.enough_cn_transitions && m.bidirectional_cn_motion &&
             m.enough_cn_range && m.enough_rg_range && m.enough_contact_patterns &&
             m.enough_cluster_diversity && m.exploration_saturated &&
             m.pattern_distribution_settled && m.enough_independent_samples &&
             m.enough_assembly_events && m.enough_for_cv_training;
}

void append_metrics_jsonl(const std::filesystem::path& path, const SamplingMetrics& metrics, int last_chunk) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::app);
    if (!out) throw std::runtime_error("Cannot write metrics JSONL: " + path.string());
    out << metrics.to_json(last_chunk) << '\n';
}

} // namespace cg
