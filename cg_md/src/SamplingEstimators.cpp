#include "SamplingEstimators.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>

namespace cg {

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

double effective_sample_size(const std::vector<double>& x, int min_blocks) {
    const auto n = x.size();
    min_blocks = std::max(min_blocks, 2);
    if (n < 2 * static_cast<std::size_t>(min_blocks)) return 0.0;

    double mean = 0.0;
    for (const double v : x) mean += v;
    mean /= static_cast<double>(n);
    double var = 0.0;
    for (const double v : x) var += (v - mean) * (v - mean);
    var /= static_cast<double>(n);
    if (var <= 0.0) return 0.0;

    double se2_max = 0.0;
    std::vector<double> means;
    for (std::size_t block = 1; n / block >= static_cast<std::size_t>(min_blocks); block *= 2) {
        const auto n_blocks = n / block;
        means.assign(n_blocks, 0.0);
        for (std::size_t k = 0; k < n_blocks; ++k) {
            double sum = 0.0;
            for (std::size_t i = k * block; i < (k + 1) * block; ++i) sum += x[i];
            means[k] = sum / static_cast<double>(block);
        }
        double m = 0.0;
        for (const double v : means) m += v;
        m /= static_cast<double>(n_blocks);
        double spread = 0.0;
        for (const double v : means) spread += (v - m) * (v - m);
        spread /= static_cast<double>(n_blocks - 1);
        se2_max = std::max(se2_max, spread / static_cast<double>(n_blocks));
    }
    return se2_max > 0.0 ? var / se2_max : 0.0;
}

double jensen_shannon_bits(const std::unordered_map<std::string, int>& a,
                           const std::unordered_map<std::string, int>& b) {
    double total_a = 0.0, total_b = 0.0;
    for (const auto& [key, count] : a) total_a += count;
    for (const auto& [key, count] : b) total_b += count;
    if (total_a <= 0.0 || total_b <= 0.0) return std::numeric_limits<double>::infinity();

    const auto share = [](const std::unordered_map<std::string, int>& table,
                          const std::string& key, double total) {
        const auto it = table.find(key);
        return it == table.end() ? 0.0 : it->second / total;
    };
    const auto term = [](double p, double m) { return p > 0.0 ? p * std::log2(p / m) : 0.0; };

    double jsd = 0.0;
    for (const auto* table : {&a, &b}) {
        for (const auto& [key, count] : *table) {
            if (table == &b && a.count(key)) continue;
            const double p = share(a, key, total_a);
            const double q = share(b, key, total_b);
            const double m = 0.5 * (p + q);
            jsd += 0.5 * term(p, m) + 0.5 * term(q, m);
        }
    }
    return jsd;
}

} // namespace cg
