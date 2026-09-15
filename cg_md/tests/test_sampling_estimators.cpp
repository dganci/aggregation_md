#include "test_framework.hpp"
#include "SamplingEstimators.hpp"

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cg;


CG_TEST(effective_sample_size_is_the_frame_count_for_uncorrelated_data) {
    std::vector<double> x(1024);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = static_cast<double>(i % 2);
    const auto n_eff = effective_sample_size(x);
    CG_CHECK(n_eff > 1000.0 && n_eff < 1030.0);
}

CG_TEST(effective_sample_size_counts_the_independent_blocks_of_a_correlated_series) {
    std::vector<double> x(1024);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = static_cast<double>((i / 128) % 2);
    const auto n_eff = effective_sample_size(x);
    CG_CHECK(n_eff > 6.5 && n_eff < 7.5);
}

CG_TEST(effective_sample_size_is_zero_when_nothing_moves) {
    CG_CHECK_EQ(effective_sample_size(std::vector<double>(200, 3.0)), 0.0);
    CG_CHECK_EQ(effective_sample_size({1.0, 2.0, 3.0}), 0.0);
}


CG_TEST(jensen_shannon_is_zero_for_the_same_frequencies_and_one_bit_when_disjoint) {
    const std::unordered_map<std::string, int> a{{"100", 10}, {"010", 30}};
    const std::unordered_map<std::string, int> same_shape{{"100", 1}, {"010", 3}};
    const std::unordered_map<std::string, int> other{{"001", 5}};
    CG_CHECK(jensen_shannon_bits(a, a) < 1e-12);
    CG_CHECK(jensen_shannon_bits(a, same_shape) < 1e-12);
    CG_CHECK(std::abs(jensen_shannon_bits(a, other) - 1.0) < 1e-12);
}

CG_TEST(jensen_shannon_matches_a_hand_calculation) {
    const std::unordered_map<std::string, int> p{{"a", 4}};
    const std::unordered_map<std::string, int> q{{"a", 2}, {"b", 2}};
    CG_CHECK(std::abs(jensen_shannon_bits(p, q) - 0.31128) < 1e-4);
    CG_CHECK(std::abs(jensen_shannon_bits(q, p) - 0.31128) < 1e-4);
}

CG_TEST(jensen_shannon_of_an_empty_half_is_not_zero) {
    const std::unordered_map<std::string, int> a{{"100", 10}};
    CG_CHECK(std::isinf(jensen_shannon_bits(a, {})));
    CG_CHECK(std::isinf(jensen_shannon_bits({}, a)));
}
