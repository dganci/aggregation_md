#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace cg {


/// Causal moving average over the last `window` values.
std::vector<double> moving_average(const std::vector<double>& x, int window);

/// Low / middle / high thirds of the observed range, as 0 / 1 / 2.
std::vector<int> discretize_three_states(const std::vector<double>& x);

/// Drops every visit shorter than `min_residence` samples and collapses the
/// rest to the sequence of states actually dwelt in.
std::vector<int> compress_states_with_residence(const std::vector<int>& states, int min_residence);

struct TransitionCounts { int transitions = 0, up = 0, down = 0; };

/// How often a state sequence changed, and in which direction.
TransitionCounts count_transitions(const std::vector<int>& states);

/// Effective number of independent samples in a correlated series, from the
/// block-averaging plateau (Flyvbjerg & Petersen 1989)
double effective_sample_size(const std::vector<double>& x, int min_blocks = 5);

/// Jensen-Shannon divergence between two frequency tables, in bits: 0 when they
/// describe the same distribution, 1 when they share no key at all.
double jensen_shannon_bits(const std::unordered_map<std::string, int>& a,
                           const std::unordered_map<std::string, int>& b);

} // namespace cg
