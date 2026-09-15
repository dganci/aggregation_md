#pragma once

#include <vector>

namespace cg {

/// Undirected adjacency list over `n` nodes (protomers), built from whichever
/// pairwise contact definition the caller supplies (e.g. COLVAR cn_i_j
/// exceeding a threshold).
using ContactGraph = std::vector<std::vector<int>>;

/// Size of the largest connected component in `graph` (0 if n == 0).
int largest_connected_component(int n, const ContactGraph& graph);

/// Sizes of every connected component in `graph`, in descending order - e.g. 5
/// protomers split into a trimer and a dimer returns {3, 2}.
std::vector<int> connected_component_sizes(int n, const ContactGraph& graph);

} // namespace cg
