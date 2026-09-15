#include "ClusterAnalysis.hpp"

#include <algorithm>
#include <queue>

namespace cg {
namespace {

std::vector<int> component_sizes_bfs(int n, const ContactGraph& graph) {
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    std::vector<int> sizes;

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
            for (int v : graph[static_cast<std::size_t>(u)]) {
                if (!seen[static_cast<std::size_t>(v)]) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    q.push(v);
                }
            }
        }
        sizes.push_back(size);
    }
    return sizes;
}

} // namespace

int largest_connected_component(int n, const ContactGraph& graph) {
    if (n <= 0) return 0;
    const auto sizes = component_sizes_bfs(n, graph);
    return *std::max_element(sizes.begin(), sizes.end());
}

std::vector<int> connected_component_sizes(int n, const ContactGraph& graph) {
    if (n <= 0) return {};
    auto sizes = component_sizes_bfs(n, graph);
    std::sort(sizes.begin(), sizes.end(), std::greater<int>());
    return sizes;
}

} // namespace cg
