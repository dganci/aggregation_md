#include "test_framework.hpp"
#include "ClusterAnalysis.hpp"

using namespace cg;

CG_TEST(largest_connected_component_finds_max_size) {
    ContactGraph graph(5);
    graph[0] = {1};
    graph[1] = {0, 2};
    graph[2] = {1};
    CG_CHECK_EQ(largest_connected_component(5, graph), 3);
}

CG_TEST(largest_connected_component_all_isolated_is_one) {
    ContactGraph graph(4);
    CG_CHECK_EQ(largest_connected_component(4, graph), 1);
}

CG_TEST(largest_connected_component_empty_graph_is_zero) {
    ContactGraph graph;
    CG_CHECK_EQ(largest_connected_component(0, graph), 0);
}

CG_TEST(connected_component_sizes_returns_all_sizes_descending) {
    ContactGraph graph(5);
    graph[0] = {1};
    graph[1] = {0, 2};
    graph[2] = {1};
    graph[3] = {4};
    graph[4] = {3};

    CG_CHECK_EQ(connected_component_sizes(5, graph), (std::vector<int>{3, 2}));
}

CG_TEST(connected_component_sizes_fully_isolated_returns_all_ones) {
    ContactGraph graph(3);
    CG_CHECK_EQ(connected_component_sizes(3, graph), (std::vector<int>{1, 1, 1}));
}

CG_TEST(connected_component_sizes_fully_connected_returns_single_component) {
    ContactGraph graph(3);
    graph[0] = {1, 2};
    graph[1] = {0, 2};
    graph[2] = {0, 1};
    CG_CHECK_EQ(connected_component_sizes(3, graph), (std::vector<int>{3}));
}
