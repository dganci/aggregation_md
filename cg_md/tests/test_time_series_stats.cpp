#include "test_framework.hpp"
#include "TimeSeriesStats.hpp"

#include <cmath>

using namespace cg;

namespace {
bool close(double a, double b, double eps = 1e-9) { return std::abs(a - b) < eps; }
} // namespace

CG_TEST(compute_time_series_stats_on_empty_series_is_all_zero) {
    const auto stats = compute_time_series_stats({});
    CG_CHECK_EQ(stats.mean, 0.0);
    CG_CHECK_EQ(stats.stddev, 0.0);
    CG_CHECK_EQ(stats.drift_per_ns, 0.0);
}

CG_TEST(compute_time_series_stats_single_point_has_no_drift_or_spread) {
    const auto stats = compute_time_series_stats({{0.0, 42.0}});
    CG_CHECK_EQ(stats.mean, 42.0);
    CG_CHECK_EQ(stats.stddev, 0.0);
    CG_CHECK_EQ(stats.drift_per_ns, 0.0);
}

CG_TEST(compute_time_series_stats_constant_series_has_zero_drift_and_spread) {
    const auto stats = compute_time_series_stats({{0.0, 300.0}, {1000.0, 300.0}, {2000.0, 300.0}});
    CG_CHECK_EQ(stats.mean, 300.0);
    CG_CHECK_EQ(stats.stddev, 0.0);
    CG_CHECK_EQ(stats.drift_per_ns, 0.0);
}

CG_TEST(compute_time_series_stats_linear_trend_recovers_known_slope_and_mean) {
    const auto stats = compute_time_series_stats({{0.0, 1.0}, {1000.0, 2.0}, {2000.0, 3.0}});
    CG_CHECK(close(stats.mean, 2.0));
    CG_CHECK(close(stats.drift_per_ns, 1.0));
    CG_CHECK(stats.stddev > 0.0);
}

CG_TEST(to_json_contains_all_three_fields) {
    TimeSeriesStats stats;
    stats.mean = 1.0;
    stats.stddev = 2.0;
    stats.drift_per_ns = 3.0;
    const auto json = stats.to_json();
    CG_CHECK(json.find("\"mean\":1") != std::string::npos);
    CG_CHECK(json.find("\"stddev\":2") != std::string::npos);
    CG_CHECK(json.find("\"drift_per_ns\":3") != std::string::npos);
}
