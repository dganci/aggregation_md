#include "test_framework.hpp"
#include "SamplingMonitor.hpp"
#include <cmath>
#include "FileUtils.hpp"

#include <exception>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cg;

namespace {
std::filesystem::path temp_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_" + name);
}
} // namespace

CG_TEST(evaluate_sampling_computes_expected_metrics_across_multiple_colvar_files) {
    const auto path_a = temp_path("colvar_a.dat");
    const auto path_b = temp_path("colvar_b.dat");

    write_lines(path_a, {
        "#! FIELDS cn_1_2 rg_com cn_total",
        "0.0 1.0 0.0",
        "0.0 1.0 0.0"
    });
    write_lines(path_b, {
        "#! FIELDS cn_1_2 rg_com cn_total",
        "30.0 1.5 30.0",
        "30.0 1.5 30.0"
    });

    SamplingRules rules;
    rules.smooth_window_frames = 1;
    rules.min_residence_frames = 1;
    rules.pair_contact_threshold = 1.0;
    rules.pattern_downsample = 1;
    rules.min_bidirectional_events = 1;

    const auto m = evaluate_sampling({path_a, path_b}, /*n_prot=*/2, /*dt_colvar_ps=*/250000.0, rules);

    CG_CHECK_EQ(m.n_frames, 4);
    CG_CHECK_EQ(m.total_time_us, 1.0);
    CG_CHECK_EQ(m.cn_min, 0.0);
    CG_CHECK_EQ(m.cn_max, 30.0);
    CG_CHECK_EQ(m.cn_range, 30.0);
    CG_CHECK_EQ(m.rg_global_min, 1.0);
    CG_CHECK_EQ(m.rg_global_max, 1.5);
    CG_CHECK_EQ(m.rg_global_range, 0.5);
    CG_CHECK_EQ(m.cn_state_transitions, 1);
    CG_CHECK_EQ(m.cn_state_up, 1);
    CG_CHECK_EQ(m.cn_state_down, 0);
    CG_CHECK(!m.bidirectional_cn_motion);
    CG_CHECK_EQ(m.unique_contact_patterns, 2);
    CG_CHECK_EQ(m.largest_cluster_unique, 2);
    CG_CHECK_EQ(m.largest_cluster_max, 2);

    const auto json = m.to_json(0);
    CG_CHECK(json.find("\"last_chunk\":0") != std::string::npos);
    CG_CHECK(json.find("\"n_frames\":4") != std::string::npos);

    std::filesystem::remove(path_a);
    std::filesystem::remove(path_b);
}

CG_TEST(evaluate_sampling_throws_on_missing_pair_contact_column) {
    const auto path = temp_path("colvar_missing_col.dat");
    write_lines(path, {
        "#! FIELDS rg_com cn_total",
        "1.0 0.0"
    });

    SamplingRules rules;
    bool threw = false;
    try {
        evaluate_sampling({path}, /*n_prot=*/2, 1000.0, rules);
    } catch (const std::exception&) {
        threw = true;
    }
    CG_CHECK(threw);

    std::filesystem::remove(path);
}

CG_TEST(exploration_saturation_distinguishes_a_still_exploring_run) {
    const auto dir = std::filesystem::temp_directory_path() / "cg_md_test_sat";
    std::filesystem::create_directories(dir);

    const std::string header = "#! FIELDS time cn_1_2 cn_1_3 cn_2_3 rg_com cn_total";
    const std::string only_12  = " 2.0 0.0 0.0 1.0 2.0";
    const std::string only_23  = " 0.0 0.0 2.0 1.0 2.0";

    SamplingRules rules;
    rules.pattern_downsample = 1;
    rules.pair_contact_threshold = 1.0;

    const auto exploring = dir / "exploring.dat";
    {
        std::vector<std::string> lines{header};
        for (int i = 0; i < 40; ++i) lines.push_back(std::to_string(i) + only_12);
        for (int i = 40; i < 80; ++i) lines.push_back(std::to_string(i) + only_23);
        write_lines(exploring, lines);
    }
    const auto m1 = evaluate_sampling({exploring}, 3, 1.0, rules);
    CG_CHECK_EQ(m1.patterns_first_half, 1);
    CG_CHECK_EQ(m1.unique_contact_patterns, 2);
    CG_CHECK(m1.pattern_growth_ratio > 1.0);
    CG_CHECK(!m1.exploration_saturated);
    CG_CHECK(std::abs(m1.pattern_jsd_halves_bits - 1.0) < 1e-9);
    CG_CHECK(!m1.pattern_distribution_settled);

    const auto settled = dir / "settled.dat";
    {
        std::vector<std::string> lines{header};
        for (int i = 0; i < 80; ++i)
            lines.push_back(std::to_string(i) + (i % 2 ? only_12 : only_23));
        write_lines(settled, lines);
    }
    const auto m2 = evaluate_sampling({settled}, 3, 1.0, rules);
    CG_CHECK_EQ(m2.patterns_first_half, 2);
    CG_CHECK_EQ(m2.unique_contact_patterns, 2);
    CG_CHECK(std::abs(m2.pattern_growth_ratio - 1.0) < 1e-9);
    CG_CHECK(m2.exploration_saturated);
    CG_CHECK(m2.pattern_jsd_halves_bits < 1e-9);
    CG_CHECK(m2.pattern_distribution_settled);

    std::filesystem::remove_all(dir);
}

CG_TEST(rg_com_is_derived_from_distances_when_the_column_is_absent) {
    const auto dir = std::filesystem::temp_directory_path() / "cg_md_test_rgderive";
    std::filesystem::create_directories(dir);
    const auto old_style = dir / "old.dat";
    write_lines(old_style, {"#! FIELDS time d_1_2 cn_1_2 cn_total",
                            "0 4.0 2.0 2.0",
                            "1 8.0 2.0 2.0"});
    SamplingRules rules;
    rules.pattern_downsample = 1;
    const auto m = evaluate_sampling({old_style}, 2, 1.0, rules);
    CG_CHECK(std::abs(m.rg_global_min - 2.0) < 1e-9);
    CG_CHECK(std::abs(m.rg_global_max - 4.0) < 1e-9);
    std::filesystem::remove_all(dir);
}


CG_TEST(assembly_events_ignore_a_grazing_contact_and_count_a_bound_one) {
    const auto dir = std::filesystem::temp_directory_path() / "cg_md_test_events";
    std::filesystem::create_directories(dir);
    const auto path = dir / "events.dat";
    std::vector<std::string> lines{"#! FIELDS time cn_1_2 rg_com cn_total"};
    int t = 0;
    const auto add = [&](int frames, double cn) {
        for (int i = 0; i < frames; ++i, t += 10)
            lines.push_back(std::to_string(t) + " " + std::to_string(cn) + " 1.0 " + std::to_string(cn));
    };
    add(20, 0.0); add(2, 5.0); add(20, 0.0); add(10, 5.0); add(20, 0.0);
    write_lines(path, lines);

    SamplingRules rules;
    rules.pattern_downsample = 1;
    rules.event_residence_ps = 30.0;
    rules.min_assembly_events = 1;
    const auto m = evaluate_sampling({path}, 2, 10.0, rules);
    CG_CHECK_EQ(m.lcc_growth_events, 2);
    CG_CHECK_EQ(m.lcc_shrink_events, 2);
    CG_CHECK_EQ(m.assembly_growth_events, 1);
    CG_CHECK_EQ(m.assembly_shrink_events, 1);
    CG_CHECK(m.enough_assembly_events);

    rules.event_residence_ps = 500.0;
    const auto none = evaluate_sampling({path}, 2, 10.0, rules);
    CG_CHECK_EQ(none.assembly_growth_events, 0);
    CG_CHECK(!none.enough_assembly_events);

    std::filesystem::remove_all(dir);
}


CG_TEST(a_run_that_only_ever_assembled_does_not_stop) {
    const auto dir = std::filesystem::temp_directory_path() / "cg_md_test_oneway";
    std::filesystem::create_directories(dir);
    const auto path = dir / "oneway.dat";
    std::vector<std::string> lines{"#! FIELDS time cn_1_2 rg_com cn_total"};
    for (int i = 0; i < 100; ++i) lines.push_back(std::to_string(i * 10) + " 0.0 3.0 0.0");
    for (int i = 100; i < 200; ++i) lines.push_back(std::to_string(i * 10) + " 5.0 1.0 5.0");
    write_lines(path, lines);

    SamplingRules rules;
    rules.min_total_us = 0.0;
    rules.min_cn_state_transitions = 0;
    rules.min_bidirectional_events = 0;
    rules.min_cn_range = 0.0;
    rules.min_rg_global_range = 0.0;
    rules.min_unique_contact_patterns = 0;
    rules.min_largest_cluster_unique = 0;
    rules.max_pattern_growth_ratio = 10.0;
    rules.max_pattern_jsd_bits = 1.0;
    rules.min_effective_samples = 0;
    rules.min_assembly_events = 1;
    rules.min_time_over_its = 0.0;
    rules.pattern_downsample = 1;
    const auto m = evaluate_sampling({path}, 2, 10.0, rules);
    CG_CHECK_EQ(m.assembly_growth_events, 1);
    CG_CHECK_EQ(m.assembly_shrink_events, 0);
    CG_CHECK(!m.enough_assembly_events);
    CG_CHECK(!m.stop);

    std::filesystem::remove_all(dir);
}
