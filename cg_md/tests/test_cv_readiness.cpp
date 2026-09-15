#include "test_framework.hpp"
#include "CvReadiness.hpp"
#include "FileUtils.hpp"
#include "SamplingMonitor.hpp"

#include <cmath>
#include <filesystem>
#include <string>

using namespace cg;

namespace {
std::filesystem::path temp_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_" + name);
}
} // namespace


CG_TEST(cv_readiness_is_read_back_and_gates_on_plateau_and_ratio) {
    const auto path = temp_path("cv_readiness.txt");
    write_lines(path, {"plateau 1", "its1_ps 120000", "lag_lo_ps 100", "lag_hi_ps 1000",
                       "time_over_its 1.7", "frames 100000", "dt_ps 2", "total_time_ps 206000"});
    const auto cv = read_cv_readiness(path);
    CG_CHECK(cv.ok);
    CG_CHECK(cv.plateau);
    CG_CHECK(std::abs(cv.its1_ps - 120000.0) < 1e-6);

    SamplingRules rules;
    rules.min_time_over_its = 10.0;
    SamplingMetrics m;
    m.total_time_us = 0.206;
    apply_cv_readiness(m, rules, cv);
    CG_CHECK(m.cv_plateau_found);
    CG_CHECK(std::abs(m.cv_time_over_its - 0.206e6 / 120000.0) < 1e-9);
    CG_CHECK(!m.enough_for_cv_training);

    m.total_time_us = 1.5;
    apply_cv_readiness(m, rules, cv);
    CG_CHECK(m.enough_for_cv_training);

    std::filesystem::remove(path);
}

CG_TEST(cv_readiness_without_a_plateau_never_passes_whatever_the_ratio) {
    const auto path = temp_path("cv_readiness_noplateau.txt");
    write_lines(path, {"plateau 0", "its1_ps 1000", "total_time_ps 5000000"});
    SamplingRules rules;
    SamplingMetrics m;
    m.total_time_us = 5.0;
    apply_cv_readiness(m, rules, read_cv_readiness(path));
    CG_CHECK(m.cv_time_over_its > 1000.0);
    CG_CHECK(!m.enough_for_cv_training);
    std::filesystem::remove(path);
}

CG_TEST(a_missing_or_broken_readiness_file_leaves_the_gate_closed) {
    SamplingRules rules;
    SamplingMetrics m;
    m.total_time_us = 5.0;
    apply_cv_readiness(m, rules, read_cv_readiness(temp_path("does_not_exist.txt")));
    CG_CHECK(!m.cv_plateau_found);
    CG_CHECK_EQ(m.cv_time_over_its, 0.0);
    CG_CHECK(!m.enough_for_cv_training);

    const auto broken = temp_path("cv_readiness_broken.txt");
    write_lines(broken, {"plateau 1", "its1_ps twelve"});
    const auto cv = read_cv_readiness(broken);
    CG_CHECK(!cv.ok);
    std::filesystem::remove(broken);

    const auto noted = temp_path("cv_readiness_note.txt");
    write_lines(noted, {"note cvgen_estimators.py not found; set CG_CVGEN_SCRIPTS"});
    CG_CHECK(!read_cv_readiness(noted).ok);
    std::filesystem::remove(noted);
}

CG_TEST(the_cv_gate_can_be_switched_off) {
    SamplingRules rules;
    rules.min_time_over_its = 0.0;
    SamplingMetrics m;
    decide_stop(m, rules);
    CG_CHECK(m.enough_for_cv_training);
    rules.min_time_over_its = 10.0;
    decide_stop(m, rules);
    CG_CHECK(!m.enough_for_cv_training);
}
