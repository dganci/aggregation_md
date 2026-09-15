#include "test_framework.hpp"
#include "ColvarTable.hpp"
#include "FileUtils.hpp"
#include "RunLedger.hpp"

#include <exception>
#include <filesystem>

using namespace cg;

namespace {
std::filesystem::path temp_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_ledger_" + name);
}
} // namespace

CG_TEST(a_fresh_ledger_records_whatever_the_first_invocation_used) {
    const auto path = temp_path("fresh.ledger");
    std::filesystem::remove(path);

    auto ledger = RunLedger::load(path);
    CG_CHECK(ledger.empty());
    ledger.require("chunk_nsteps", "50000", "--adaptive-chunk-us");
    ledger.save(path);

    auto reloaded = RunLedger::load(path);
    CG_CHECK(!reloaded.empty());
    reloaded.require("chunk_nsteps", "50000", "--adaptive-chunk-us");

    std::filesystem::remove(path);
}

CG_TEST(the_ledger_refuses_a_changed_segment_length_mid_run) {
    const auto path = temp_path("drift.ledger");
    std::filesystem::remove(path);

    auto first = RunLedger::load(path);
    first.require("chunk_nsteps", "50000", "--adaptive-chunk-us");
    first.save(path);

    auto second = RunLedger::load(path);
    bool threw = false;
    try {
        second.require("chunk_nsteps", "100000", "--adaptive-chunk-us");
    } catch (const std::exception&) {
        threw = true;
    }
    CG_CHECK(threw);

    std::filesystem::remove(path);
}

CG_TEST(truncating_a_plumed_file_drops_the_interrupted_attempts_segment) {
    const auto path = temp_path("hills.dat");
    write_lines(path, {
        "#! FIELDS time cv sigma height biasf",
        "#! SET multivariate false",
        "10.0 1.0 0.1 0.8 10",
        "20.0 1.1 0.1 0.7 10",
        "10.0 1.2 0.1 0.6 10",
        "20.0 1.3 0.1 0.5 10"
    });

    const auto removed = truncate_plumed_file_to_segment(path, 1);
    CG_CHECK_EQ(removed, static_cast<std::size_t>(2));

    const auto table = read_colvars({path});
    CG_CHECK_EQ(static_cast<int>(table.rows.size()), 2);
    CG_CHECK(table.t_last_ps == 20.0);
    const auto text = read_text(path);
    CG_CHECK(text.find("#! FIELDS") != std::string::npos);
    CG_CHECK(text.find("#! SET multivariate") != std::string::npos);

    CG_CHECK_EQ(truncate_plumed_file_to_segment(path, 1), static_cast<std::size_t>(0));

    std::filesystem::remove(path);
}

CG_TEST(a_segment_rewind_keeps_every_completed_batch) {
    const auto path = temp_path("hills3.dat");
    write_lines(path, {
        "#! FIELDS time cv sigma height biasf",
        "5.0 1.0 0.1 0.8 10", "10.0 1.1 0.1 0.8 10",
        "5.0 1.2 0.1 0.7 10", "10.0 1.3 0.1 0.7 10",
        "5.0 1.4 0.1 0.6 10"
    });

    CG_CHECK_EQ(truncate_plumed_file_to_segment(path, 2), static_cast<std::size_t>(1));
    CG_CHECK_EQ(count_plumed_records(path), static_cast<std::size_t>(4));

    std::filesystem::remove(path);
}

CG_TEST(truncating_a_missing_file_is_a_no_op) {
    CG_CHECK_EQ(truncate_plumed_file_to_segment(temp_path("does_not_exist.dat"), 1),
                static_cast<std::size_t>(0));
    CG_CHECK_EQ(truncate_plumed_file_to(temp_path("does_not_exist.dat"), 5),
                static_cast<std::size_t>(0));
    CG_CHECK_EQ(count_plumed_records(temp_path("does_not_exist.dat")),
                static_cast<std::size_t>(0));
}

CG_TEST(rewinding_by_record_count_is_exact_where_a_time_threshold_is_not) {
    const auto path = temp_path("hills_count.dat");
    write_lines(path, {
        "#! FIELDS time cv sigma height biasf",
        "0.0 1.0 0.1 0.8 10",
        "5.0 1.1 0.1 0.7 10",
        "10.0 1.2 0.1 0.6 10",
        "15.0 1.3 0.1 0.5 10"
    });
    CG_CHECK_EQ(count_plumed_records(path), static_cast<std::size_t>(4));

    CG_CHECK_EQ(truncate_plumed_file_to(path, 2), static_cast<std::size_t>(2));
    CG_CHECK_EQ(count_plumed_records(path), static_cast<std::size_t>(2));

    const auto table = read_colvars({path});
    CG_CHECK_EQ(static_cast<int>(table.rows.size()), 2);
    CG_CHECK(table.t_last_ps == 5.0);
    CG_CHECK(read_text(path).find("#! FIELDS") != std::string::npos);

    CG_CHECK_EQ(truncate_plumed_file_to(path, 2), static_cast<std::size_t>(0));
    CG_CHECK_EQ(truncate_plumed_file_to(path, 99), static_cast<std::size_t>(0));

    std::filesystem::remove(path);
}

CG_TEST(the_ledger_can_carry_bookkeeping_that_legitimately_changes) {
    const auto path = temp_path("records.ledger");
    std::filesystem::remove(path);

    auto ledger = RunLedger::load(path);
    ledger.set("after_000_HILLS", "12");
    ledger.save(path);

    auto reloaded = RunLedger::load(path);
    CG_CHECK(reloaded.get("after_000_HILLS") == "12");
    CG_CHECK(reloaded.get("after_001_HILLS").empty());
    reloaded.set("after_000_HILLS", "24");
    CG_CHECK(reloaded.get("after_000_HILLS") == "24");

    std::filesystem::remove(path);
}
