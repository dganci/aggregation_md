#include "test_framework.hpp"
#include "ColvarTable.hpp"
#include "FileUtils.hpp"

#include <cmath>
#include <exception>
#include <filesystem>

using namespace cg;

namespace {
std::filesystem::path temp_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_colvar_" + name);
}
bool close(double a, double b, double eps = 1e-9) { return std::abs(a - b) < eps; }
} // namespace

CG_TEST(read_colvars_drops_the_duplicated_frame_at_a_batch_boundary) {
    const auto a = temp_path("batch_a.dat");
    const auto b = temp_path("batch_b.dat");

    write_lines(a, {"#! FIELDS time cn_total rg_com",
                    "0.0 1.0 2.0",
                    "10.0 2.0 2.1",
                    "20.0 3.0 2.2"});
    write_lines(b, {"#! FIELDS time cn_total rg_com",
                    "20.0 3.0 2.2",
                    "30.0 4.0 2.3",
                    "40.0 5.0 2.4"});

    const auto table = read_colvars({a, b});

    CG_CHECK(table.has_time);
    CG_CHECK_EQ(static_cast<int>(table.rows.size()), 5);
    CG_CHECK_EQ(static_cast<int>(table.dropped_non_monotonic_rows), 1);
    CG_CHECK(close(table.t_first_ps, 0.0));
    CG_CHECK(close(table.t_last_ps, 40.0));
    CG_CHECK(close(table.dt_ps, 10.0));
    CG_CHECK(close(table.total_time_ps(999.0), 50.0));

    const auto cn = column_values(table, "cn_total");
    CG_CHECK_EQ(static_cast<int>(cn.size()), 5);
    CG_CHECK(close(cn[2], 3.0));
    CG_CHECK(close(cn[3], 4.0));

    std::filesystem::remove(a);
    std::filesystem::remove(b);
}

CG_TEST(read_colvars_falls_back_to_frame_counting_without_a_time_column) {
    const auto path = temp_path("no_time.dat");
    write_lines(path, {"#! FIELDS cn_total rg_com", "1.0 2.0", "2.0 2.1"});

    const auto table = read_colvars({path});
    CG_CHECK(!table.has_time);
    CG_CHECK_EQ(static_cast<int>(table.rows.size()), 2);
    CG_CHECK(close(table.total_time_ps(250.0), 500.0));

    std::filesystem::remove(path);
}

CG_TEST(column_series_uses_the_real_time_axis_when_available) {
    const auto path = temp_path("offset.dat");
    write_lines(path, {"#! FIELDS time rg1", "100.0 1.0", "110.0 1.5", "120.0 2.0"});

    const auto table = read_colvars({path});
    const auto series = column_series(table, "rg1", 1.0);

    CG_CHECK_EQ(static_cast<int>(series.size()), 3);
    CG_CHECK(close(series.front().first, 100.0));
    CG_CHECK(close(series.back().first, 120.0));
    CG_CHECK(close(series.back().second, 2.0));

    std::filesystem::remove(path);
}

CG_TEST(read_colvars_rejects_a_header_mismatch_between_batches) {
    const auto a = temp_path("hdr_a.dat");
    const auto b = temp_path("hdr_b.dat");
    write_lines(a, {"#! FIELDS time cn_total", "0.0 1.0"});
    write_lines(b, {"#! FIELDS time rg_com", "10.0 1.0"});

    bool threw = false;
    try {
        read_colvars({a, b});
    } catch (const std::exception&) {
        threw = true;
    }
    CG_CHECK(threw);

    std::filesystem::remove(a);
    std::filesystem::remove(b);
}


namespace {

std::filesystem::path write_colvar(const std::string& name, const std::string& body) {
    const auto p = std::filesystem::temp_directory_path() / ("cg_md_test_" + name + ".dat");
    write_text(p, "#! FIELDS time cn_total\n" + body);
    return p;
}

} // namespace

CG_TEST(concatenated_chunks_with_restarted_plumed_clocks_accumulate_time) {
    const auto a = write_colvar("clk_a", "0 1\n10 2\n20 3\n30 4\n");
    const auto b = write_colvar("clk_b", "0 5\n10 6\n20 7\n30 8\n");

    const auto t = read_colvars({a, b});

    CG_CHECK_EQ(t.rows.size(), std::size_t{8});
    CG_CHECK_EQ(t.dropped_non_monotonic_rows, std::size_t{0});

    CG_CHECK_EQ(t.rows[4][0], 40.0);
    CG_CHECK_EQ(t.rows[7][0], 70.0);

    CG_CHECK_EQ(t.rows[4][1], 5.0);

    CG_CHECK_EQ(t.total_time_ps(10.0), 80.0);

    std::filesystem::remove(a);
    std::filesystem::remove(b);
}

CG_TEST(absolute_time_files_are_still_concatenated_unshifted) {
    const auto a = write_colvar("abs_a", "0 1\n10 2\n20 3\n");
    const auto b = write_colvar("abs_b", "20 3\n30 4\n40 5\n");

    const auto t = read_colvars({a, b});

    CG_CHECK_EQ(t.rows.size(), std::size_t{5});
    CG_CHECK_EQ(t.dropped_non_monotonic_rows, std::size_t{1});
    CG_CHECK_EQ(t.rows[0][0], 0.0);
    CG_CHECK_EQ(t.rows[4][0], 40.0);

    std::filesystem::remove(a);
    std::filesystem::remove(b);
}

CG_TEST(three_restarted_chunks_chain_correctly) {
    const auto a = write_colvar("clk3_a", "0 1\n5 2\n");
    const auto b = write_colvar("clk3_b", "0 3\n5 4\n");
    const auto c = write_colvar("clk3_c", "0 5\n5 6\n");

    const auto t = read_colvars({a, b, c});

    CG_CHECK_EQ(t.rows.size(), std::size_t{6});
    CG_CHECK_EQ(t.dropped_non_monotonic_rows, std::size_t{0});
    CG_CHECK_EQ(t.rows[2][0], 10.0);
    CG_CHECK_EQ(t.rows[4][0], 20.0);
    CG_CHECK_EQ(t.rows[5][0], 25.0);

    std::filesystem::remove(a);
    std::filesystem::remove(b);
    std::filesystem::remove(c);
}

CG_TEST(a_clock_restart_inside_one_file_is_rebased_too) {
    const auto path = write_colvar("inline_restart", "0 1\n10 2\n20 3\n0 4\n10 5\n20 6\n");

    const auto t = read_colvars({path});

    CG_CHECK_EQ(t.rows.size(), std::size_t{6});
    CG_CHECK_EQ(t.dropped_non_monotonic_rows, std::size_t{0});
    CG_CHECK_EQ(t.rows[3][0], 30.0);
    CG_CHECK_EQ(t.rows[5][0], 50.0);
    CG_CHECK_EQ(t.rows[3][1], 4.0);

    std::filesystem::remove(path);
}

CG_TEST(a_malformed_number_names_the_file_it_came_from) {
    const auto path = write_colvar("malformed", "0 1\n10 not-a-number\n");

    std::string message;
    try {
        read_colvars({path});
    } catch (const std::exception& e) {
        message = e.what();
    }
    CG_CHECK(message.find("not-a-number") != std::string::npos);
    CG_CHECK(message.find("malformed") != std::string::npos);

    std::filesystem::remove(path);
}
