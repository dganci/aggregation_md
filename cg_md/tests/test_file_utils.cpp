#include "test_framework.hpp"
#include "FileUtils.hpp"

#include <filesystem>
#include <string>

using namespace cg;

namespace {

std::filesystem::path scratch(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_" + name);
}

} // namespace

CG_TEST(read_tail_returns_the_whole_file_when_it_is_smaller_than_the_window) {
    const auto p = scratch("tail_small.txt");
    write_text(p, "short\n");
    CG_CHECK_EQ(read_tail(p, 4096), std::string("short\n"));
    std::filesystem::remove(p);
}

CG_TEST(read_tail_returns_only_the_last_bytes_of_a_larger_file) {
    const auto p = scratch("tail_large.txt");
    write_text(p, std::string(10000, 'a') + "TAIL");
    const auto tail = read_tail(p, 8);
    CG_CHECK_EQ(tail.size(), std::size_t{8});
    CG_CHECK_EQ(tail.substr(4), std::string("TAIL"));
    std::filesystem::remove(p);
}

CG_TEST(read_tail_finds_a_marker_that_a_whole_file_read_would_also_find) {
    const auto p = scratch("tail_marker.log");
    write_text(p, std::string(500000, 'x') + "\nFinished mdrun on rank 0 Mon Aug 3\n");
    CG_CHECK(read_tail(p, 64 * 1024).find("Finished mdrun") != std::string::npos);
    std::filesystem::remove(p);
}

CG_TEST(read_tail_of_a_missing_file_is_empty_rather_than_throwing) {
    CG_CHECK(read_tail(scratch("tail_absent.log"), 1024).empty());
}

CG_TEST(read_tail_of_an_empty_file_is_empty) {
    const auto p = scratch("tail_empty.log");
    write_text(p, "");
    CG_CHECK(read_tail(p, 1024).empty());
    std::filesystem::remove(p);
}
