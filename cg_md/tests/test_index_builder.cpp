#include "test_framework.hpp"
#include "IndexBuilder.hpp"
#include "FileUtils.hpp"

#include <filesystem>
#include <vector>

using namespace cg;

CG_TEST(format_index_line_wraps_at_requested_width) {
    CG_CHECK_EQ(format_index_line({1, 2, 3}, 15), "1 2 3");
    CG_CHECK_EQ(format_index_line({1, 2, 3, 4}, 2), "1 2\n3 4");
    CG_CHECK_EQ(format_index_line({}, 15), "");
}

CG_TEST(read_index_and_write_index_round_trip) {
    const auto path = std::filesystem::temp_directory_path() / "cg_md_test_index.ndx";

    IndexGroups groups;
    groups["Protein1"] = {1, 2, 3};
    groups["System"] = {1, 2, 3, 4, 5};

    write_index(groups, path);
    const auto read_back = read_index(path);

    CG_CHECK_EQ(read_back.at("Protein1"), (std::vector<int>{1, 2, 3}));
    CG_CHECK_EQ(read_back.at("System"), (std::vector<int>{1, 2, 3, 4, 5}));

    std::filesystem::remove(path);
}
