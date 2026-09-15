#include "test_framework.hpp"
#include "StringUtils.hpp"

using namespace cg;

CG_TEST(trim_removes_surrounding_whitespace) {
    CG_CHECK_EQ(trim("  hello  "), "hello");
    CG_CHECK_EQ(trim("\t\nhi\r\n"), "hi");
    CG_CHECK_EQ(trim(""), "");
    CG_CHECK_EQ(trim("   "), "");
}

CG_TEST(starts_with_and_ends_with_basic_cases) {
    CG_CHECK(starts_with("martinize2", "martin"));
    CG_CHECK(!starts_with("martinize2", "2martin"));
    CG_CHECK(ends_with("desmin_head.pdb", ".pdb"));
    CG_CHECK(!ends_with("desmin_head.pdb", ".itp"));
    CG_CHECK(starts_with("abc", "abc"));
    CG_CHECK(ends_with("abc", ""));
    CG_CHECK(!starts_with("ab", "abc"));
}

CG_TEST(split_ws_splits_on_any_whitespace_and_ignores_runs) {
    const auto tokens = split_ws("  a  b\tc\n d ");
    CG_CHECK_EQ(tokens.size(), std::size_t{4});
    CG_CHECK_EQ(tokens[0], "a");
    CG_CHECK_EQ(tokens[1], "b");
    CG_CHECK_EQ(tokens[2], "c");
    CG_CHECK_EQ(tokens[3], "d");
}

CG_TEST(join_inserts_separator_between_elements_only) {
    CG_CHECK_EQ(join({"a", "b", "c"}, ","), "a,b,c");
    CG_CHECK_EQ(join({"only"}, ","), "only");
    CG_CHECK_EQ(join({}, ","), "");
}

CG_TEST(shell_quote_escapes_embedded_single_quotes) {
    CG_CHECK_EQ(shell_quote("simple"), "'simple'");
    CG_CHECK_EQ(shell_quote("it's"), "'it'\\''s'");
    CG_CHECK_EQ(shell_quote(""), "''");
}

CG_TEST(basename_without_ext_strips_directory_and_extension) {
    CG_CHECK_EQ(basename_without_ext("/a/b/desmin_head.pdb"), "desmin_head");
    CG_CHECK_EQ(basename_without_ext("desmin_head.pdb"), "desmin_head");
    CG_CHECK_EQ(basename_without_ext("noext"), "noext");
    CG_CHECK_EQ(basename_without_ext("C:\\a\\b\\file.itp"), "file");
}

CG_TEST(zero_padded_pads_to_requested_width) {
    CG_CHECK_EQ(zero_padded(7), "007");
    CG_CHECK_EQ(zero_padded(42, 2), "42");
    CG_CHECK_EQ(zero_padded(0, 4), "0000");
}
