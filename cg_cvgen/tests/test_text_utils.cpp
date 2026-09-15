#include "test_framework.hpp"
#include "TextUtils.hpp"

#include <exception>
#include <vector>

using namespace cgcv;

CG_TEST(shell_quote_wraps_and_escapes) {
    CG_CHECK_EQ(shell_quote("plain"), "'plain'");
    CG_CHECK_EQ(shell_quote("it's"), "'it'\\''s'");
    CG_CHECK_EQ(shell_quote(""), "''");
}

CG_TEST(json_escape_handles_special_characters) {
    CG_CHECK_EQ(json_escape("a\\b"), "a\\\\b");
    CG_CHECK_EQ(json_escape("a\"b"), "a\\\"b");
    CG_CHECK_EQ(json_escape("a\nb"), "a\\nb");
    CG_CHECK_EQ(json_escape("plain"), "plain");
}

CG_TEST(join_inserts_separator_between_elements_only) {
    CG_CHECK_EQ(join({"a", "b", "c"}, ","), "a,b,c");
    CG_CHECK_EQ(join({"only"}, ","), "only");
    CG_CHECK_EQ(join({}, ","), "");
}

CG_TEST(parse_int_list_accepts_comma_and_semicolon_separators) {
    CG_CHECK_EQ(parse_int_list("5,7,10"), (std::vector<int>{5, 7, 10}));
    CG_CHECK_EQ(parse_int_list("5;7;10"), (std::vector<int>{5, 7, 10}));
    CG_CHECK_EQ(parse_int_list(" 5 , 7 , 10 "), (std::vector<int>{5, 7, 10}));
}

CG_TEST(parse_int_list_rejects_empty_and_non_positive_values) {
    bool threw_empty = false;
    try { parse_int_list(""); } catch (const std::exception&) { threw_empty = true; }
    CG_CHECK(threw_empty);

    bool threw_zero = false;
    try { parse_int_list("0,1"); } catch (const std::exception&) { threw_zero = true; }
    CG_CHECK(threw_zero);

    bool threw_negative = false;
    try { parse_int_list("-1"); } catch (const std::exception&) { threw_negative = true; }
    CG_CHECK(threw_negative);
}

CG_TEST(parse_layers_is_an_alias_of_parse_int_list) {
    CG_CHECK_EQ(parse_layers("32,16"), (std::vector<int>{32, 16}));
}

CG_TEST(to_json_array_formats_as_json) {
    CG_CHECK_EQ(to_json_array({1, 2, 3}), "[1, 2, 3]");
    CG_CHECK_EQ(to_json_array({}), "[]");
}
