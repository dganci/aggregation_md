#pragma once

#include <string>
#include <vector>

namespace cg {

std::string trim(const std::string& s);
bool starts_with(const std::string& s, const std::string& prefix);
bool ends_with(const std::string& s, const std::string& suffix);
std::vector<std::string> split_ws(const std::string& s);
std::string join(const std::vector<std::string>& xs, const std::string& sep);
std::string shell_quote(const std::string& s);
std::string basename_without_ext(const std::string& path);
std::string to_string_fixed(double x, int precision = 6);
std::string zero_padded(int value, int width = 3);

} // namespace cg
