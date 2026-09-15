#pragma once

#include <string>
#include <vector>

namespace cg {

/// Small, dependency-free string helpers shared across the codebase.
std::string trim(const std::string& s);
bool starts_with(const std::string& s, const std::string& prefix);
bool ends_with(const std::string& s, const std::string& suffix);
std::vector<std::string> split_ws(const std::string& s);
/// Splits a comma-separated option value, trimming each item and dropping empty
/// ones.
std::vector<std::string> split_csv(const std::string& csv);
std::string join(const std::vector<std::string>& xs, const std::string& sep);
/// `value` repeated `n` times, comma-separated - PLUMED asks for one value per
/// component even where every component takes the same one (COEFFICIENTS,
/// POWERS, GRID_BIN).
std::string repeat_csv(const std::string& value, std::size_t n);
/// Wraps `s` in single quotes for POSIX shells, escaping embedded quotes (used
/// to build every external command line; see Shell::commandString).
std::string shell_quote(const std::string& s);
std::string basename_without_ext(const std::string& path);
std::string to_string_fixed(double x, int precision = 6);
/// Zero-pads `value` to `width` digits, e.g. zero_padded(7) == "007" (used for
/// adaptive-sampling chunk file names).
std::string zero_padded(int value, int width = 3);

} // namespace cg
