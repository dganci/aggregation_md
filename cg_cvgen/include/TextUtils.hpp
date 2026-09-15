#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cgcv {

/// Wraps `value` in single quotes for POSIX shells, escaping embedded quotes.
std::string shell_quote(const std::string& value);
/// Escapes backslash/quote/newline/CR/tab for embedding `value` in a JSON
/// string.
std::string json_escape(const std::string& value);
std::string join(const std::vector<std::string>& values, const std::string& sep);
/// Parses a comma- or semicolon-separated list of positive integers, e.g.
/// "5,7,10" or "5;7;10" -> {5, 7, 10}.
std::vector<int> parse_int_list(const std::string& value);
/// Alias of parse_int_list() used for --hidden-layers (kept as a distinct name
/// for readability at call sites; same validation rules apply).
std::vector<int> parse_layers(const std::string& value);
std::string to_json_array(const std::vector<int>& values);
/// Locates scripts/cvgen_backend.py relative to the current working directory
/// or the cg_cvgen executable's directory (searching up to two parent
/// directories, to work whether the executable lives directly in the project
/// root or in a build/ subdirectory).
std::filesystem::path default_backend_path(const char* argv0);

} // namespace cgcv
