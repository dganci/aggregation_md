#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cgcv {

std::string shell_quote(const std::string& value);
std::string json_escape(const std::string& value);
std::string join(const std::vector<std::string>& values, const std::string& sep);
std::vector<int> parse_int_list(const std::string& value);
std::vector<int> parse_layers(const std::string& value);
std::string to_json_array(const std::vector<int>& values);
std::filesystem::path default_backend_path(const char* argv0);

} // namespace cgcv
