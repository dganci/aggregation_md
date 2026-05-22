#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

std::string read_text(const std::filesystem::path& path);
std::vector<std::string> read_lines(const std::filesystem::path& path);
void write_text(const std::filesystem::path& path, const std::string& text);
void write_lines(const std::filesystem::path& path, const std::vector<std::string>& lines);
void copy_overwrite(const std::filesystem::path& from, const std::filesystem::path& to);
void remove_if_exists(const std::filesystem::path& path);

} // namespace cg
