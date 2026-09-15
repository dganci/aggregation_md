#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

/// Small filesystem helpers used everywhere else in cg_md instead of raw
/// <fstream> calls, so every read/write failure raises the same
/// std::runtime_error("Cannot read/write <path>") message and directories are
/// created on demand before writing.
std::string read_text(const std::filesystem::path& path);

/// The last `max_bytes` of `path` (the whole file if it is smaller), for
/// looking near the end of a file that may be very large.
std::string read_tail(const std::filesystem::path& path, std::size_t max_bytes);
std::vector<std::string> read_lines(const std::filesystem::path& path);
void write_text(const std::filesystem::path& path, const std::string& text);
void write_lines(const std::filesystem::path& path, const std::vector<std::string>& lines);
void copy_overwrite(const std::filesystem::path& from, const std::filesystem::path& to);
void remove_if_exists(const std::filesystem::path& path);

} // namespace cg
