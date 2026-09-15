#pragma once

#include <filesystem>
#include <string>

namespace cgcv {

/// Throws std::runtime_error("Missing <label>: <path>") if `path` is not a
/// regular file.
void require_file(const std::filesystem::path& path, const std::string& label);
void ensure_dir(const std::filesystem::path& path);
/// Creates `path`'s parent directory if needed, then overwrites `path` with
/// `text`.
void write_text(const std::filesystem::path& path, const std::string& text);
std::string read_text(const std::filesystem::path& path);

} // namespace cgcv
