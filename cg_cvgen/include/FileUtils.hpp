#pragma once

#include <filesystem>
#include <string>

namespace cgcv {

void require_file(const std::filesystem::path& path, const std::string& label);
void ensure_dir(const std::filesystem::path& path);
void write_text(const std::filesystem::path& path, const std::string& text);
std::string read_text(const std::filesystem::path& path);

} // namespace cgcv
