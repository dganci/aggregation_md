#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

inline constexpr const char* kHillsFile = "HILLS";
inline constexpr const char* kColvarFile = "COLVAR";

void require_file(const std::filesystem::path& path, const std::string& label, bool dry_run);
void copy_if_distinct(const std::filesystem::path& src, const std::filesystem::path& dst, bool dry_run);
std::vector<std::string> walker_dirs(int n);

} // namespace cg
