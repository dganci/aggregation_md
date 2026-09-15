#include "MetadynamicsInternal.hpp"
#include "FileUtils.hpp"

#include <stdexcept>

namespace cg {

void require_file(const std::filesystem::path& path, const std::string& label, bool dry_run) {
    if (!dry_run && !std::filesystem::exists(path)) throw std::runtime_error("Missing " + label + ": " + path.string());
}

void copy_if_distinct(const std::filesystem::path& src, const std::filesystem::path& dst, bool dry_run) {
    if (dry_run) return;
    if (std::filesystem::exists(dst) && std::filesystem::equivalent(src, dst)) return;
    copy_overwrite(src, dst);
}

std::vector<std::string> walker_dirs(int n) {
    std::vector<std::string> dirs;
    dirs.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) dirs.push_back("walker" + std::to_string(i));
    return dirs;
}

} // namespace cg
