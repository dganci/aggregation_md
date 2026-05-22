#pragma once

#include <filesystem>
#include <string>

namespace cg {

std::string find_first_itp(const std::filesystem::path& topology_dir,
                           const std::string& prefix);

void patch_martini_topology(const std::filesystem::path& topology_path,
                            const std::string& molecule_name,
                            int n_prot,
                            const std::string& itp_filename);

} // namespace cg
