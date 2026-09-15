#pragma once

#include <filesystem>
#include <string>

namespace cg {

/// Returns the first `*.itp` filename in `topology_dir` starting with `prefix`
/// (used to find the .itp martinize2 generated for the protomer).
std::string find_first_itp(const std::filesystem::path& topology_dir,
                           const std::string& prefix);

/// Rewrites the martinize2-generated .top file at `topology_path` so it
/// includes the Martini force-field/ion/solvent .itp files from `martini_dir`
/// plus the protomer's own `itp_filename`.
void patch_martini_topology(const std::filesystem::path& topology_path,
                            const std::filesystem::path& martini_dir,
                            const std::string& molecule_name,
                            int n_prot,
                            const std::string& itp_filename,
                            const std::string& main_itp = "martini_v3.0.0.itp");

/// Copies the topology at `src` to `dst` keeping exactly ONE molecule entry:
/// the first one, with a count of 1.
void write_single_molecule_topology(const std::filesystem::path& src,
                                    const std::filesystem::path& dst);

} // namespace cg
