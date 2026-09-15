#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace cg {

class Shell;
struct Config;

using IndexGroups = std::map<std::string, std::vector<int>>;

/// Wraps `atom_ids` onto lines of at most `width` entries, GROMACS .ndx style.
std::string format_index_line(const std::vector<int>& atom_ids, int width = 15);
IndexGroups read_index(const std::filesystem::path& ndx_path);
void write_index(const IndexGroups& groups, const std::filesystem::path& ndx_path);

/// Runs `gmx make_ndx` on `gro_file`, then adds one "ProteinN" group per
/// protomer (assuming protomers are laid out back-to-back with
/// Config::atoms_per_prot atoms each), plus aggregate "Protein"/"AllProteins"
/// and "Solvent_and_ions" groups.
void generate_index(const Shell& sh, const Config& cfg, const std::filesystem::path& gro_file);

} // namespace cg
