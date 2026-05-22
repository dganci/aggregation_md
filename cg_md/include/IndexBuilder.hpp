#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace cg {

class Shell;
struct Config;

using IndexGroups = std::map<std::string, std::vector<int>>;

std::string format_index_line(const std::vector<int>& atom_ids, int width = 15);
IndexGroups read_index(const std::filesystem::path& ndx_path);
void write_index(const IndexGroups& groups, const std::filesystem::path& ndx_path);
void generate_index(const Shell& sh, const Config& cfg, const std::filesystem::path& gro_file);

} // namespace cg
