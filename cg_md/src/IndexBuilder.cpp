#include "IndexBuilder.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace cg {
namespace {

std::vector<int> sorted_unique(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

std::vector<int> atom_range(int first, int count) {
    std::vector<int> atoms(static_cast<std::size_t>(count));
    std::iota(atoms.begin(), atoms.end(), first);
    return atoms;
}

} // namespace

std::string format_index_line(const std::vector<int>& atom_ids, int width) {
    std::ostringstream out;
    for (std::size_t i = 0; i < atom_ids.size(); i += static_cast<std::size_t>(width)) {
        if (i) out << '\n';
        const auto stop = std::min(atom_ids.size(), i + static_cast<std::size_t>(width));
        for (auto j = i; j < stop; ++j) out << (j == i ? "" : " ") << atom_ids[j];
    }
    return out.str();
}

IndexGroups read_index(const std::filesystem::path& ndx_path) {
    std::ifstream in(ndx_path);
    if (!in) throw std::runtime_error("Cannot read index file " + ndx_path.string());

    IndexGroups groups;
    std::string current;
    for (std::string line; std::getline(in, line);) {
        line = trim(line);
        if (line.empty()) continue;
        if (starts_with(line, "[") && ends_with(line, "]")) {
            current = trim(line.substr(1, line.size() - 2));
            groups[current];
        } else if (!current.empty()) {
            for (const auto& x : split_ws(line)) groups[current].push_back(std::stoi(x));
        }
    }
    return groups;
}

void write_index(const IndexGroups& groups, const std::filesystem::path& ndx_path) {
    std::ostringstream out;
    for (const auto& [name, atoms] : groups) {
        out << "[ " << name << " ]\n" << format_index_line(atoms) << "\n\n";
    }
    write_text(ndx_path, out.str());
}

void generate_index(const Shell& sh, const Config& cfg, const std::filesystem::path& gro_file) {
    const auto ndx_path = cfg.resultDir() / "index.ndx";
    const auto input_path = cfg.resultDir() / "tmp_ndx_input.txt";
    write_text(input_path, "q\n");

    sh.run({"/bin/bash", "-lc", shell_quote(cfg.gmx) +
        " make_ndx -f " + shell_quote(gro_file.string()) +
        " -o " + shell_quote(ndx_path.string()) +
        " < " + shell_quote(input_path.string())});

    if (sh.dryRun()) {
        std::cerr << "[dry-run] Skipping index parsing because " << ndx_path << " was not generated.\n";
        return;
    }

    auto groups = read_index(ndx_path);
    std::vector<int> protein;
    protein.reserve(static_cast<std::size_t>(cfg.n_prot * cfg.atoms_per_prot));

    for (int i = 0, first = 1; i < cfg.n_prot; ++i, first += cfg.atoms_per_prot) {
        auto atoms = atom_range(first, cfg.atoms_per_prot);
        protein.insert(protein.end(), atoms.begin(), atoms.end());
        groups["Protein" + std::to_string(i + 1)] = std::move(atoms);
    }

    protein = sorted_unique(std::move(protein));
    groups["Protein"] = protein;
    groups["AllProteins"] = protein;

    if (const auto it = groups.find("System"); it != groups.end()) {
        const std::set<int> protein_set(protein.begin(), protein.end());
        std::vector<int> solvent;
        std::copy_if(it->second.begin(), it->second.end(), std::back_inserter(solvent),
                     [&](int a) { return !protein_set.count(a); });
        groups["Solvent_and_ions"] = sorted_unique(std::move(solvent));
    }

    write_index(groups, ndx_path);
    remove_if_exists(input_path);
}

} // namespace cg
