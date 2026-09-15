#include "PdbUtils.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace cg {

const std::vector<std::string> kKnownModifiedResidues = {"SEP", "TPO", "PTR"};

namespace {

bool is_atom_record(const std::string& line) {
    return starts_with(line, "ATOM") || starts_with(line, "HETATM");
}

} // namespace

std::vector<std::string> detect_modified_residues(const std::filesystem::path& pdb_path) {
    std::vector<std::string> found;
    for (const auto& line : read_lines(pdb_path)) {
        if (line.size() < 20) continue;
        if (!is_atom_record(line)) continue;

        const auto resname = trim(line.substr(17, 3));
        for (const auto& known : kKnownModifiedResidues) {
            if (resname != known) continue;
            if (std::find(found.begin(), found.end(), resname) == found.end()) found.push_back(resname);
            break;
        }
    }
    return found;
}

std::string detect_modified_residue(const std::filesystem::path& pdb_path) {
    for (const auto& line : read_lines(pdb_path)) {
        if (line.size() < 20) continue;
        if (!is_atom_record(line)) continue;

        const auto resname = trim(line.substr(17, 3));
        for (const auto& known : kKnownModifiedResidues)
            if (resname == known) return resname;
    }
    return {};
}

int count_pdb_residues(const std::filesystem::path& pdb_path) {
    std::string last_key;
    int count = 0;
    for (const auto& line : read_lines(pdb_path)) {
        if (starts_with(line, "ENDMDL")) break;
        if (line.size() < 27) continue;
        if (!is_atom_record(line)) continue;

        const auto key = line.substr(21, 6);
        if (key != last_key) {
            ++count;
            last_key = key;
        }
    }
    return count;
}

int count_pdb_atoms(const std::filesystem::path& pdb_path) {
    int count = 0;
    for (const auto& line : read_lines(pdb_path)) {
        if (starts_with(line, "ENDMDL")) break;
        if (is_atom_record(line)) ++count;
    }
    return count;
}

std::vector<std::string> pdb_chain_ids(const std::filesystem::path& pdb_path) {
    std::vector<std::string> chains;
    for (const auto& line : read_lines(pdb_path)) {
        if (starts_with(line, "ENDMDL")) break;
        if (line.size() < 22) continue;
        if (!is_atom_record(line)) continue;

        const auto id = trim(line.substr(21, 1));
        if (id.empty()) continue;
        if (std::find(chains.begin(), chains.end(), id) == chains.end()) chains.push_back(id);
    }
    return chains;
}

double pdb_max_molecule_dmax(const std::filesystem::path& pdb_path, int atoms_per_molecule) {
    if (atoms_per_molecule <= 1) return 0.0;

    std::vector<std::array<double, 3>> xyz;
    for (const auto& line : read_lines(pdb_path)) {
        if (starts_with(line, "ENDMDL")) break;
        if (!is_atom_record(line) || line.size() < 54) continue;
        try {
            xyz.push_back({std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                           std::stod(line.substr(46, 8))});
        } catch (const std::exception&) {
            continue;
        }
    }

    const std::size_t per = static_cast<std::size_t>(atoms_per_molecule);
    double worst = 0.0;
    for (std::size_t base = 0; base + per <= xyz.size(); base += per) {
        for (std::size_t i = base; i < base + per; ++i) {
            for (std::size_t j = i + 1; j < base + per; ++j) {
                const double dx = xyz[i][0] - xyz[j][0];
                const double dy = xyz[i][1] - xyz[j][1];
                const double dz = xyz[i][2] - xyz[j][2];
                worst = std::max(worst, dx * dx + dy * dy + dz * dz);
            }
        }
    }
    return std::sqrt(worst);
}

std::pair<int, int> pdb_residue_number_range(const std::filesystem::path& pdb_path) {
    int lo = 0, hi = 0;
    bool seen = false;
    for (const auto& line : read_lines(pdb_path)) {
        if (starts_with(line, "ENDMDL")) break;
        if (!is_atom_record(line) || line.size() < 26) continue;
        int resid = 0;
        try {
            resid = std::stoi(trim(line.substr(22, 4)));
        } catch (const std::exception&) {
            continue;
        }
        if (!seen) { lo = hi = resid; seen = true; continue; }
        lo = std::min(lo, resid);
        hi = std::max(hi, resid);
    }
    return {lo, hi};
}

} // namespace cg
