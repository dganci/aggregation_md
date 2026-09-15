#include "Preparer.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "GroUtils.hpp"
#include "GromacsDriver.hpp"
#include "PdbUtils.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"
#include "TopologyEditor.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace cg {

void Preparer::write_packmol_input() {
    const auto r = cfg_.packmol_cluster_radius_A;
    const auto sources = cfg_.packmolSourcePaths();
    std::ostringstream out;
    out << "seed " << cfg_.seed << '\n'
        << "tolerance " << cfg_.packmol_tolerance_A << '\n'
        << "filetype pdb\n";
    if (sources.size() == 1) {
        out << "structure " << sources.front().string() << '\n'
            << "  number " << cfg_.n_prot << '\n'
            << "  inside sphere 0.0 0.0 0.0 " << r << '\n'
            << "  radius " << cfg_.packmol_radius_A << '\n'
            << "end structure\n";
    } else {
        for (const auto& src : sources) {
            out << "structure " << src.string() << '\n'
                << "  number 1\n"
                << "  inside sphere 0.0 0.0 0.0 " << r << '\n'
                << "  radius " << cfg_.packmol_radius_A << '\n'
                << "end structure\n";
        }
    }
    out << "output " << cfg_.packRawPath().string() << '\n';
    if (sh_.dryRun()) return;
    write_text(cfg_.packInputPath(), out.str());
}

void Preparer::run_packmol() {
    sh_.runShell(cfg_.packmol + " < " + shell_quote(cfg_.packInputPath().string()));
}

void Preparer::clean_packmol_pdb() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping Packmol PDB cleaning because " << cfg_.packRawPath() << " was not generated.\n";
        return;
    }

    std::ifstream in(cfg_.packRawPath());
    if (!in) throw std::runtime_error("Cannot read " + cfg_.packRawPath().string());

    std::ostringstream out;
    for (std::string line; std::getline(in, line);) {
        if (starts_with(line, "CONECT")) continue;
        if (starts_with(line, "HETATM")) line.replace(0, 6, "ATOM  ");
        out << line << '\n';
    }
    write_text(cfg_.packCleanPath(), out.str());
}

} // namespace cg
