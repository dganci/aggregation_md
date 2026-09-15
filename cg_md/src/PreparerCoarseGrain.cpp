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
namespace {

std::string parent_residue_for_modification(const std::string& resname) {
    if (resname == "SEP") return "SER";
    if (resname == "TPO") return "THR";
    if (resname == "PTR") return "TYR";
    return resname;
}

std::filesystem::path write_phospho_repaired_pdb(const std::filesystem::path& src,
                                                  const std::filesystem::path& dst) {
    auto lines = read_lines(src);
    for (auto& line : lines) {
        if (line.size() < 20) continue;
        if (!(starts_with(line, "ATOM") || starts_with(line, "HETATM"))) continue;

        const auto resname = trim(line.substr(17, 3));
        const auto parent = parent_residue_for_modification(resname);
        if (parent != resname) line.replace(17, 3, parent);
    }
    write_lines(dst, lines);
    return dst;
}
} // namespace

void Preparer::coarse_grain() {
    const std::string ss = !cfg_.ss_string.empty()
        ? cfg_.ss_string
        : std::string(static_cast<std::size_t>(cfg_.seq_length), 'C');

    std::filesystem::path martinize_input = cfg_.pdbPath();
    if (cfg_.is_phospho) {
        martinize_input = cfg_.pdbPath().parent_path() / (cfg_.protomer_name + "_phosphorepaired.pdb");
        write_phospho_repaired_pdb(cfg_.pdbPath(), martinize_input);
    }

    std::vector<std::string> cmd = {
        cfg_.martinize, "-f", martinize_input.string(), "-x", cfg_.cgPath().string(),
        "-o", cfg_.topologyPath().string(), "-name", cfg_.protomer_name,
        "-ff", "martini3IDP"
    };

    if (cfg_.use_dssp) {
        cmd.push_back("-dssp");
    } else {
        cmd.insert(cmd.end(), {"-ss", ss});
    }

    if (cfg_.is_phospho) {
        cmd.insert(cmd.end(), {"-bonds-from", "both"});
    }

    const auto ff_dir = cfg_.martinizeFfDirPath();
    const auto map_dir = cfg_.martinizeMapDirPath();
    if (std::filesystem::is_directory(ff_dir)) cmd.insert(cmd.end(), {"-ff-dir", ff_dir.string()});
    if (std::filesystem::is_directory(map_dir)) cmd.insert(cmd.end(), {"-map-dir", map_dir.string()});

    cmd.insert(cmd.end(), {"-p", "backbone", "-cys", "auto"});

    if (!cfg_.merge_chains.empty()) cmd.insert(cmd.end(), {"-merge", cfg_.merge_chains});

    if (cfg_.use_elastic) {
        cmd.insert(cmd.end(), {"-elastic", "-ef", "700.0", "-el", "0", "-eu", "0.85"});
        if (!trim(cfg_.elastic_units).empty()) {
            cmd.insert(cmd.end(), {"-eunit", trim(cfg_.elastic_units)});
        }
    }

    if (!cfg_.modify.empty()) {
        cmd.insert(cmd.end(), {"-modify", cfg_.modify});
    }

    sh_.run(cmd, cfg_.project_dir);
    move_generated_itps();
    verify_bead_count();
    verify_lambda_coverage();
}

void Preparer::move_generated_itps() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping ITP move after martinize2.\n";
        return;
    }

    const auto prefix = cfg_.protomer_name + "_";
    for (const auto& entry : std::filesystem::directory_iterator(cfg_.project_dir)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".itp") continue;
        if (!starts_with(entry.path().filename().string(), prefix)) continue;
        const auto dst = cfg_.topologyPath().parent_path() / entry.path().filename();
        remove_if_exists(dst);
        std::filesystem::rename(entry.path(), dst);
    }
}

} // namespace cg
