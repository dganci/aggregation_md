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

Preparer::Preparer(Config& cfg, Shell& sh, GromacsDriver& gmx) : cfg_(cfg), sh_(sh), gmx_(gmx) {}

void Preparer::run_insane() {
    sh_.run({cfg_.insane, "-o", cfg_.beadsPath().string(), "-p", cfg_.topologyPath().string(),
             "-f", cfg_.packCleanPath().string(), "-sol", "W", "-salt", std::to_string(cfg_.salt_M),
             "-pbc", "optimal", "-d", std::to_string(cfg_.box_margin_nm)});
    normalize_ions();
}

void Preparer::make_compact_box() {
    sh_.run({cfg_.gmx, "editconf", "-f", cfg_.packCleanPath().string(), "-o", cfg_.boxedPath().string(),
             "-bt", cfg_.box_type, "-d", std::to_string(cfg_.box_margin_nm), "-c"});
}

void Preparer::solvate_compact_box() {
    if (!sh_.dryRun() && !std::filesystem::exists(cfg_.martiniWaterPath())) {
        throw std::runtime_error("Missing Martini water coordinate file: " + cfg_.martiniWaterPath().string() +
                                 "\nProvide it with --martini-water-gro PATH");
    }

    sh_.run({cfg_.gmx, "solvate", "-cp", cfg_.boxedPath().string(), "-cs", cfg_.martiniWaterPath().string(),
             "-radius", std::to_string(cfg_.solvate_radius_nm), "-o", cfg_.solvatedPath().string(),
             "-p", cfg_.topologyPath().string()});
    if (!sh_.dryRun()) copy_overwrite(cfg_.solvatedPath(), cfg_.beadsPath());
}

std::string Preparer::make_ndx_and_find_water_group(const std::filesystem::path& gro_file,
                                                    const std::filesystem::path& ndx_path,
                                                    const std::filesystem::path& log_path) const {
    sh_.runShell("printf 'q\\n' | " + shell_quote(cfg_.gmx) +
                 " make_ndx -f " + shell_quote(gro_file.string()) +
                 " -o " + shell_quote(ndx_path.string()) +
                 " > " + shell_quote(log_path.string()),
                 cfg_.resultDir());

    if (sh_.dryRun()) return "W";

    for (const auto& line : read_lines(log_path)) {
        const auto tokens = split_ws(line);
        if (tokens.size() >= 2 && tokens[1] == "W") return tokens[0];
    }
    throw std::runtime_error("Cannot identify W group for genion.");
}

void Preparer::run_genion(const std::filesystem::path& tpr_path,
                          const std::filesystem::path& ndx_path,
                          const std::string& water_group,
                          bool include_salt) {
    std::ostringstream cmd;
    cmd << "printf '" << water_group << "\\n' | " << shell_quote(cfg_.gmx)
        << " genion -s " << shell_quote(tpr_path.string())
        << " -n " << shell_quote(ndx_path.string())
        << " -o " << shell_quote(cfg_.beadsPath().string())
        << " -p " << shell_quote(cfg_.topologyPath().string())
        << " -pname NA+ -nname CL- -neutral";
    if (include_salt && cfg_.salt_M > 0.0) cmd << " -conc " << cfg_.salt_M;
    sh_.runShell(cmd.str());
    normalize_ions();
}

void Preparer::add_ions_to_solvated_system() {
    if (!cfg_.add_ions) {
        normalize_ions();
        return;
    }

    const auto ndx = cfg_.resultDir() / "genion_tmp.ndx";
    const auto log = cfg_.resultDir() / "make_ndx.out";
    const auto water = make_ndx_and_find_water_group(cfg_.beadsPath(), ndx, log);
    gmx_.grompp(cfg_.systemDir() / "emin.mdp", cfg_.beadsPath(), cfg_.ionsTprPath(), {}, {}, ndx);
    run_genion(cfg_.ionsTprPath(), ndx, water, true);
    if (!sh_.dryRun()) { remove_if_exists(ndx); remove_if_exists(log); }
}

void Preparer::prepare_solvated_system() {
    if (cfg_.solvation_mode == "insane") {
        run_insane();
        patch_topology();
        verify_periodic_margin();
        verify_pmf_wall();
        return;
    }

    patch_topology();
    make_compact_box();
    solvate_compact_box();
    add_ions_to_solvated_system();
    verify_periodic_margin();
    verify_pmf_wall();
}

void Preparer::normalize_ions() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping GRO ion name normalization.\n";
        return;
    }

    const auto stats = normalize_gro_ion_names(cfg_.beadsPath());
    if (stats.total()) {
        std::cerr << "Normalized GRO ion names in " << cfg_.beadsPath() << ": "
                  << stats.atom_names << " atom name field(s), "
                  << stats.residue_names << " residue name field(s).\n";
    } else {
        std::cerr << "GRO ion names already match Martini ion names.\n";
    }
}

void Preparer::neutralize_if_needed() {
    if (!cfg_.is_phospho || !cfg_.neutralize_if_phospho) return;

    const auto ndx = cfg_.resultDir() / "genion_tmp.ndx";
    const auto log = cfg_.resultDir() / "make_ndx.out";
    const auto water = make_ndx_and_find_water_group(cfg_.beadsPath(), ndx, log);
    run_genion(cfg_.resultDir() / "em.tpr", ndx, water, false);
    if (!sh_.dryRun()) { remove_if_exists(ndx); remove_if_exists(log); }
}

void Preparer::patch_topology() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping topology patch because martinize2 was not executed.\n";
        return;
    }

    const auto itp_filename = find_first_itp(cfg_.topologyPath().parent_path(), cfg_.protomer_name);

    patch_martini_topology(
        cfg_.topologyPath(),
        cfg_.project_dir / "martini_v300",
        cfg_.protomer_name,
        cfg_.n_prot,
        itp_filename,
        cfg_.martiniMainItp()
    );
}

void Preparer::prepare_em() {
    normalize_ions();
    gmx_.grompp(cfg_.systemDir() / "emin.mdp", cfg_.beadsPath(), cfg_.resultDir() / "em.tpr");
}

} // namespace cg
