#include "RelaxRunner.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "GromacsDriver.hpp"
#include "GroUtils.hpp"
#include "MdpWriter.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"
#include "TopologyEditor.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace cg {

RelaxRunner::RelaxRunner(Config& cfg, Shell& sh, GromacsDriver& gmx) : cfg_(cfg), sh_(sh), gmx_(gmx) {}

std::filesystem::path RelaxRunner::dir() const {
    return cfg_.project_dir / "runs" / ("relax_" + cfg_.protomer_name);
}

std::filesystem::path RelaxRunner::topology() const {
    return cfg_.topologyPath().parent_path() / ("relax_" + cfg_.protomer_name + ".top");
}

void RelaxRunner::build_system() {
    if (!sh_.dryRun()) std::filesystem::create_directories(dir());

    const auto cg = cfg_.cgPath();
    if (!sh_.dryRun() && !std::filesystem::exists(cg))
        throw std::runtime_error("Relaxation needs the coarse-grained protomer: " + cg.string() +
                                 ". Run --stage prepare first.");

    const auto topo = topology();
    if (!sh_.dryRun() && cfg_.solvation_mode != "insane")
        write_single_molecule_topology(cfg_.topologyPath(), topo);

    if (cfg_.solvation_mode == "insane") {
        sh_.run({cfg_.insane, "-o", (dir() / "solvated.gro").string(),
                 "-p", topo.string(), "-f", cg.string(),
                 "-sol", "W", "-salt", to_string_fixed(cfg_.salt_M, 3),
                 "-pbc", "optimal", "-d", to_string_fixed(cfg_.box_margin_nm, 3)});

        if (!sh_.dryRun()) {
            patch_martini_topology(topo, cfg_.project_dir / "martini_v300", cfg_.protomer_name,
                                   /*n_prot=*/1,
                                   find_first_itp(cfg_.topologyPath().parent_path(),
                                                  cfg_.protomer_name),
                                   cfg_.martiniMainItp());
        }
    } else {
        sh_.run({cfg_.gmx, "editconf", "-f", cg.string(), "-o", (dir() / "boxed.gro").string(),
                 "-bt", cfg_.box_type, "-d", to_string_fixed(cfg_.box_margin_nm, 3), "-c"});

        sh_.run({cfg_.gmx, "solvate", "-cp", (dir() / "boxed.gro").string(),
                 "-cs", cfg_.martiniWaterPath().string(),
                 "-radius", to_string_fixed(cfg_.solvate_radius_nm, 3),
                 "-o", (dir() / "solvated.gro").string(), "-p", topo.string()});
    }

    if (!sh_.dryRun()) normalize_gro_ion_names(dir() / "solvated.gro");
}

void RelaxRunner::equilibrate() {
    const auto topo = topology();
    const auto run_gmx = [&](const std::string& mdp, const std::string& in, const std::string& out,
                             const std::string& cpt) {
        std::vector<std::string> cmd = {cfg_.gmx, "grompp", "-f", (cfg_.systemDir() / mdp).string(),
                                        "-c", (dir() / in).string(), "-p", topo.string(),
                                        "-o", (dir() / (out + ".tpr")).string(),
                                        "-po", (dir() / (out + "_mdout.mdp")).string()};
        if (!cpt.empty()) cmd.insert(cmd.end(), {"-t", (dir() / cpt).string()});
        sh_.run(cmd);
        gmx_.mdrun(dir() / out);
    };

    run_gmx("emin.mdp", "solvated.gro", "em", "");
    run_gmx("nvt.mdp", "em.gro", "nvt", "");
    run_gmx("npt.mdp", "nvt.gro", "npt", "nvt.cpt");
}

void RelaxRunner::produce() {
    const auto topo = topology();
    const auto steps = static_cast<std::int64_t>(cfg_.relax_us * 1'000'000.0 / cfg_.md_dt_ps + 0.5);
    if (steps <= 0) throw std::runtime_error("--relax-us is shorter than one --md-dt-ps step");

    const auto mdp = write_relax_mdp(cfg_, steps);

    sh_.run({cfg_.gmx, "grompp", "-f", mdp.string(), "-c", (dir() / "npt.gro").string(),
             "-p", topo.string(), "-t", (dir() / "npt.cpt").string(),
             "-o", (dir() / "md.tpr").string(),
             "-po", (dir() / "md_mdout.mdp").string()});
    gmx_.mdrun(dir() / "md");
}

void RelaxRunner::run() {
    std::cerr << "\n============================================================\n"
              << "RELAXING A SINGLE PROTOMER (" << cfg_.relax_us << " us) BEFORE PACKING\n"
              << "============================================================\n";
    build_system();
    equilibrate();
    produce();
    gmx_.trjconv(dir() / "md.tpr", dir() / "md.xtc", dir() / "md_prot.xtc",
                 "1\n", {"-pbc", "whole"});
    if (cfg_.distinct_copies) write_distinct_copies();
    else select_representative();
}

} // namespace cg
