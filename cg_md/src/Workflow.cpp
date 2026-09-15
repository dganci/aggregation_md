#include "Workflow.hpp"
#include "FileUtils.hpp"
#include "MdpWriter.hpp"
#include "PdbUtils.hpp"
#include "Reporter.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
namespace cg {

Workflow::Workflow(Config cfg)
    : cfg_(std::move(cfg)),
      sh_(cfg_.dry_run, true),
      gmx_(cfg_, sh_),
      preparer_(cfg_, sh_, gmx_),
      simulation_(cfg_, sh_, gmx_),
      adaptive_(cfg_, sh_, gmx_),
      metadynamics_(cfg_, sh_, gmx_),
      relax_(cfg_, sh_, gmx_) {
#ifdef _WIN32
    _putenv_s("GMX_MAXBACKUP", "-1");
#else
    setenv("GMX_MAXBACKUP", "-1", 1);
#endif
}

void Workflow::run_steps(std::initializer_list<Step> steps) {
    for (const auto& step : steps) step();
}
void Workflow::write_all_mdps() {
    write_em_mdp(cfg_);
    write_nvt_mdp(cfg_);
    write_npt_mdp(cfg_);
    write_md_mdp(cfg_);
}

void Workflow::report_existing_run() {
    if (cfg_.metadynamics) { run_report(cfg_, sh_, cfg_.metadDir(), "metadynamics"); return; }
    if (cfg_.adaptive) { run_report(cfg_, sh_, cfg_.resultDir(), "adaptive"); return; }

    if (std::filesystem::exists(cfg_.metadDir() / "diagnostics.jsonl") ||
        std::filesystem::exists(cfg_.metadDir() / "fes.dat")) {
        run_report(cfg_, sh_, cfg_.metadDir(), "metadynamics");
        return;
    }

    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(cfg_.resultDir(), ec)) {
        if (ec) break;
        if (starts_with(entry.path().filename().string(), "COLVAR_chunk_")) {
            run_report(cfg_, sh_, cfg_.resultDir(), "adaptive");
            return;
        }
    }

    run_report(cfg_, sh_, cfg_.resultDir(), "production");
}

void Workflow::build_and_equilibrate() {
    run_steps({
        [this] { preparer_.coarse_grain(); },
        [this] { if (cfg_.relax) relax_.run(); },
        [this] { preparer_.write_packmol_input(); },
        [this] { preparer_.run_packmol(); },
        [this] { preparer_.clean_packmol_pdb(); },
        [this] { preparer_.prepare_solvated_system(); },
        [this] { preparer_.prepare_em(); }
    });

    if (cfg_.solvation_mode == "insane" && cfg_.is_phospho && cfg_.neutralize_if_phospho) {
        preparer_.neutralize_if_needed();
        preparer_.prepare_em();
    }

    run_steps({
        [this] { simulation_.run_em(); },
        [this] { simulation_.analyze_em(); },
        [this] { simulation_.prepare_nvt(); },
        [this] { simulation_.run_nvt(); },
        [this] { simulation_.analyze_nvt(); },
        [this] { simulation_.prepare_npt(); },
        [this] { simulation_.run_npt(); },
        [this] { simulation_.analyze_npt(); }
    });
}

void Workflow::run_production() {
    if (cfg_.adaptive) { adaptive_.run(); return; }
    if (cfg_.metadynamics) { metadynamics_.run(); return; }
    run_steps({
        [this] { simulation_.prepare_production(); },
        [this] { simulation_.run_production(); },
        [this] { simulation_.analyze_production(); },
        [this] { simulation_.center_trajectory(); }
    });
}

bool Workflow::already_equilibrated() const {
    return std::filesystem::exists(cfg_.resultDir() / "npt.gro") &&
           std::filesystem::exists(cfg_.resultDir() / "npt.cpt");
}

void Workflow::run_all() {
    write_all_mdps();

    if (already_equilibrated()) {
        std::cerr << "npt.gro and npt.cpt already exist in " << cfg_.resultDir()
                  << " - this run was already prepared and equilibrated.\n"
                  << "Skipping coarse-graining/Packmol/solvation/EM/NVT/NPT and resuming"
                     " production directly.\n";
    } else {
        build_and_equilibrate();
    }

    run_production();
}

void Workflow::run() {
    ensure_dirs();
    cfg_.print();
    check_inputs();

    if (cfg_.stage == "all") {
        run_all();
        return;
    }

    const std::vector<Step> metad_steps = {[this] { write_all_mdps(); },
                                           [this] { metadynamics_.run(); }};

    const std::unordered_map<std::string, std::vector<Step>> stages = {
        {"cg", {[this] { preparer_.coarse_grain(); }}},
        {"prepare", {[this] { preparer_.coarse_grain(); }, [this] { preparer_.write_packmol_input(); },
                     [this] { preparer_.run_packmol(); }, [this] { preparer_.clean_packmol_pdb(); },
                     [this] { preparer_.prepare_solvated_system(); }, [this] { write_all_mdps(); }}},
        {"em", {[this] { write_all_mdps(); }, [this] { preparer_.prepare_em(); },
                [this] { simulation_.run_em(); }, [this] { simulation_.analyze_em(); }}},
        {"nvt", {[this] { write_all_mdps(); }, [this] { simulation_.prepare_nvt(); }, [this] { simulation_.run_nvt(); },
                 [this] { simulation_.analyze_nvt(); }}},
        {"npt", {[this] { write_all_mdps(); }, [this] { simulation_.prepare_npt(); },
                 [this] { simulation_.run_npt(); }, [this] { simulation_.analyze_npt(); }}},
        {"production", {[this] { write_all_mdps(); }, [this] { simulation_.prepare_production(); },
                        [this] { simulation_.run_production(); }, [this] { simulation_.analyze_production(); }}},
        {"adaptive", {[this] { write_all_mdps(); }, [this] { adaptive_.run(); }}},
        {"metadynamics", metad_steps},
        {"metad", metad_steps},
        {"pmf", metad_steps},
        {"relax", {[this] { write_all_mdps(); }, [this] { relax_.run(); }}},
        {"fes", {[this] { metadynamics_.finalize(); }}},
        {"center", {[this] { simulation_.center_trajectory(); }}},
        {"report", {[this] { report_existing_run(); }}},
        {"contact-map", {[this] {
             const auto adaptive_traj = cfg_.resultDir() / "md_all_center.xtc";
             const bool have_adaptive = std::filesystem::exists(adaptive_traj);
             run_contact_map(cfg_, sh_,
                             have_adaptive ? adaptive_traj : cfg_.resultDir() / "md_center.xtc",
                             have_adaptive ? cfg_.resultDir() / "md_all_center.gro"
                                           : cfg_.resultDir() / "md_center.gro",
                             cfg_.resultDir() / "report");
         }}}
    };

    const auto it = stages.find(cfg_.stage);
    if (it == stages.end()) throw std::runtime_error("Unknown stage: " + cfg_.stage);
    for (const auto& step : it->second) step();
}

} // namespace cg
