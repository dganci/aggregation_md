#include "Workflow.hpp"
#include "MdpWriter.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cg {

Workflow::Workflow(Config cfg) : cfg_(std::move(cfg)), sh_(cfg_.dry_run, true) {
#ifdef _WIN32
    _putenv_s("GMX_MAXBACKUP", "-1");
#else
    setenv("GMX_MAXBACKUP", "-1", 1);
#endif
}

void Workflow::run_steps(std::initializer_list<Step> steps) {
    for (const auto step : steps) (this->*step)();
}

void Workflow::run_all() {
    run_steps({
        &Workflow::coarse_grain,
        &Workflow::write_packmol_input,
        &Workflow::run_packmol,
        &Workflow::clean_packmol_pdb,
        &Workflow::prepare_solvated_system,
        &Workflow::prepare_em
    });

    if (cfg_.solvation_mode == "insane" && cfg_.is_phospho && cfg_.neutralize_if_phospho) {
        neutralize_if_needed();
        prepare_em();
    }

    run_steps({
        &Workflow::run_em,
        &Workflow::analyze_em,
        &Workflow::prepare_nvt,
        &Workflow::run_nvt,
        &Workflow::prepare_npt,
        &Workflow::run_npt,
        &Workflow::analyze_npt
    });

    if (cfg_.adaptive) {
        run_adaptive_production();
    } else if (cfg_.metadynamics) {
        run_metadynamics_workflow();
    } else {
        run_steps({&Workflow::prepare_production, &Workflow::run_production, &Workflow::center_trajectory});
    }
}

void Workflow::run() {
    ensure_dirs();
    cfg_.print();
    check_inputs();

    if (cfg_.stage == "all") {
        run_all();
        return;
    }

    const std::unordered_map<std::string, std::vector<Step>> stages = {
        {"prepare", {&Workflow::coarse_grain, &Workflow::write_packmol_input, &Workflow::run_packmol,
                     &Workflow::clean_packmol_pdb, &Workflow::prepare_solvated_system}},
        {"em", {&Workflow::write_all_mdps, &Workflow::prepare_em, &Workflow::run_em, &Workflow::analyze_em}},
        {"nvt", {&Workflow::write_all_mdps, &Workflow::prepare_nvt, &Workflow::run_nvt}},
        {"npt", {&Workflow::write_all_mdps, &Workflow::prepare_npt, &Workflow::run_npt, &Workflow::analyze_npt}},
        {"production", {&Workflow::write_all_mdps, &Workflow::prepare_production, &Workflow::run_production}},
        {"adaptive", {&Workflow::write_all_mdps, &Workflow::run_adaptive_production}},
        {"metadynamics", {&Workflow::write_all_mdps, &Workflow::run_metadynamics_workflow}},
        {"metad", {&Workflow::write_all_mdps, &Workflow::run_metadynamics_workflow}},
        {"fes", {&Workflow::finalize_metadynamics}},
        {"center", {&Workflow::center_trajectory}}
    };

    const auto it = stages.find(cfg_.stage);
    if (it == stages.end()) throw std::runtime_error("Unknown stage: " + cfg_.stage);
    for (const auto step : it->second) (this->*step)();
}

void Workflow::write_all_mdps() {
    write_em_mdp(cfg_);
    write_nvt_mdp(cfg_);
    write_npt_mdp(cfg_);
    write_md_mdp(cfg_);
    if (cfg_.metadynamics) write_metad_mdp(cfg_);
}

} // namespace cg
