#include "SimulationRunner.hpp"
#include "Config.hpp"
#include "GromacsDriver.hpp"
#include "IndexBuilder.hpp"
#include "PlumedWriter.hpp"
#include "Reporter.hpp"
#include "RunDiagnostics.hpp"
#include "Shell.hpp"

#include <exception>
#include <iostream>

namespace cg {
namespace {

template <class Fn>
void run_diagnostics_safely(bool dry_run, const char* label, Fn&& fn) {
    if (dry_run) return;
    try {
        fn();
    } catch (const std::exception& e) {
        std::cerr << "Warning: " << label << " diagnostics failed (ignored, purely observational): "
                  << e.what() << '\n';
    }
}

} // namespace

SimulationRunner::SimulationRunner(Config& cfg, Shell& sh, GromacsDriver& gmx) : cfg_(cfg), sh_(sh), gmx_(gmx) {}

void SimulationRunner::run_em() {
    gmx_.mdrun(cfg_.resultDir() / "em", cfg_.resultDir() / "em.tpr", {}, 1);
}

void SimulationRunner::analyze_em() {
    gmx_.energy(cfg_.resultDir() / "em.edr", cfg_.resultDir() / "potential.xvg", "Potential");
}

void SimulationRunner::prepare_nvt() {
    generate_index(sh_, cfg_, cfg_.resultDir() / "em.gro");
    gmx_.grompp(cfg_.systemDir() / "nvt.mdp", cfg_.resultDir() / "em.gro", cfg_.resultDir() / "nvt.tpr", {}, {}, cfg_.resultDir() / "index.ndx");
}

void SimulationRunner::run_nvt() {
    gmx_.mdrun(cfg_.resultDir() / "nvt");
}

void SimulationRunner::analyze_nvt() {
    run_diagnostics_safely(sh_.dryRun(), "NVT", [&] {
        RunDiagnostics diag;
        diag.energy = compute_energy_diagnostics(gmx_, cfg_.resultDir() / "nvt.edr",
                                                  cfg_.resultDir() / "diagnostics_scratch",
                                                  /*include_pressure=*/false, /*include_density=*/false);
        append_diagnostics_jsonl(cfg_.resultDir() / "diagnostics.jsonl", diag, "nvt");
    });
}

void SimulationRunner::prepare_npt() {
    generate_index(sh_, cfg_, cfg_.resultDir() / "nvt.gro");
    gmx_.grompp(cfg_.systemDir() / "npt.mdp", cfg_.resultDir() / "nvt.gro", cfg_.resultDir() / "npt.tpr",
                cfg_.resultDir() / "nvt.cpt", cfg_.resultDir() / "nvt.gro", cfg_.resultDir() / "index.ndx");
}

void SimulationRunner::run_npt() {
    gmx_.mdrun(cfg_.resultDir() / "npt");
}

void SimulationRunner::analyze_npt() {
    run_diagnostics_safely(sh_.dryRun(), "NPT", [&] {
        RunDiagnostics diag;
        diag.energy = compute_energy_diagnostics(gmx_, cfg_.resultDir() / "npt.edr",
                                                  cfg_.resultDir() / "diagnostics_scratch",
                                                  /*include_pressure=*/true, /*include_density=*/true);
        append_diagnostics_jsonl(cfg_.resultDir() / "diagnostics.jsonl", diag, "npt");
    });
}

void SimulationRunner::prepare_production() {
    generate_index(sh_, cfg_, cfg_.resultDir() / "npt.gro");
    write_plumed_dat(cfg_);
    gmx_.grompp(cfg_.systemDir() / "md.mdp", cfg_.resultDir() / "npt.gro", cfg_.resultDir() / "md.tpr",
                cfg_.resultDir() / "npt.cpt", {}, cfg_.resultDir() / "index.ndx");
}

void SimulationRunner::run_production() {
    gmx_.mdrun(cfg_.resultDir() / "md", {}, cfg_.resultDir() / "plumed.dat");
}

void SimulationRunner::analyze_production() {
    run_diagnostics_safely(sh_.dryRun(), "production", [&] {
        RunDiagnostics diag;
        diag.energy = compute_energy_diagnostics(gmx_, cfg_.resultDir() / "md.edr",
                                                  cfg_.resultDir() / "diagnostics_scratch",
                                                  /*include_pressure=*/true, /*include_density=*/true);
        const double dt_colvar_ps = cfg_.md_dt_ps * static_cast<double>(cfg_.plumed_stride);
        diag.structure = compute_structural_diagnostics({cfg_.resultDir() / "COLVAR"}, cfg_.n_prot,
                                                         dt_colvar_ps, cfg_.adaptive_pair_contact_threshold);
        diag.has_structure = true;
        append_diagnostics_jsonl(cfg_.resultDir() / "diagnostics.jsonl", diag, "production");
    });

    run_report(cfg_, sh_, cfg_.resultDir(), "production");
}

void SimulationRunner::center_trajectory() {
    const auto tpr = cfg_.resultDir() / "md.tpr";
    const auto xtc = cfg_.resultDir() / "md.xtc";
    const auto nojump = cfg_.resultDir() / "md_nojump.xtc";
    const auto centered_xtc = cfg_.resultDir() / "md_center.xtc";
    const auto centered_gro = cfg_.resultDir() / "md_center.gro";

    gmx_.trjconv(tpr, xtc, nojump, "0\n", {"-pbc nojump"});

    gmx_.trjconv(tpr, nojump, centered_xtc,
                 "1\n0\n", {"-center", "-pbc mol", "-ur compact"});

    gmx_.trjconv(tpr, centered_xtc, centered_gro,
                 "0\n", {"-dump", "0"});
}

} // namespace cg
