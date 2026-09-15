#pragma once

#include <filesystem>

namespace cg {

struct Config;
class Shell;
class GromacsDriver;

/// The fixed-length legs of the pipeline once a solvated system exists: energy
/// minimization, NVT/NPT equilibration, and (fixed-length) production, plus the
/// post-processing trjconv pass that centers the final trajectory.
class SimulationRunner {
public:
    SimulationRunner(Config& cfg, Shell& sh, GromacsDriver& gmx);

    void run_em();
    void analyze_em();

    void prepare_nvt();
    void run_nvt();
    /// Purely observational: appends temperature mean/stddev/drift for the
    /// completed NVT segment to diagnostics.jsonl.
    void analyze_nvt();

    void prepare_npt();
    void run_npt();
    void analyze_npt();

    void prepare_production();
    void run_production();
    /// Purely observational energy + structural diagnostics for the completed
    /// production segment (see analyze_nvt()'s no-throw note).
    void analyze_production();

    void center_trajectory();

private:
    Config& cfg_;
    Shell& sh_;
    GromacsDriver& gmx_;
};

} // namespace cg
