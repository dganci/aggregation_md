#pragma once

#include "AdaptiveSampler.hpp"
#include "Config.hpp"
#include "GromacsDriver.hpp"
#include "MetadynamicsRunner.hpp"
#include "Preparer.hpp"
#include "RelaxRunner.hpp"
#include "Shell.hpp"
#include "SimulationRunner.hpp"

#include <functional>
#include <initializer_list>

namespace cg {

/// Top-level orchestrator: owns the Config and Shell for one invocation, and
/// dispatches each pipeline stage to the collaborator that implements it
/// (Preparer, SimulationRunner, AdaptiveSampler, MetadynamicsRunner).
class Workflow {
public:
    explicit Workflow(Config cfg);
    void run();

private:
    using Step = std::function<void()>;

    Config cfg_;
    Shell sh_;
    GromacsDriver gmx_;
    Preparer preparer_;
    SimulationRunner simulation_;
    AdaptiveSampler adaptive_;
    MetadynamicsRunner metadynamics_;
    RelaxRunner relax_;

    void run_all();
    /// Coarse-graining -> (relaxation) -> Packmol -> solvation -> EM/NVT/NPT.
    void build_and_equilibrate();
    /// The production leg `--stage all` always runs: adaptive, metadynamics or
    /// fixed-length, whichever this Config selects.
    void run_production();
    bool already_equilibrated() const;
    void run_steps(std::initializer_list<Step> steps);
    /// `--stage report`: renders the offline figures for whichever kind of run
    /// this project directory actually holds, detected from the output on disk
    /// when --adaptive/--metadynamics were not repeated on the command line.
    void report_existing_run();

    void ensure_dirs();
    /// Writes a residuetypes.dat next to this run's output that lists
    /// SEP/TPO/PTR as Protein, seeded from the GROMACS installation's own copy.
    void ensure_residuetypes();
    /// Fails fast with a clear message for project-state problems that would
    /// otherwise surface much later as a cryptic martinize2/GROMACS error —
    /// most importantly a missing PDB, a PDB carrying a modified residue
    /// (SEP/TPO/PTR) without --phospho.
    void check_inputs();
    void write_all_mdps();
};

} // namespace cg
