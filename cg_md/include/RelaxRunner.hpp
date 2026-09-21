#pragma once

#include <filesystem>

namespace cg {

struct Config;
class Shell;
class GromacsDriver;

/// Produces a relaxed, representative starting conformation for one protomer,
/// in its own directory, and writes it where Packmol will pick it up.
class RelaxRunner {
public:
    RelaxRunner(Config& cfg, Shell& sh, GromacsDriver& gmx);

    /// Runs the whole relaxation and installs the representative structure.
    void run();

private:
    Config& cfg_;
    Shell& sh_;
    GromacsDriver& gmx_;

    std::filesystem::path dir() const;
    /// The single-copy topology, written beside the project topology so its
    /// relative `#include` of the protomer .itp resolves.
    std::filesystem::path topology() const;
    /// Single-copy solvated, ionized system built from the coarse-grained
    /// protomer, with its own one-molecule topology.
    void build_system();
    /// Writes runs/relax_<protomer>/index.ndx with a "Solvent_and_ions" group
    /// and returns its path (empty when the thermostat does not need one).
    ///
    /// --thermostat protein-solvent makes MdpWriter emit
    /// `tc-grps = Protein Solvent_and_ions`, and Solvent_and_ions is NOT a
    /// group GROMACS derives on its own - it only exists in an index file.
    /// The relaxation's grompp passed no -n, so every entry died on the NVT
    /// step with "group not found", before a single ns of anything. It did not
    /// show up earlier because the manifest then used --thermostat legacy,
    /// whose `Protein W ION` groups ARE derived from residue names.
    ///
    /// IndexBuilder::generate_index() cannot be reused here: it lays out
    /// Config::n_prot protomers back-to-back, and the relaxation box holds
    /// exactly one.
    std::filesystem::path write_relax_index() const;
    void equilibrate();
    void produce();
    /// Clusters the protein-only trajectory and installs the centroid of the
    /// most populated cluster as the new Packmol input structure.
    void select_representative();
    /// Writes --n-prot distinct starting conformations, one per copy, as evenly
    /// spaced frames of the relaxation trajectory.
    void write_distinct_copies();
};

} // namespace cg
