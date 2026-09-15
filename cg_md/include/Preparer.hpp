#pragma once

#include <filesystem>
#include <string>

namespace cg {

struct Config;
class Shell;
class GromacsDriver;

/// Takes an all-atom PDB to a solvated, ionized, energy-minimization-ready
/// Martini coarse-grained system.
class Preparer {
public:
    Preparer(Config& cfg, Shell& sh, GromacsDriver& gmx);

    void coarse_grain();
    void write_packmol_input();
    void run_packmol();
    void clean_packmol_pdb();
    void prepare_solvated_system();
    /// For --phospho runs with --solvation-mode insane: insane's built-in
    /// neutralization does not always account for the extra charge carried by
    /// phosphorylated residues, so a corrective genion pass is run here.
    void neutralize_if_needed();
    /// Rewrites NA/CL residue and atom names in the working .gro file to the
    /// NA+/CL- spelling Martini's ion topology expects.
    void normalize_ions();
    /// Builds em.tpr from the current beads .gro; shared by the "all" and "em"
    /// stages, and by the extra phospho-neutralization pass.
    void prepare_em();

private:
    Config& cfg_;
    Shell& sh_;
    GromacsDriver& gmx_;

    void move_generated_itps();
    /// Cross-checks the bead count martinize2 actually produced for ONE
    /// protomer against --atoms-per-prot, immediately after coarse-graining.
    void verify_bead_count() const;
    /// With --martini-lambda-pw != 1.0, refuses to continue unless the lambda
    /// file rescales EVERY bead type this protomer actually uses.
    void verify_lambda_coverage() const;
    /// Refuses to continue when a protomer is long enough to reach its own
    /// periodic image: shortest box vector < Dmax + 2 * rcoulomb.
    void verify_periodic_margin() const;

    /// The PMF's UPPER_WALLS confines the pair to a finite volume, and the
    /// plateau it produces is the unbound reference state.
    void verify_pmf_wall() const;
    void run_insane();
    void make_compact_box();
    void solvate_compact_box();
    void add_ions_to_solvated_system();
    void patch_topology();

    std::string make_ndx_and_find_water_group(const std::filesystem::path& gro_file,
                                              const std::filesystem::path& ndx_path,
                                              const std::filesystem::path& log_path) const;
    void run_genion(const std::filesystem::path& tpr_path,
                    const std::filesystem::path& ndx_path,
                    const std::string& water_group,
                    bool include_salt);
};

} // namespace cg
