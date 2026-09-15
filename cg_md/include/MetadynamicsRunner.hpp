#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

struct Config;
class Shell;
class GromacsDriver;

/// Metadynamics production driven by a pre-trained PyTorch collective variable
/// (see PlumedWriter::write_metad_plumed_dat), optionally replicated across
/// multiple walkers via `mdrun -multidir`, and optionally split into
/// consecutive batches.
class MetadynamicsRunner {
public:
    MetadynamicsRunner(Config& cfg, Shell& sh, GromacsDriver& gmx);

    /// Prepares inputs, runs every outstanding batch, and finalizes (sum_hills
    /// + centering + report).
    void run();
    /// Sum_hills + trajectory centering + report only; used by the standalone
    /// "fes" stage to re-finalize an already-completed metadynamics run.
    void finalize();

private:
    Config& cfg_;
    Shell& sh_;
    GromacsDriver& gmx_;

    /// Directories mdrun runs in *according to the current options*: the
    /// metadynamics directory itself for a single walker, one walkerN
    /// subdirectory each otherwise.
    std::vector<std::filesystem::path> run_dirs() const;

    /// The same directories, but discovered from disk.
    std::vector<std::filesystem::path> existing_run_dirs() const;

    std::string batch_name(int batch) const;
    /// True when `batch` finished in every walker directory and left behind the
    /// .gro/.cpt the next batch needs to continue from.
    bool batch_completed(int batch) const;
    int first_pending_batch() const;

    /// Loads sigma / grid bounds / the trained feature list from the CV
    /// parameter pickle, and cross-checks them against the CLI options.
    void load_cv_parameters();
    void prepare_common();
    void prepare_batch(int batch);
    void run_batch(int batch);

    /// The PLUMED output files a batch appends to (HILLS, the CV COLVAR, the
    /// monitoring COLVAR).
    std::vector<std::string> plumed_outputs() const;
    /// Records how long each of them was when `batch` finished, so a later
    /// re-run of the *next* batch can rewind them exactly.
    void record_plumed_lengths(const std::filesystem::path& dir, int batch) const;
    /// Drops whatever an interrupted earlier attempt at `batch` appended.
    void rewind_plumed_outputs(const std::filesystem::path& dir, int batch) const;

    /// Checks PLUMED's own log for how many walkers it joined, after the first
    /// batch.
    void check_walker_sharing() const;

    /// PLUMED backup files per run directory, in run_dirs() order.
    std::vector<std::size_t> plumed_backup_counts() const;

    /// Throws if batch k > 0 restarted from a flat bias instead of resuming the
    /// accumulated one - a new bck.* file, or a PLUMED output that got shorter.
    void check_bias_continuity(int batch, const std::vector<std::size_t>& before) const;

    /// Purely observational energy + structural diagnostics for the completed
    /// metadynamics run (one record per walker).
    void analyze();
};

} // namespace cg
