#pragma once

#include "TimeSeriesStats.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

class GromacsDriver;

/// Standard, well-established MD "is this run behaving" checks.
struct EnergyDiagnostics {
    TimeSeriesStats temperature_K;
    TimeSeriesStats total_energy_kJ_per_mol;
    TimeSeriesStats pressure_bar;
    TimeSeriesStats density_kg_per_m3;
    bool has_pressure = false;
    bool has_density = false;

    std::string to_json() const;
};

/// Per-protomer radius-of-gyration statistics and the distribution of
/// inter-protomer cluster sizes observed across frames.
struct StructuralDiagnostics {
    std::vector<TimeSeriesStats> per_chain_rg_nm;
    TimeSeriesStats global_rg_nm;
    TimeSeriesStats total_contacts;
    /// Size of the largest oligomer, and the number of distinct oligomers, per
    /// frame.
    TimeSeriesStats largest_cluster;
    TimeSeriesStats n_clusters;

    /// index i (0-based) -> number of (frame, cluster) pairs in which some
    /// oligomer had size i+1.
    std::vector<int> cluster_size_histogram;
    /// index i (0-based) -> number of frames whose *largest* oligomer had size
    /// i+1.
    std::vector<int> largest_cluster_histogram;

    /// Fraction of frames in which each protomer pair was in contact, in the
    /// same i<j order as pair_labels ("cn_i_j").
    std::vector<double> pair_contact_occupancy;
    std::vector<std::string> pair_labels;

    int n_frames = 0;
    double total_time_ps = 0.0;

    std::string to_json() const;
};

struct RunDiagnostics {
    EnergyDiagnostics energy;
    StructuralDiagnostics structure;
    bool has_structure = false;

    std::string to_json(const std::string& segment_label) const;
};

/// Runs `gmx energy` once (see GromacsDriver::energy()) on `edr_path` for
/// Temperature and Total Energy, plus Pressure/Density when requested, writing
/// scratch .xvg files under `scratch_dir` and reducing each series via
/// compute_time_series_stats().
EnergyDiagnostics compute_energy_diagnostics(const GromacsDriver& gmx,
                                             const std::filesystem::path& edr_path,
                                             const std::filesystem::path& scratch_dir,
                                             bool include_pressure,
                                             bool include_density);

/// Reduces the concatenation of `colvar_paths` (see ColvarTable.hpp) to
/// per-chain Rg statistics (`dt_colvar_ps` converts row index to simulated
/// time, so per-chain Rg drift is expressed per ns like every other
/// TimeSeriesStats) and a histogram of connected-component ("oligomer") sizes
/// observed across all frames, using the same pairwise-contact definition
/// (cn_i_j >= pair_contact_threshold) as SamplingMonitor's evaluate_sampling().
StructuralDiagnostics compute_structural_diagnostics(const std::vector<std::filesystem::path>& colvar_paths,
                                                     int n_prot,
                                                     double dt_colvar_ps,
                                                     double pair_contact_threshold);

/// Appends one JSON-Lines record (tagging it with `segment_label`, e.g. "nvt",
/// "npt", "chunk_003", "metad") to `path`, creating parent directories as
/// needed.
void append_diagnostics_jsonl(const std::filesystem::path& path,
                              const RunDiagnostics& diag,
                              const std::string& segment_label);

} // namespace cg
