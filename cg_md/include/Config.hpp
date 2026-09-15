#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace cg {

/// Placeholder value of Config::elastic_units, refused by validate().
inline const std::string kDefaultElasticUnits = "90:102";

/// Every tunable parameter of the cg_md pipeline, plus the derived paths used
/// to lay out a project directory.
struct Config {
    std::filesystem::path project_dir = ".";

    std::string gmx = "gmx_mpi";
    std::string martinize = "martinize2";
    std::string insane = "insane";
    std::string packmol = "packmol";
    std::string plumed = "plumed";
    std::string mpirun = "mpirun";
    std::string python = "python3";

    int n_prot = 10;
    std::string protomer_name = "desmin_head";
    int seq_length = 108;
    /// Set by --phospho.
    bool is_phospho = false;
    int atoms_per_prot = 237;

    /// Explicit per-residue secondary-structure string passed to martinize2 via
    /// `-ss` when is_phospho is false.
    std::string ss_string = "";

    /// Let martinize2 derive secondary structure from the input structure
    /// (`-dssp`) instead of using --ss-string.
    bool use_dssp = false;

    bool relax = false;
    /// Length of the single-protomer relaxation run, in us.
    double relax_us = 1.0;
    /// RMSD cut-off (nm) for the gromos clustering that picks the
    /// representative conformation.
    double relax_cluster_cutoff_nm = 0.25;
    /// Fraction of the relaxation trajectory discarded as equilibration before
    /// clustering.
    double relax_discard_frac = 0.1;
    /// How many frames the relaxation trajectory should hold, whatever
    /// --relax-us and --nst-xtc are.
    int relax_frames = 2000;

    /// Extra martinize2/vermouth force-field and mapping directories, passed
    /// via `-ff-dir`/`-map-dir`.
    std::filesystem::path martinize_ff_dir = "force_fields";
    std::filesystem::path martinize_map_dir = "mappings";

    /// Passed through verbatim as martinize2 `-modify RESNAME-RESID` (e.g. to
    /// apply a single one-off residue modification outside the --phospho path).
    std::string modify = "";

    /// Comma-separated PDB chain ids to merge into ONE coarse-grained molecule
    /// (martinize2's `-merge`), e.g. "A,B".
    std::string merge_chains = "";

    double packmol_box_A = 145.0;
    double packmol_cluster_radius_A = 120.0;
    double packmol_tolerance_A = 2.5;
    double packmol_radius_A = 5.0;

    std::string solvation_mode = "gromacs";
    std::string box_type = "dodecahedron";
    double box_margin_nm = 1.5;
    double solvate_radius_nm = 0.21;
    std::filesystem::path martini_water_gro = "system/water.gro";

    /// Scaling factor on the protein-water interactions of the Martini force
    /// field (lambda_PW).
    double martini_lambda_pw = 1.0;

    /// Give each of the --n-prot copies a DIFFERENT starting conformation,
    /// taken as evenly spaced frames from the relaxation trajectory, instead of
    /// replicating one representative --n-prot times.
    bool distinct_copies = false;
    bool add_ions = true;

    double temperature_K = 300.0;
    double salt_M = 0.15;

    double epsilon_r = 15.0;
    /// Non-bonded cut-offs (nm).
    double rcoulomb_nm = 1.1;
    double rvdw_nm = 1.1;
    int nstlist = 20;

    /// Verlet buffer.
    double verlet_buffer_tolerance = -1.0;
    double rlist_nm = 1.35;

    /// Master seed for every stochastic step in the pipeline.
    int seed = -1;

    double em_emstep = 0.01;
    int em_nsteps = 50000;

    double nvt_dt_ps = 0.01;
    std::int64_t nvt_nsteps = 100000;

    double npt_dt_ps = 0.01;
    std::int64_t npt_nsteps = 1000000;

    double md_dt_ps = 0.01;
    double md_total_us = 10.0;

    int plumed_stride = 100;
    double contact_r0_nm = 0.55;
    /// RATIONAL switching-function exponents for the pairwise inter-protomer
    /// COORDINATION descriptors.
    int contact_nn = 6;
    int contact_mm = 12;
    /// Trajectory write-out interval, in steps, for every dynamics stage
    /// (nstxout-compressed / nstenergy / nstlog).
    int nst_xtc = 1000;
    int nst_energy = 1000;
    int nst_log = 1000;

    /// Precision of the compressed trajectory, in 1/nm.
    int xtc_precision = 100;
    int ntomp = 0;

    /// GROMACS thread pinning (`mdrun -pin/-pinoffset/-pinstride`).
    std::string pin = "auto";
    int pinoffset = -1;
    int pinstride = -1;

    bool metadynamics = false;
    int metad_walkers = 1;
    int metad_nodes = 3;
    int metad_pace = 500;
    int metad_print_stride = 100;
    double metad_height = 0.8;
    double metad_biasfactor = 10.0;
    std::string metad_sigma;
    std::filesystem::path metad_model;
    std::filesystem::path metad_cv_params;
    bool metad_sum_hills = true;
    bool metad_center = true;

    /// Total metadynamics production time.
    double metad_total_us = 0.0;
    /// Length of one metadynamics batch.
    double metad_chunk_us = 0.0;

    /// Feed the CV model permutation-invariant (sorted) descriptors rather than
    /// the raw index-ordered ones.
    bool metad_permutation_invariant = true;

    /// Ordered list of COLVAR column names fed to the PYTORCH_MODEL CV.
    std::string metad_feature_cols;

    /// Bias grid.
    bool metad_grid = true;
    std::string metad_grid_min;
    std::string metad_grid_max;
    int metad_grid_bin = 200;
    std::string metad_grid_file = "grid.dat";
    int metad_grid_wstride = 0;

    /// Restraining walls at the grid boundary.
    double metad_wall_kt = 50.0;
    double metad_wall_margin_frac = 0.05;

    /// c(t) reweighting factor (METAD's rbias/rct components).
    bool metad_calc_rct = true;
    int metad_rct_ustride = 10;
    /// Accumulated bias work - a cheap well-tempered convergence diagnostic.
    bool metad_calc_work = true;

    /// Filename of the auxiliary, monitoring-only PRINT emitted alongside the
    /// CV COLVAR (raw physical descriptors + rbias, for offline reweighting).
    std::string metad_monitor_file = "COLVAR_monitor";

    /// Emit fes_<n>.dat every `metad_fes_stride` hills in the finalization
    /// sum_hills pass, so the FES convergence with simulated time can be
    /// plotted.
    int metad_fes_stride = 0;

    bool pmf = false;
    /// Flat-bottom upper wall on d_1_2, in nm.
    double pmf_wall_nm = 6.0;
    double pmf_wall_kappa = 2000.0;
    /// Gaussian width along d_1_2, in nm.
    double pmf_sigma_nm = 0.05;
    /// Bound-state cut-off used when integrating the association constant.
    double pmf_bound_cutoff_nm = 0.0;
    /// Fraction of the wall distance beyond which the PMF is taken to have
    /// plateaued; the reference (unbound) level is averaged over
    /// [pmf_plateau_frac * wall, 0.95 * wall], i.e. inside the wall but past
    /// the interaction range.
    double pmf_plateau_frac = 0.75;

    bool adaptive = false;
    double adaptive_chunk_us = 0.25;
    double adaptive_max_total_us = 5.0;
    double adaptive_min_total_us = 0.5;
    int adaptive_min_cn_state_transitions = 10;
    int adaptive_min_bidirectional_events = 1;
    double adaptive_min_cn_range = 20.0;
    double adaptive_min_rg_global_range = 0.5;
    int adaptive_min_unique_contact_patterns = 10;
    int adaptive_min_largest_cluster_unique = 3;
    double adaptive_pair_contact_threshold = 1.0;
    /// The convergence gates, as opposed to the exploration gates above.
    double adaptive_max_pattern_growth = 1.10;
    double adaptive_max_pattern_jsd = 0.05;
    int adaptive_min_effective_samples = 20;
    int adaptive_min_assembly_events = 20;
    double adaptive_event_residence_ps = 100.0;
    double adaptive_min_time_over_its = 10.0;
    int adaptive_smooth_window_frames = 50;
    int adaptive_min_residence_frames = 5;
    int adaptive_pattern_downsample = 10;

    /// Keep md_all.xtc and md_all_whole.xtc, the two intermediate stages of the
    /// concatenate -> whole -> centre chain AdaptiveSampler runs at the end of
    /// a chunked production.
    bool keep_intermediate_trajectories = false;

    bool report = true;
    /// Contact cut-off (nm) and frame stride for the residue-residue contact
    /// map (scripts/contact_map.py).
    double contact_map_cutoff_nm = 0.6;
    int contact_map_stride = 1;
    std::filesystem::path report_script;

    std::string thermostat_mode = "system";
    /// Residue ranges that each get their own elastic network, in the input
    /// PDB's own numbering (martinize2 -eunit, comma separated, ends included).
    std::string elastic_units = kDefaultElasticUnits;
    bool use_elastic = true;
    bool neutralize_if_phospho = true;

    std::string run_tag = "";

    bool dry_run = false;
    std::string stage = "all";

    std::string systemName() const;

    /// The parameters that define WHAT this trajectory is, as canonical "key =
    /// value" lines.
    std::string identityString() const;

    /// --run-tag when given, else the first 8 hex digits of a 64-bit FNV-1a
    /// hash of identityString().
    std::string runId() const;

    /// The identity is FROZEN at the first call to runId(), and fromArgs()
    /// calls it right after validate(): every Config built from the command
    /// line is therefore born with its id already fixed.
    mutable std::string frozen_run_id;

    /// systemName() + "_" + runId(), the leaf of runDir().
    std::string runName() const;
    std::int64_t mdNsteps() const;
    std::int64_t adaptiveChunkNsteps() const;
    int adaptiveMaxChunks() const;

    /// Total metadynamics length in us (metad_total_us, falling back to
    /// md_total_us) and the corresponding step count.
    double metadTotalUs() const;
    std::int64_t metadNsteps() const;
    /// Steps per metadynamics batch, and how many batches cover metadNsteps().
    std::int64_t metadChunkNsteps() const;
    int metadNumChunks() const;
    /// The length actually simulated: metadNumChunks() * metadChunkNsteps(),
    /// which is >= metadNsteps() and equal to it whenever the total is a whole
    /// multiple of the chunk.
    std::int64_t effectiveMetadNsteps() const;

    std::filesystem::path runDir() const;
    /// Intermediate structures for this run (coarse-grained, relaxed, packed,
    /// solvated) and its topology.
    std::filesystem::path prepDir() const;
    std::filesystem::path pdbPath() const;
    std::filesystem::path cgPath() const;
    /// Where --stage relax writes the representative conformation.
    std::filesystem::path cgRelaxedPath() const;
    /// The structure Packmol actually replicates: the relaxed representative if
    /// --stage relax has produced one, else the raw coarse-grained input.
    std::filesystem::path packmolSourcePath() const;
    /// Path of the i-th distinct starting conformer, for --distinct-copies.
    std::filesystem::path cgRelaxedCopyPath(int i) const;
    /// One structure per copy when --distinct-copies produced them, otherwise a
    /// single-element list holding packmolSourcePath().
    std::vector<std::filesystem::path> packmolSourcePaths() const;
    std::filesystem::path topologyPath() const;
    std::filesystem::path packRawPath() const;
    std::filesystem::path packCleanPath() const;
    std::filesystem::path packInputPath() const;
    std::filesystem::path boxedPath() const;
    std::filesystem::path solvatedPath() const;
    std::filesystem::path ionsTprPath() const;
    std::filesystem::path martiniWaterPath() const;
    /// Name of the main Martini .itp the topology includes; depends on
    /// martini_lambda_pw.
    std::string martiniMainItp() const;
    std::filesystem::path martinizeFfDirPath() const;
    std::filesystem::path martinizeMapDirPath() const;
    std::filesystem::path beadsPath() const;
    std::filesystem::path systemDir() const;
    std::filesystem::path resultDir() const;
    std::filesystem::path metadDir() const;
    std::filesystem::path metadModelPath() const;
    std::filesystem::path metadCvParamsPath() const;

    /// Parses argv into a Config and calls validate() before returning.
    static Config fromArgs(int argc, char** argv);
    /// Throws std::runtime_error on any internally-inconsistent option value.
    void validate() const;
    void print() const;
};

} // namespace cg
