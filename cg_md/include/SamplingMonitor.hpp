#pragma once

#include "CvReadiness.hpp"
#include "SamplingEstimators.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

/// Thresholds an adaptive-sampling chunk sequence must clear before
/// AdaptiveSampler stops early.
struct SamplingRules {
    double min_total_us = 0.5;
    int min_cn_state_transitions = 10;
    int min_bidirectional_events = 1;
    double min_cn_range = 20.0;
    double min_rg_global_range = 0.5;
    /// Contact patterns are counted twice: over the first half of the
    /// trajectory and over all of it.
    double max_pattern_growth_ratio = 1.10;
    /// Jensen-Shannon divergence, in bits, between the contact-pattern
    /// frequencies of the two halves.
    double max_pattern_jsd_bits = 0.05;
    /// Independent samples of cn_total, from the block-averaging plateau.
    int min_effective_samples = 20;
    /// Growth and shrink events of the largest oligomer, each way, after the
    /// residence filter below.
    int min_assembly_events = 20;
    /// How long a change in the largest oligomer must persist to count as an
    /// event, in ps.
    double event_residence_ps = 100.0;
    /// Sampled time over the slowest implied timescale; 0 disables the gate.
    double min_time_over_its = 10.0;
    int min_unique_contact_patterns = 10;
    int min_largest_cluster_unique = 3;
    double pair_contact_threshold = 1.0;
    int smooth_window_frames = 50;
    int min_residence_frames = 5;
    int pattern_downsample = 10;
};

struct SamplingMetrics {
    double total_time_us = 0.0;
    int n_frames = 0;
    double cn_min = 0.0, cn_max = 0.0, cn_range = 0.0;
    double rg_global_min = 0.0, rg_global_max = 0.0, rg_global_range = 0.0;
    int patterns_first_half = 0;
    double pattern_growth_ratio = 0.0;
    /// See SamplingRules::max_pattern_jsd_bits.
    double pattern_jsd_halves_bits = 0.0;
    /// Contact patterns holding at least 1% of the sampled frames.
    int populated_contact_patterns = 0;
    int cn_state_transitions = 0, cn_state_up = 0, cn_state_down = 0;
    int unique_contact_patterns = 0, largest_cluster_unique = 0, largest_cluster_max = 0;

    /// Mean size of the largest oligomer over the sampled frames, and how often
    /// it grew / shrank between consecutive sampled frames, with no residence
    /// filter: every touch counts.
    double largest_cluster_mean = 0.0;
    int lcc_growth_events = 0;
    int lcc_shrink_events = 0;

    /// The same growth / shrink counts after SamplingRules::event_residence_ps:
    /// a new largest-oligomer size has to persist that long before the change
    /// is an event.
    int assembly_growth_events = 0;
    int assembly_shrink_events = 0;

    /// Independent samples in the cn_total series - see
    /// SamplingRules::min_effective_samples.
    double cn_effective_samples = 0.0;

    /// From scripts/cv_readiness.py, when AdaptiveSampler could run it.
    double cv_its1_ps = 0.0;
    bool cv_plateau_found = false;
    double cv_time_over_its = 0.0;

    /// Frames discarded as duplicated batch-boundary repeats when the chunk
    /// COLVARs were concatenated (see ColvarTable).
    int duplicate_frames_dropped = 0;

    bool enough_time = false;
    bool enough_cn_transitions = false;
    bool bidirectional_cn_motion = false;
    bool enough_cn_range = false;
    bool enough_rg_range = false;
    bool enough_contact_patterns = false;
    bool exploration_saturated = false;
    bool pattern_distribution_settled = false;
    bool enough_independent_samples = false;
    bool enough_assembly_events = false;
    bool enough_for_cv_training = false;
    bool enough_cluster_diversity = false;
    bool stop = false;

    std::string to_json(int last_chunk) const;
};

/// Reads the concatenation of `colvar_paths` (PLUMED COLVAR files, all sharing
/// the same #! FIELDS header) and computes the aggregation-sampling metrics
/// used to decide whether adaptive production has run long enough.
SamplingMetrics evaluate_sampling(const std::vector<std::filesystem::path>& colvar_paths,
                                  int n_prot,
                                  double dt_colvar_ps,
                                  const SamplingRules& rules);

/// Sets every check_* flag and `stop` from the numeric fields of `m`.
void decide_stop(SamplingMetrics& m, const SamplingRules& rules);

/// Appends one line of `metrics.to_json(last_chunk)` to `path` (JSON Lines).
void append_metrics_jsonl(const std::filesystem::path& path,
                          const SamplingMetrics& metrics,
                          int last_chunk);

} // namespace cg
