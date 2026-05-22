#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

struct SamplingRules {
    double min_total_us = 0.5;
    int min_cn_state_transitions = 10;
    int min_bidirectional_events = 1;
    double min_cn_range = 20.0;
    double min_rg_global_range = 0.5;
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
    int cn_state_transitions = 0, cn_state_up = 0, cn_state_down = 0;
    int unique_contact_patterns = 0, largest_cluster_unique = 0, largest_cluster_max = 0;
    bool enough_time = false;
    bool enough_cn_transitions = false;
    bool bidirectional_cn_motion = false;
    bool enough_cn_range = false;
    bool enough_rg_range = false;
    bool enough_contact_patterns = false;
    bool enough_cluster_diversity = false;
    bool stop = false;

    std::string to_json(int last_chunk) const;
};

SamplingMetrics evaluate_sampling(const std::vector<std::filesystem::path>& colvar_paths,
                                  int n_prot,
                                  double dt_colvar_ps,
                                  const SamplingRules& rules);

void append_metrics_jsonl(const std::filesystem::path& path,
                          const SamplingMetrics& metrics,
                          int last_chunk);

} // namespace cg
