#pragma once

#include <filesystem>
#include <string>

namespace cg {

struct SamplingMetrics;
struct SamplingRules;

/// What scripts/cv_readiness.py reports, as read back from its key-value output
/// file.
struct CvReadiness {
    bool ok = false;
    bool plateau = false;
    double its1_ps = 0.0;
    double total_time_ps = 0.0;
    std::string note;
};

/// Parses the `key value` lines cv_readiness.py writes.
CvReadiness read_cv_readiness(const std::filesystem::path& path);

/// Records a CvReadiness result in `m` and re-evaluates the stop decision.
void apply_cv_readiness(SamplingMetrics& m, const SamplingRules& rules, const CvReadiness& cv);

} // namespace cg
