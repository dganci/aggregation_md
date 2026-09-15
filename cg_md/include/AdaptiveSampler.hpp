#pragma once

#include "SamplingMonitor.hpp"

#include <filesystem>
#include <vector>

namespace cg {

struct Config;
class Shell;
class GromacsDriver;

/// Runs production MD in fixed-length chunks, evaluating aggregation-relevant
/// sampling metrics (see SamplingMonitor) after every chunk and stopping as
/// soon as the configured convergence criteria are met, or once
/// adaptive_max_total_us is reached.
class AdaptiveSampler {
public:
    AdaptiveSampler(Config& cfg, Shell& sh, GromacsDriver& gmx);
    void run();

private:
    Config& cfg_;
    Shell& sh_;
    GromacsDriver& gmx_;

    SamplingRules rules() const;
    /// Runs scripts/cv_readiness.py over `colvars` and reads its result; a
    /// missing script or a failed run comes back with ok = false.
    CvReadiness cv_readiness(const std::vector<std::filesystem::path>& colvars) const;
    std::filesystem::path chunk_base(int chunk) const;
    std::filesystem::path chunk_tpr(int chunk) const;
    std::filesystem::path chunk_colvar(int chunk) const;
    std::filesystem::path chunk_plumed(int chunk) const;
    bool chunk_completed(int chunk) const;
    std::vector<int> completed_chunks() const;
    std::vector<std::filesystem::path> completed_colvars() const;
    /// Renders the offline figures/summary for the completed run.
    void report() const;

    /// Joins the completed chunk trajectories into one continuous .xtc and
    /// un-breaks the molecules, so trajectory-based analysis (residue contact
    /// maps, shape descriptors, clustering) can run on an adaptive run the same
    /// way it can on a fixed-length one.
    void center_chunks() const;
};

} // namespace cg
