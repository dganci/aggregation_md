#pragma once

#include <filesystem>
#include <map>
#include <utility>
#include <vector>

namespace cg {

/// Mean, (population) standard deviation, and linear drift of a scalar time
/// series - the standard "is this run behaving" reduction applied to GROMACS
/// energy-file observables (temperature, pressure, density, total energy).
struct TimeSeriesStats {
    double mean = 0.0;
    double stddev = 0.0;
    /// Ordinary-least-squares slope of value against time, in value-units per
    /// nanosecond (time in the input series is expected in ps, GROMACS's own
    /// .xvg convention).
    double drift_per_ns = 0.0;

    std::string to_json() const;
};

/// Reads a multi-column .xvg and keys each series by its "@ sN legend" name,
/// normalised (lowercased, spaces/hyphens/underscores removed) so that the name
/// gmx ENERGY PRINTS matches the name it was SELECTED by.
std::map<std::string, std::vector<std::pair<double, double>>>
read_xvg_by_legend(const std::filesystem::path& path);

/// The key read_xvg_by_legend() files a series under, for a gmx energy term.
std::string normalise_legend(const std::string& name);

/// Reduces `series` (pairs of (time_ps, value)) to mean/stddev/drift.
TimeSeriesStats compute_time_series_stats(const std::vector<std::pair<double, double>>& series);

} // namespace cg
