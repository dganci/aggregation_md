#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cg {

/// PLUMED's leading time column, as read_colvars() and the truncation helpers
/// in ColvarTruncate.cpp both look for it.
inline constexpr const char* kTimeField = "time";

/// In-memory representation of one or more concatenated PLUMED COLVAR files
/// sharing the same `#! FIELDS` header (as produced across AdaptiveSampler's
/// chunk restarts, MetadynamicsRunner's batch restarts, or a single
/// fixed-length production run).
struct ColvarTable {
    std::vector<std::string> fields;
    std::unordered_map<std::string, std::size_t> col;
    std::vector<std::vector<double>> rows;

    /// True when the files carried PLUMED's automatic leading "time" column
    /// (always the case for real PLUMED output; synthetic test fixtures may
    /// omit it).
    bool has_time = false;

    /// Simulated time of the first / last retained frame, in ps.
    double t_first_ps = 0.0;
    double t_last_ps = 0.0;

    /// Median spacing between consecutive retained frames, in ps (0 when
    /// has_time is false or fewer than two frames were read).
    double dt_ps = 0.0;

    /// Number of rows dropped by the monotonic-time filter in read_colvars() -
    /// i.e. duplicated boundary frames between consecutive batches, or rows
    /// from a re-run segment overlapping one already read.
    std::size_t dropped_non_monotonic_rows = 0;

    /// Values of the leading "time" column, one per row (empty when has_time is
    /// false).
    std::vector<double> times() const;

    /// Total simulated time spanned by the table, in ps.
    double total_time_ps(double fallback_dt_ps) const;
};

/// Reads and concatenates `paths` in order.
ColvarTable read_colvars(const std::vector<std::filesystem::path>& paths);

/// Returns the values of column `name`, one per row of `table`, in row order.
std::vector<double> column_values(const ColvarTable& table, const std::string& name);

/// Like column_values(), but returns (time_ps, value) pairs suitable for
/// compute_time_series_stats().
std::vector<std::pair<double, double>> column_series(const ColvarTable& table,
                                                     const std::string& name,
                                                     double fallback_dt_ps);

/// Number of data (non-comment, non-blank) rows in a PLUMED output file, or 0
/// if it does not exist.
std::size_t count_plumed_records(const std::filesystem::path& path);

/// Rewrites a PLUMED output file (COLVAR, HILLS, ...) in place, keeping every
/// comment/header line and only the first `keep_records` data rows.
std::size_t truncate_plumed_file_to(const std::filesystem::path& path, std::size_t keep_records);

/// Fallback for runs that predate the record-count ledger: keeps the first
/// `keep_segments` segments and drops everything after them, where a new
/// segment begins wherever the leading time column goes BACKWARDS.
std::size_t truncate_plumed_file_to_segment(const std::filesystem::path& path,
                                            std::size_t keep_segments);

} // namespace cg
