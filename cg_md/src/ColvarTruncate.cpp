#include "ColvarTable.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cg {

namespace {

bool is_data_line(const std::string& line) {
    const auto trimmed = trim(line);
    return !trimmed.empty() && trimmed.front() != '#';
}

} // namespace

std::size_t count_plumed_records(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) return 0;
    std::size_t n = 0;
    for (const auto& line : read_lines(path)) n += is_data_line(line) ? 1 : 0;
    return n;
}

std::size_t truncate_plumed_file_to(const std::filesystem::path& path, std::size_t keep_records) {
    if (!std::filesystem::exists(path)) return 0;

    const auto lines = read_lines(path);
    std::vector<std::string> kept;
    kept.reserve(lines.size());

    std::size_t seen = 0, removed = 0;
    for (const auto& line : lines) {
        if (!is_data_line(line)) { kept.push_back(line); continue; }
        if (seen++ < keep_records) kept.push_back(line);
        else ++removed;
    }

    if (removed) write_lines(path, kept);
    return removed;
}

std::size_t truncate_plumed_file_to_segment(const std::filesystem::path& path,
                                            std::size_t keep_segments) {
    if (!std::filesystem::exists(path)) return 0;

    const auto lines = read_lines(path);
    std::vector<std::string> kept;
    kept.reserve(lines.size());

    bool time_is_first_field = false;
    bool saw_fields = false;
    std::size_t segment = 0, removed = 0;
    bool have_last = false;
    double last_t = 0.0;

    for (const auto& line : lines) {
        const auto trimmed = trim(line);
        if (trimmed.empty() || trimmed.front() == '#') {
            if (starts_with(trimmed, "#!")) {
                const auto tokens = split_ws(trimmed);
                if (tokens.size() >= 3 && tokens[1] == "FIELDS") {
                    time_is_first_field = tokens[2] == kTimeField;
                    saw_fields = true;
                }
            }
            kept.push_back(line);
            continue;
        }

        if (!saw_fields || !time_is_first_field) {
            kept.push_back(line);
            continue;
        }

        const auto tokens = split_ws(trimmed);
        double t = 0.0;
        try {
            t = std::stod(tokens.front());
        } catch (const std::exception&) {
            kept.push_back(line);
            continue;
        }

        if (have_last && t < last_t) ++segment;
        have_last = true;
        last_t = t;

        if (segment < keep_segments) kept.push_back(line);
        else ++removed;
    }

    if (removed) write_lines(path, kept);
    return removed;
}

std::vector<std::pair<double, double>> column_series(const ColvarTable& table,
                                                     const std::string& name,
                                                     double fallback_dt_ps) {
    const auto values = column_values(table, name);
    std::vector<std::pair<double, double>> series;
    series.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double t = table.has_time ? table.rows[i][0] : static_cast<double>(i) * fallback_dt_ps;
        series.emplace_back(t, values[i]);
    }
    return series;
}

} // namespace cg
