#include "ColvarTable.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <stdexcept>

namespace cg {
namespace {


double parse_field(const std::string& token, const std::filesystem::path& path) {
    try {
        return std::stod(token);
    } catch (const std::exception&) {
        throw std::runtime_error("Malformed number \"" + token + "\" in COLVAR file: " + path.string());
    }
}

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    const auto mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    return v[mid];
}

} // namespace

std::vector<double> ColvarTable::times() const {
    if (!has_time) return {};
    std::vector<double> out;
    out.reserve(rows.size());
    for (const auto& row : rows) out.push_back(row[0]);
    return out;
}

double ColvarTable::total_time_ps(double fallback_dt_ps) const {
    if (!has_time || rows.empty()) return static_cast<double>(rows.size()) * fallback_dt_ps;
    const double spacing = dt_ps > 0.0 ? dt_ps : fallback_dt_ps;
    return (t_last_ps - t_first_ps) + spacing;
}

ColvarTable read_colvars(const std::vector<std::filesystem::path>& paths) {
    ColvarTable table;
    bool have_last_time = false;
    double last_time = 0.0;
    double spacing = 0.0;
    double offset = 0.0;

    for (const auto& path : paths) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot read COLVAR file: " + path.string());

        std::vector<std::string> local_fields;
        bool have_fields = false;
        for (std::string line; std::getline(in, line);) {
            line = trim(line);
            if (line.empty()) continue;

            if (starts_with(line, "#!")) {
                const auto tokens = split_ws(line);
                if (tokens.size() >= 3 && tokens[1] == "FIELDS") {
                    local_fields.assign(tokens.begin() + 2, tokens.end());
                    have_fields = true;
                    if (table.fields.empty()) {
                        table.fields = local_fields;
                        for (std::size_t i = 0; i < table.fields.size(); ++i) table.col[table.fields[i]] = i;
                        table.has_time = !table.fields.empty() && table.fields.front() == kTimeField;
                    } else if (local_fields != table.fields) {
                        throw std::runtime_error("COLVAR header mismatch in file: " + path.string());
                    }
                }
                continue;
            }

            if (!have_fields && table.fields.empty())
                throw std::runtime_error("COLVAR file has no '#! FIELDS' header: " + path.string());

            const auto tokens = split_ws(line);
            if (tokens.size() != table.fields.size())
                throw std::runtime_error("COLVAR row has " + std::to_string(tokens.size()) +
                                         " columns, expected " + std::to_string(table.fields.size()) +
                                         " in file: " + path.string());

            std::vector<double> row;
            row.reserve(tokens.size());
            for (const auto& t : tokens) row.push_back(parse_field(t, path));

            if (table.has_time) {
                const double slack = spacing > 0.0 ? spacing * 0.5 : 0.0;
                if (have_last_time && row[0] + offset < last_time - slack) {
                    offset = (last_time + spacing) - row[0];
                }
                row[0] += offset;

                const double t = row[0];
                if (have_last_time && !(t > last_time)) {
                    ++table.dropped_non_monotonic_rows;
                    continue;
                }
                if (have_last_time && spacing <= 0.0) spacing = t - last_time;
                last_time = t;
                have_last_time = true;
            }

            table.rows.push_back(std::move(row));
        }
    }

    if (table.rows.empty()) throw std::runtime_error("No COLVAR data rows found.");

    if (table.has_time) {
        table.t_first_ps = table.rows.front()[0];
        table.t_last_ps = table.rows.back()[0];
        std::vector<double> gaps;
        gaps.reserve(table.rows.size());
        for (std::size_t i = 1; i < table.rows.size(); ++i) {
            const double d = table.rows[i][0] - table.rows[i - 1][0];
            if (std::isfinite(d) && d > 0.0) gaps.push_back(d);
        }
        table.dt_ps = median(std::move(gaps));
    }

    return table;
}

std::vector<double> column_values(const ColvarTable& table, const std::string& name) {
    const auto it = table.col.find(name);
    if (it == table.col.end()) throw std::runtime_error("Required COLVAR column not found: " + name);

    std::vector<double> out;
    out.reserve(table.rows.size());
    for (const auto& row : table.rows) out.push_back(row[it->second]);
    return out;
}

} // namespace cg
