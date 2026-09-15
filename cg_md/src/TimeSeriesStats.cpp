#include "TimeSeriesStats.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cg {

std::string TimeSeriesStats::to_json() const {
    std::ostringstream o;
    o << std::setprecision(10)
      << "{\"mean\":" << mean << ",\"stddev\":" << stddev << ",\"drift_per_ns\":" << drift_per_ns << "}";
    return o.str();
}

std::string normalise_legend(const std::string& name) {
    std::string out;
    for (const char c : name) {
        if (c == ' ' || c == '-' || c == '_') continue;
        out += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return out;
}

std::map<std::string, std::vector<std::pair<double, double>>>
read_xvg_by_legend(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read .xvg file: " + path.string());

    std::map<std::size_t, std::string> legend;
    std::vector<std::vector<std::pair<double, double>>> cols;
    std::size_t rows = 0;

    for (std::string line; std::getline(in, line);) {
        std::size_t i = 0;
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        if (i == line.size()) continue;

        if (line[i] == '@') {
            const auto s = line.find(" s", i);
            const auto q1 = line.find('"');
            const auto q2 = line.rfind('"');
            if (s == std::string::npos || q1 == std::string::npos || q2 <= q1) continue;
            if (line.find("legend", s) == std::string::npos) continue;
            try {
                legend[static_cast<std::size_t>(std::stoul(line.substr(s + 2)))] =
                    line.substr(q1 + 1, q2 - q1 - 1);
            } catch (const std::exception&) {  }
            continue;
        }
        if (line[i] == '#') continue;

        std::istringstream ss(line);
        double t = 0.0;
        if (!(ss >> t)) continue;
        std::size_t c = 0;
        for (double v = 0.0; ss >> v; ++c) {
            if (c >= cols.size()) cols.emplace_back();
            cols[c].emplace_back(t, v);
        }
        ++rows;
    }

    if (!rows) throw std::runtime_error("No data rows found in .xvg file: " + path.string());
    if (legend.empty())
        throw std::runtime_error(".xvg file " + path.string() + " carries no '@ sN legend' lines, "
                                 "so its columns cannot be identified by name.");

    std::map<std::string, std::vector<std::pair<double, double>>> out;
    for (const auto& [index, name] : legend)
        if (index < cols.size()) out[normalise_legend(name)] = cols[index];
    return out;
}

TimeSeriesStats compute_time_series_stats(const std::vector<std::pair<double, double>>& series) {
    TimeSeriesStats stats;
    if (series.empty()) return stats;

    double sum_y = 0.0;
    for (const auto& [t, y] : series) { (void)t; sum_y += y; }
    stats.mean = sum_y / static_cast<double>(series.size());

    double sum_sq = 0.0;
    for (const auto& [t, y] : series) { (void)t; sum_sq += (y - stats.mean) * (y - stats.mean); }
    stats.stddev = std::sqrt(sum_sq / static_cast<double>(series.size()));

    if (series.size() < 2) return stats;

    double sum_t = 0.0;
    for (const auto& [t, y] : series) { (void)y; sum_t += t; }
    const double mean_t = sum_t / static_cast<double>(series.size());

    double num = 0.0, den = 0.0;
    for (const auto& [t, y] : series) {
        num += (t - mean_t) * (y - stats.mean);
        den += (t - mean_t) * (t - mean_t);
    }
    if (den > 1e-12) stats.drift_per_ns = (num / den) * 1000.0;
    return stats;
}

} // namespace cg
