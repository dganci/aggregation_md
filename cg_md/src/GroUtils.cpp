#include "GroUtils.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace cg {
namespace {

std::string field5(const std::string& value) {
    std::ostringstream out;
    out << std::setw(5) << value;
    return out.str();
}

std::string ion_name(const std::string& residue, const std::string& atom) {
    if (residue == "NA" || residue == "NA+" || atom == "NA" || atom == "NA+") return "NA+";
    if (residue == "CL" || residue == "CL-" || atom == "CL" || atom == "CL-") return "CL-";
    return {};
}

} // namespace

GroNormalizationStats normalize_gro_ion_names(const std::filesystem::path& gro_path) {
    auto lines = read_lines(gro_path);
    if (lines.size() < 3) throw std::runtime_error("Invalid GRO file: " + gro_path.string());

    const auto natoms = static_cast<std::size_t>(std::stoul(trim(lines[1])));
    const auto first = std::size_t{2};
    const auto last = first + natoms;
    if (lines.size() <= last) throw std::runtime_error("GRO atom count exceeds file length: " + gro_path.string());

    GroNormalizationStats stats;
    for (std::size_t i = first; i < last; ++i) {
        auto& line = lines[i];
        if (line.size() < 20) continue;

        const auto residue = trim(line.substr(5, 5));
        const auto atom = trim(line.substr(10, 5));
        const auto name = ion_name(residue, atom);
        if (name.empty()) continue;

        if (atom != name) {
            line.replace(10, 5, field5(name));
            ++stats.atom_names;
        }
        if ((residue == "NA" || residue == "NA+" || residue == "CL" || residue == "CL-") && residue != name) {
            line.replace(5, 5, field5(name));
            ++stats.residue_names;
        }
    }

    if (stats.total()) write_lines(gro_path, lines);
    return stats;
}

std::array<std::array<double, 3>, 3> read_gro_box_vectors(const std::filesystem::path& gro_path) {
    const auto lines = read_lines(gro_path);
    for (auto it = lines.rbegin(); it != lines.rend(); ++it) {
        std::istringstream in(*it);
        double v[9] = {};
        if (!(in >> v[0] >> v[1] >> v[2])) continue;
        for (int i = 3; i < 9; ++i) {
            if (!(in >> v[i])) { v[i] = 0.0; }
        }
        return {{{v[0], v[3], v[4]}, {v[5], v[1], v[6]}, {v[7], v[8], v[2]}}};
    }
    throw std::runtime_error("Could not read a box line from " + gro_path.string());
}

double min_image_distance_nm(const std::array<std::array<double, 3>, 3>& box) {
    double best = std::numeric_limits<double>::max();
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            for (int k = -2; k <= 2; ++k) {
                if (i == 0 && j == 0 && k == 0) continue;
                double len2 = 0.0;
                for (int d = 0; d < 3; ++d) {
                    const double c = i * box[0][d] + j * box[1][d] + k * box[2][d];
                    len2 += c * c;
                }
                best = std::min(best, len2);
            }
        }
    }
    return std::sqrt(best);
}

bool periodic_margin_ok(double dmax_nm, double image_nm, double rcut_nm) {
    return image_nm >= dmax_nm + 2.0 * rcut_nm;
}

} // namespace cg
