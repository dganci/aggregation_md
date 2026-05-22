#include "GroUtils.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <iomanip>
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

} // namespace cg
