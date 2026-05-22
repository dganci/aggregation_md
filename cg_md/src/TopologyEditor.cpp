#include "TopologyEditor.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace cg {
namespace {

std::vector<std::string> martini_includes(const std::string& itp) {
    return {
        "#include \"../../martini_v300/martini_v300/martini_v3.0.0.itp\"",
        "#include \"../../martini_v300/martini_v300/martini_v3.0.0_ions_v1.itp\"",
        "#include \"../../martini_v300/martini_v300/martini_v3.0.0_solvents_v1.itp\"",
        "#include \"" + itp + "\""
    };
}

bool is_replaced_include(const std::string& line, const std::vector<std::string>& includes) {
    if (line == "#include \"martini.itp\"") return true;
    return std::find(includes.begin(), includes.end(), line) != includes.end();
}

std::string molecule_line(const std::string& molecule, const std::string& itp, int n_prot) {
    const auto name = itp.empty() ? molecule : basename_without_ext(itp);
    std::ostringstream out;
    out << name << (name.size() < 20 ? std::string(20 - name.size(), ' ') : " ") << n_prot;
    return out.str();
}

} // namespace

std::string find_first_itp(const std::filesystem::path& topology_dir, const std::string& prefix) {
    for (const auto& entry : std::filesystem::directory_iterator(topology_dir)) {
        const auto name = entry.path().filename().string();
        if (entry.is_regular_file() && starts_with(name, prefix) && ends_with(name, ".itp")) return name;
    }
    throw std::runtime_error("No ITP found in " + topology_dir.string() + " with prefix " + prefix);
}

void patch_martini_topology(const std::filesystem::path& topology_path,
                            const std::string& molecule_name,
                            int n_prot,
                            const std::string& itp_filename) {
    const auto includes = martini_includes(itp_filename);
    const auto lines = read_lines(topology_path);
    std::vector<std::string> out = includes;

    bool in_molecules = false;
    bool replaced = false;
    for (const auto& line : lines) {
        const auto t = trim(line);
        if (is_replaced_include(t, includes)) continue;

        if (starts_with(t, "[")) in_molecules = starts_with(t, "[ molecules ]");
        out.push_back(line);

        if (!in_molecules || replaced || t.empty() || starts_with(t, ";") || starts_with(t, "[")) continue;
        out.back() = molecule_line(molecule_name, itp_filename, n_prot);
        replaced = true;
    }

    write_lines(topology_path, out);
}

} // namespace cg
