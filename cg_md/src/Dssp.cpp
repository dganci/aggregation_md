#include "Dssp.hpp"
#include "Config.hpp"
#include "Shell.hpp"

namespace cg {
namespace {

char martinize_ss(char c) {
    switch (c) {
        case 'H': case 'G': case 'I': return 'H';
        case 'E': case 'B': return 'E';
        default: return 'C';
    }
}

} // namespace

std::string run_mkdssp(const Shell& sh, const Config& cfg) {
    return sh.runCapture({cfg.mkdssp, cfg.pdbPath().string()}, {}, true).output;
}

std::string dssp_to_martinize_ss(const std::string& dssp_output, int seq_length) {
    bool in_table = false;
    std::string ss;
    ss.reserve(static_cast<std::size_t>(seq_length));

    for (std::size_t first = 0; first < dssp_output.size();) {
        auto last = dssp_output.find('\n', first);
        if (last == std::string::npos) last = dssp_output.size();
        const auto line = dssp_output.substr(first, last - first);

        if (line.find("#  RESIDUE") != std::string::npos || line.find("# RESIDUE") != std::string::npos) {
            in_table = true;
        } else if (in_table && line.size() > 16) {
            const char aa = line.size() > 13 ? line[13] : ' ';
            if (aa != '!' && aa != '*') ss.push_back(martinize_ss(line[16]));
        }
        first = last + 1;
    }

    if (static_cast<int>(ss.size()) < seq_length) return {};
    ss.resize(static_cast<std::size_t>(seq_length));
    return ss;
}

} // namespace cg
