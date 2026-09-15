#include "GromacsDriver.hpp"
#include "Config.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <filesystem>
#include <sstream>
#include <vector>

namespace cg {

void append_pinning(const Config& cfg, std::vector<std::string>& cmd) {
    if (cfg.pin == "auto") return;
    cmd.insert(cmd.end(), {"-pin", cfg.pin});
    if (cfg.pin != "on") return;
    if (cfg.pinoffset >= 0) cmd.insert(cmd.end(), {"-pinoffset", std::to_string(cfg.pinoffset)});
    if (cfg.pinstride >= 0) cmd.insert(cmd.end(), {"-pinstride", std::to_string(cfg.pinstride)});
}

namespace {

std::filesystem::path mdout_for(const std::filesystem::path& tpr) {
    return tpr.parent_path() / (tpr.stem().string() + "_mdout.mdp");
}

} // namespace

GromacsDriver::GromacsDriver(Config& cfg, Shell& sh) : cfg_(cfg), sh_(sh) {}

void GromacsDriver::grompp(const std::filesystem::path& mdp,
                           const std::filesystem::path& coordinates,
                           const std::filesystem::path& output,
                           const std::filesystem::path& checkpoint,
                           const std::filesystem::path& restraints,
                           const std::filesystem::path& index) const {
    std::vector<std::string> cmd = {cfg_.gmx, "grompp", "-f", mdp.string(), "-c", coordinates.string(),
                                    "-p", cfg_.topologyPath().string(), "-o", output.string(),
                                    "-po", mdout_for(output).string()};
    if (!checkpoint.empty()) cmd.insert(cmd.end(), {"-t", checkpoint.string()});
    if (!restraints.empty()) cmd.insert(cmd.end(), {"-r", restraints.string()});
    if (!index.empty()) cmd.insert(cmd.end(), {"-n", index.string()});
    sh_.run(cmd);
}

void GromacsDriver::mdrun(const std::filesystem::path& deffnm,
                          const std::filesystem::path& tpr,
                          const std::filesystem::path& plumed,
                          int ntomp_override) const {
    std::vector<std::string> cmd = {cfg_.gmx, "mdrun", "-v", "-deffnm", deffnm.string()};
    if (!tpr.empty()) cmd.insert(cmd.end(), {"-s", tpr.string()});
    if (!plumed.empty()) cmd.insert(cmd.end(), {"-plumed", plumed.string()});

    const int ntomp = ntomp_override > 0 ? ntomp_override : cfg_.ntomp;
    if (ntomp > 0) cmd.insert(cmd.end(), {"-ntomp", std::to_string(ntomp)});
    append_pinning(cfg_, cmd);
    sh_.run(cmd);
}

void GromacsDriver::energy(const std::filesystem::path& edr,
                           const std::filesystem::path& xvg,
                           const std::vector<std::string>& observables) const {
    std::string selection;
    for (const auto& obs : observables) selection += obs + "\\n";
    sh_.runShell("printf '" + selection + "0\\n' | " + shell_quote(cfg_.gmx) +
                 " energy -f " + shell_quote(edr.string()) +
                 " -o " + shell_quote(xvg.string()));
}

void GromacsDriver::energy(const std::filesystem::path& edr,
                           const std::filesystem::path& xvg,
                           const std::string& observable) const {
    energy(edr, xvg, std::vector<std::string>{observable});
}

void GromacsDriver::trjconv(const std::filesystem::path& tpr,
                            const std::filesystem::path& in,
                            const std::filesystem::path& out,
                            const std::string& groups,
                            const std::vector<std::string>& options) const {
    std::ostringstream cmd;
    cmd << "printf '" << groups << "' | " << shell_quote(cfg_.gmx)
        << " trjconv -s " << shell_quote(tpr.string())
        << " -f " << shell_quote(in.string())
        << " -o " << shell_quote(out.string());
    for (const auto& opt : options) cmd << ' ' << opt;
    sh_.runShell(cmd.str());
}

} // namespace cg
