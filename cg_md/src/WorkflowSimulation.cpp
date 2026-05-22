#include "Workflow.hpp"
#include "IndexBuilder.hpp"
#include "PlumedWriter.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <filesystem>
#include <sstream>

namespace cg {

void Workflow::grompp(const std::filesystem::path& mdp,
                      const std::filesystem::path& coordinates,
                      const std::filesystem::path& output,
                      const std::filesystem::path& checkpoint,
                      const std::filesystem::path& restraints,
                      const std::filesystem::path& index) {
    std::vector<std::string> cmd = {cfg_.gmx, "grompp", "-f", mdp.string(), "-c", coordinates.string(),
                                    "-p", cfg_.topologyPath().string(), "-o", output.string()};
    if (!checkpoint.empty()) cmd.insert(cmd.end(), {"-t", checkpoint.string()});
    if (!restraints.empty()) cmd.insert(cmd.end(), {"-r", restraints.string()});
    if (!index.empty()) cmd.insert(cmd.end(), {"-n", index.string()});
    sh_.run(cmd);
}

void Workflow::mdrun(const std::filesystem::path& deffnm,
                     const std::filesystem::path& tpr,
                     const std::filesystem::path& plumed,
                     int ntomp_override) {
    std::vector<std::string> cmd = {cfg_.gmx, "mdrun", "-v", "-deffnm", deffnm.string()};
    if (!tpr.empty()) cmd.insert(cmd.end(), {"-s", tpr.string()});
    if (!plumed.empty()) cmd.insert(cmd.end(), {"-plumed", plumed.string()});

    const int ntomp = ntomp_override > 0 ? ntomp_override : cfg_.ntomp;
    if (ntomp > 0) cmd.insert(cmd.end(), {"-ntomp", std::to_string(ntomp)});
    sh_.run(cmd);
}

void Workflow::energy(const std::filesystem::path& edr,
                      const std::filesystem::path& xvg,
                      const std::string& observable) {
    sh_.runShell("printf '" + observable + "\\n0\\n' | " + shell_quote(cfg_.gmx) +
                 " energy -f " + shell_quote(edr.string()) +
                 " -o " + shell_quote(xvg.string()));
}

void Workflow::trjconv(const std::filesystem::path& tpr,
                       const std::filesystem::path& in,
                       const std::filesystem::path& out,
                       const std::string& groups,
                       const std::vector<std::string>& options) {
    std::ostringstream cmd;
    cmd << "printf '" << groups << "' | " << shell_quote(cfg_.gmx)
        << " trjconv -s " << shell_quote(tpr.string())
        << " -f " << shell_quote(in.string())
        << " -o " << shell_quote(out.string());
    for (const auto& opt : options) cmd << ' ' << opt;
    sh_.runShell(cmd.str());
}

void Workflow::prepare_em() {
    normalize_ions();
    grompp(cfg_.systemDir() / "emin.mdp", cfg_.beadsPath(), cfg_.resultDir() / "em.tpr");
}

void Workflow::run_em() {
    mdrun(cfg_.resultDir() / "em", cfg_.resultDir() / "em.tpr", {}, 1);
}

void Workflow::analyze_em() {
    energy(cfg_.resultDir() / "em.edr", cfg_.resultDir() / "potential.xvg", "Potential");
}

void Workflow::prepare_nvt() {
    generate_index(sh_, cfg_, cfg_.resultDir() / "em.gro");
    grompp(cfg_.systemDir() / "nvt.mdp", cfg_.resultDir() / "em.gro", cfg_.resultDir() / "nvt.tpr", {}, {}, cfg_.resultDir() / "index.ndx");
}

void Workflow::run_nvt() {
    mdrun(cfg_.resultDir() / "nvt");
}

void Workflow::prepare_npt() {
    generate_index(sh_, cfg_, cfg_.resultDir() / "nvt.gro");
    grompp(cfg_.systemDir() / "npt.mdp", cfg_.resultDir() / "nvt.gro", cfg_.resultDir() / "npt.tpr",
           cfg_.resultDir() / "nvt.cpt", cfg_.resultDir() / "nvt.gro", cfg_.resultDir() / "index.ndx");
}

void Workflow::run_npt() {
    mdrun(cfg_.resultDir() / "npt");
}

void Workflow::analyze_npt() {
    energy(cfg_.resultDir() / "npt.edr", cfg_.resultDir() / "pressure.xvg", "Pressure");
    energy(cfg_.resultDir() / "npt.edr", cfg_.resultDir() / "density.xvg", "Density");
}

void Workflow::prepare_production() {
    generate_index(sh_, cfg_, cfg_.resultDir() / "npt.gro");
    write_plumed_dat(cfg_);
    grompp(cfg_.systemDir() / "md.mdp", cfg_.resultDir() / "npt.gro", cfg_.resultDir() / "md.tpr",
           cfg_.resultDir() / "npt.cpt", {}, cfg_.resultDir() / "index.ndx");
}

void Workflow::run_production() {
    mdrun(cfg_.resultDir() / "md", {}, cfg_.resultDir() / "plumed.dat");
}

void Workflow::center_trajectory() {
    const auto tpr = cfg_.resultDir() / "md.tpr";
    const auto xtc = cfg_.resultDir() / "md.xtc";
    const auto nojump = cfg_.resultDir() / "md_nojump.xtc";
    const auto centered_xtc = cfg_.resultDir() / "md_center.xtc";
    const auto centered_gro = cfg_.resultDir() / "md_center.gro";

    trjconv(tpr, xtc, nojump, "0\n", {"-pbc nojump"});

    trjconv(tpr, nojump, centered_xtc,
            "1\n0\n", {"-center", "-pbc mol", "-ur compact"});

    trjconv(tpr, centered_xtc, centered_gro,
            "0\n", {"-dump", "0"});
}

} // namespace cg
