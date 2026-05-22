#include "Workflow.hpp"
#include "FileUtils.hpp"
#include "IndexBuilder.hpp"
#include "MdpWriter.hpp"
#include "PlumedWriter.hpp"
#include "StringUtils.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace cg {
namespace {

void require_file(const std::filesystem::path& path, const std::string& label, bool dry_run) {
    if (!dry_run && !std::filesystem::exists(path)) throw std::runtime_error("Missing " + label + ": " + path.string());
}

void copy_if_distinct(const std::filesystem::path& src, const std::filesystem::path& dst, bool dry_run) {
    if (dry_run) return;
    if (std::filesystem::exists(dst) && std::filesystem::equivalent(src, dst)) return;
    copy_overwrite(src, dst);
}

std::size_t csv_size(const std::string& csv) {
    std::size_t n = 0;
    std::stringstream ss(csv);
    for (std::string item; std::getline(ss, item, ',');)
        if (!trim(item).empty()) ++n;
    return n;
}

std::string load_sigma(const Shell& sh, const Config& cfg) {
    if (!cfg.metad_sigma.empty()) return cfg.metad_sigma;

    const auto params = cfg.metadCvParamsPath();
    require_file(params, "CV parameter pickle; pass --metad-sigma if unavailable", sh.dryRun());
    if (sh.dryRun()) return "SIGMA_FROM_CV_PARAMS";

    const std::string code =
        "import pickle,sys; "
        "p=pickle.load(open(sys.argv[1],'rb')); "
        "v=p['sigma']; "
        "v=v.tolist() if hasattr(v,'tolist') else list(v); "
        "print(','.join(format(float(x),'.12g') for x in v))";

    auto out = trim(sh.runCapture({cfg.python, "-c", code, params.string()}).output);
    if (out.empty()) throw std::runtime_error("Cannot extract sigma from " + params.string());
    return out;
}

std::vector<std::string> walker_dirs(int n) {
    std::vector<std::string> dirs;
    dirs.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) dirs.push_back("walker" + std::to_string(i));
    return dirs;
}

} // namespace

void Workflow::prepare_metadynamics() {
    const auto npt_gro = cfg_.resultDir() / "npt.gro";
    const auto npt_cpt = cfg_.resultDir() / "npt.cpt";
    const auto model_src = cfg_.metadModelPath();
    const auto metad_dir = cfg_.metadDir();
    const auto index_src = cfg_.resultDir() / "index.ndx";
    const auto index_dst = metad_dir / "index.ndx";
    const auto model_dst = metad_dir / "CVs_torchscript.pt";
    const auto plumed_dst = metad_dir / "plumed.dat";

    require_file(npt_gro, "NPT coordinates", sh_.dryRun());
    require_file(npt_cpt, "NPT checkpoint", sh_.dryRun());
    require_file(model_src, "TorchScript CV model", sh_.dryRun());

    cfg_.metad_sigma = load_sigma(sh_, cfg_);
    if (!sh_.dryRun() && csv_size(cfg_.metad_sigma) != static_cast<std::size_t>(cfg_.metad_nodes))
        throw std::runtime_error("--metad-sigma must contain exactly --metad-nodes values");

    std::filesystem::create_directories(metad_dir);
    generate_index(sh_, cfg_, npt_gro);
    copy_if_distinct(index_src, index_dst, sh_.dryRun());
    copy_if_distinct(model_src, model_dst, sh_.dryRun());

    write_metad_mdp(cfg_);
    grompp(cfg_.systemDir() / "md_metad.mdp", npt_gro, metad_dir / "md.tpr", npt_cpt, {}, index_src);
    write_metad_plumed_dat(cfg_, plumed_dst, "index.ndx", "CVs_torchscript.pt", "COLVAR");

    if (cfg_.metad_walkers == 1) return;

    for (const auto& dir_name : walker_dirs(cfg_.metad_walkers)) {
        const auto wdir = metad_dir / dir_name;
        std::filesystem::create_directories(wdir);
        copy_if_distinct(metad_dir / "md.tpr", wdir / "md.tpr", sh_.dryRun());
        copy_if_distinct(index_dst, wdir / "index.ndx", sh_.dryRun());
        copy_if_distinct(model_dst, wdir / "CVs_torchscript.pt", sh_.dryRun());
        copy_if_distinct(plumed_dst, wdir / "plumed.dat", sh_.dryRun());
    }
}

void Workflow::run_metadynamics() {
    std::vector<std::string> cmd;
    const int ntomp = cfg_.ntomp;

    if (cfg_.metad_walkers == 1) {
        cmd = {cfg_.gmx, "mdrun", "-plumed", "plumed.dat", "-v", "-deffnm", "md"};
    } else {
        cmd = {cfg_.mpirun, "-np", std::to_string(cfg_.metad_walkers), cfg_.gmx, "mdrun",
               "-plumed", "plumed.dat", "-multidir"};
        for (const auto& dir : walker_dirs(cfg_.metad_walkers)) cmd.push_back(dir);
        cmd.insert(cmd.end(), {"-v", "-deffnm", "md"});
    }

    if (ntomp > 0) cmd.insert(cmd.end() - 2, {"-ntomp", std::to_string(ntomp)});
    sh_.run(cmd, cfg_.metadDir());
}

void Workflow::finalize_metadynamics() {
    const auto metad_dir = cfg_.metadDir();

    if (cfg_.metad_sum_hills) {
        const auto hills = cfg_.metad_walkers > 1 ? metad_dir / "walker0" / "HILLS" : metad_dir / "HILLS";
        if (std::filesystem::exists(hills) || sh_.dryRun()) {
            sh_.run({cfg_.plumed, "sum_hills", "--hills", hills.string(),
                     "--outfile", (metad_dir / "fes.dat").string(), "--mintozero"});
        } else {
            std::cerr << "Skipping FES: missing HILLS file " << hills << '\n';
        }
    }

    if (!cfg_.metad_center) return;

    const auto centers = cfg_.metad_walkers == 1 ? std::vector<std::filesystem::path>{metad_dir}
                                                : [&] {
                                                      std::vector<std::filesystem::path> dirs;
                                                      for (const auto& d : walker_dirs(cfg_.metad_walkers)) dirs.push_back(metad_dir / d);
                                                      return dirs;
                                                  }();

    for (const auto& dir : centers) {
        const auto tpr = dir / "md.tpr";
        const auto xtc = dir / "md.xtc";
        const auto gro = dir / "md.gro";
        if (std::filesystem::exists(xtc) || sh_.dryRun())
            trjconv(tpr, xtc, dir / "md_center.xtc", "1\n0\n", {"-center", "-pbc mol", "-ur compact"});
        if (std::filesystem::exists(gro) || sh_.dryRun()) {
            trjconv(tpr, gro, dir / "md_center.gro", "1\n1\n", {"-center", "-pbc mol"});
            sh_.run({cfg_.gmx, "editconf", "-f", (dir / "md.gro").string(), "-o", (dir / "md.pdb").string()});
        }
    }
}

void Workflow::run_metadynamics_workflow() {
    prepare_metadynamics();
    run_metadynamics();
    finalize_metadynamics();
}

} // namespace cg
