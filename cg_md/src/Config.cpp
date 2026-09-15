#include "Config.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/stat.h>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>


namespace cg {
namespace {

constexpr double ps_per_us = 1'000'000.0;

template <class T>
std::int64_t steps_from_us(T us, double dt_ps) {
    return static_cast<std::int64_t>(static_cast<double>(us) * ps_per_us / dt_ps + 0.5);
}


std::string executable_build_time() {
    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (ec) return "unknown";
    struct stat st {};
    if (::stat(exe.c_str(), &st) != 0) return "unknown";
    char buf[32];
    if (std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", std::localtime(&st.st_mtime)) == 0)
        return "unknown";
    return buf;
}

} // namespace

std::string Config::systemName() const { return std::to_string(n_prot) + "x" + protomer_name; }

std::filesystem::path Config::cgRelaxedCopyPath(int i) const {
    return prepDir() / (protomer_name + "_cg_relaxed_" + std::to_string(i) + ".pdb");
}

std::vector<std::filesystem::path> Config::packmolSourcePaths() const {
    if (distinct_copies) {
        std::vector<std::filesystem::path> out;
        for (int i = 0; i < n_prot; ++i) {
            const auto p = cgRelaxedCopyPath(i);
            if (!std::filesystem::exists(p)) { out.clear(); break; }
            out.push_back(p);
        }
        if (!out.empty()) return out;
    }
    return {packmolSourcePath()};
}

std::filesystem::path Config::packmolSourcePath() const {
    return std::filesystem::exists(cgRelaxedPath()) ? cgRelaxedPath() : cgPath();
}
std::filesystem::path Config::topologyPath() const { return prepDir() / (protomer_name + "_topo.top"); }
std::filesystem::path Config::packRawPath() const { return prepDir() / (systemName() + "_packed_raw.pdb"); }
std::filesystem::path Config::packCleanPath() const { return prepDir() / (systemName() + "_packed.pdb"); }
std::filesystem::path Config::packInputPath() const { return prepDir() / (systemName() + ".inp"); }
std::filesystem::path Config::boxedPath() const { return prepDir() / (systemName() + "_boxed.gro"); }
std::filesystem::path Config::solvatedPath() const { return prepDir() / (systemName() + "_solvated.gro"); }
std::filesystem::path Config::ionsTprPath() const { return resultDir() / "ions.tpr"; }
std::string Config::martiniMainItp() const {
    if (std::abs(martini_lambda_pw - 1.0) < 1e-9) return "martini_v3.0.0.itp";
    return "martini_v3.0.0_lpw" + std::to_string(static_cast<int>(martini_lambda_pw * 100 + 0.5)) + ".itp";
}

std::filesystem::path Config::martiniWaterPath() const { return martini_water_gro.is_absolute() ? martini_water_gro : project_dir / martini_water_gro; }
std::filesystem::path Config::martinizeFfDirPath() const { return martinize_ff_dir.is_absolute() ? martinize_ff_dir : project_dir / martinize_ff_dir; }
std::filesystem::path Config::martinizeMapDirPath() const { return martinize_map_dir.is_absolute() ? martinize_map_dir : project_dir / martinize_map_dir; }
std::filesystem::path Config::beadsPath() const { return prepDir() / (systemName() + "_beads.gro"); }
std::filesystem::path Config::systemDir() const { return runDir() / "system"; }
std::filesystem::path Config::resultDir() const { return runDir(); }
std::filesystem::path Config::metadDir() const { return resultDir() / (pmf ? "pmf" : "metadynamics"); }
std::filesystem::path Config::metadModelPath() const {
    if (!metad_model.empty()) return metad_model.is_absolute() ? metad_model : project_dir / metad_model;
    return metadDir() / "CVs" / protomer_name / "CVs_torchscript.pt";
}

std::filesystem::path Config::metadCvParamsPath() const {
    if (!metad_cv_params.empty()) return metad_cv_params.is_absolute() ? metad_cv_params : project_dir / metad_cv_params;
    return metadDir() / "CVs" / protomer_name / "cv_params.pkl";
}

std::int64_t Config::mdNsteps() const { return steps_from_us(md_total_us, md_dt_ps); }
std::int64_t Config::adaptiveChunkNsteps() const { return steps_from_us(adaptive_chunk_us, md_dt_ps); }
int Config::adaptiveMaxChunks() const {
    if (adaptive_chunk_us <= 0.0) return 0;
    return static_cast<int>(std::ceil(adaptive_max_total_us / adaptive_chunk_us - 1e-9));
}

double Config::metadTotalUs() const { return metad_total_us > 0.0 ? metad_total_us : md_total_us; }
std::int64_t Config::metadNsteps() const { return steps_from_us(metadTotalUs(), md_dt_ps); }

std::int64_t Config::metadChunkNsteps() const {
    if (metad_chunk_us <= 0.0) return metadNsteps();
    return std::min(steps_from_us(metad_chunk_us, md_dt_ps), metadNsteps());
}

int Config::metadNumChunks() const {
    const auto total = metadNsteps();
    const auto chunk = metadChunkNsteps();
    if (total <= 0 || chunk <= 0) return 0;
    return static_cast<int>((total + chunk - 1) / chunk);
}

std::int64_t Config::effectiveMetadNsteps() const {
    return static_cast<std::int64_t>(metadNumChunks()) * metadChunkNsteps();
}


void Config::print() const {
    std::cerr
        << "Configuration\n"
        << "  cg_md built       = " << executable_build_time()
        << "   (rebuild if this predates your last source edit)\n"
        << "  project_dir       = " << project_dir << '\n'
        << "  system            = " << systemName() << '\n'
        << "  run               = " << runName()
        << (run_tag.empty() ? "  (id derived from the parameters)" : "  (--run-tag)") << '\n'
        << "  n_prot            = " << n_prot << '\n'
        << "  protomer_name     = " << protomer_name << '\n'
        << "  seed              = " << seed << (seed < 0 ? " (random, non-reproducible)" : "") << '\n'
        << "  seq_length        = " << seq_length << '\n'
        << "  atoms_per_prot    = " << atoms_per_prot << '\n'
        << "  is_phospho        = " << (is_phospho ? "true" : "false") << '\n'
        << "  md_dt_ps          = " << md_dt_ps << '\n'
        << "  md_total_us       = " << md_total_us << '\n'
        << "  md_nsteps         = " << mdNsteps() << '\n'
        << "  plumed_stride     = " << plumed_stride << '\n'
        << "  thermostat_mode   = " << thermostat_mode << '\n'
        << "  epsilon_r         = " << epsilon_r << (epsilon_r == 15.0 ? " (Martini standard water)" : " (NON-STANDARD for Martini)") << '\n'
        << "  rcoulomb/rvdw_nm  = " << rcoulomb_nm << " / " << rvdw_nm << '\n'
        << "  solvation_mode    = " << solvation_mode << '\n'
        << "  packmol_R_A       = " << packmol_cluster_radius_A << '\n'
        << "  box_type          = " << box_type << '\n'
        << "  box_margin_nm     = " << box_margin_nm << '\n'
        << "  solvate_radius_nm = " << solvate_radius_nm << '\n'
        << "  martini_water_gro = " << martiniWaterPath() << '\n'
        << "  add_ions          = " << (add_ions ? "true" : "false") << '\n'
        << "  pdb_path          = " << pdbPath() << '\n'
        << "  topology_path     = " << topologyPath() << '\n'
        << "  result_dir        = " << resultDir() << '\n'
        << "  modify            = " << (modify.empty() ? "<off>" : modify) << '\n'
        << "  adaptive          = " << (adaptive ? "true" : "false") << '\n'
        << "  adaptive_chunk_us = " << adaptive_chunk_us << '\n'
        << "  adaptive_max_us   = " << adaptive_max_total_us << '\n'
        << "  metadynamics      = " << (metadynamics ? "true" : "false") << '\n';
    if (is_phospho) {
        std::cerr
            << "  martinize_ff_dir  = " << martinizeFfDirPath() << '\n'
            << "  martinize_map_dir = " << martinizeMapDirPath() << '\n';
    }
    if (metadynamics) {
        std::cerr
            << "  metad_dir         = " << metadDir() << '\n'
            << "  metad_model       = " << metadModelPath() << '\n'
            << "  metad_cv_params   = " << metadCvParamsPath() << '\n'
            << "  metad_walkers     = " << metad_walkers << '\n'
            << "  metad_total_us    = " << metadTotalUs()
            << (effectiveMetadNsteps() != metadNsteps()
                    ? "  -> " + std::to_string(effectiveMetadNsteps() * md_dt_ps / 1e6) +
                          " us (rounded up to a whole number of batches)"
                    : "")
            << '\n'
            << "  metad_nsteps      = " << effectiveMetadNsteps() << '\n'
            << "  metad_chunk_us    = " << (metad_chunk_us > 0.0 ? std::to_string(metad_chunk_us) : "<single batch>") << '\n'
            << "  metad_batches     = " << metadNumChunks() << '\n'
            << "  metad_grid        = " << (metad_grid ? "true" : "false") << '\n'
            << "  metad_calc_rct    = " << (metad_calc_rct ? "true" : "false") << '\n'
            << "  metad_sigma       = " << (metad_sigma.empty() ? "<from cv_params.pkl>" : metad_sigma) << '\n';
    }
}

} // namespace cg
