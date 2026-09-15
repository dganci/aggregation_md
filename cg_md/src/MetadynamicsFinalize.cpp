#include "MetadynamicsRunner.hpp"
#include "MetadynamicsInternal.hpp"
#include "ColvarTable.hpp"
#include "Config.hpp"
#include "GromacsDriver.hpp"
#include "PlumedLog.hpp"
#include "FileUtils.hpp"
#include "Reporter.hpp"
#include "RunLedger.hpp"
#include "RunDiagnostics.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
namespace cg {
namespace {

std::string last_batch_basename(const std::filesystem::path& dir) {
    std::string best;
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(dir, ec)) {
        if (ec) break;
        if (entry.path().extension() != ".tpr") continue;
        const auto stem = entry.path().stem().string();
        if (!starts_with(stem, "md_")) continue;
        const auto suffix = stem.substr(3);
        if (suffix.empty() || suffix.find_first_not_of("0123456789") != std::string::npos) continue;
        if (best.empty() || stem > best) best = stem;
    }
    return best.empty() ? std::string("md") : best;
}
} // namespace

void MetadynamicsRunner::analyze() {
    if (sh_.dryRun()) return;

    const auto dirs = existing_run_dirs();
    const double dt_colvar_ps = cfg_.md_dt_ps * static_cast<double>(cfg_.plumed_stride);
    const int batches = std::max(1, cfg_.metadNumChunks());

    for (std::size_t i = 0; i < dirs.size(); ++i) {
        const auto& dir = dirs[i];
        const auto label = dirs.size() == 1 ? std::string("metad") : "metad_walker" + std::to_string(i);
        try {
            RunDiagnostics diag;

            std::filesystem::path last_edr;
            for (int b = 0; b < batches; ++b) {
                const auto edr = dir / (batch_name(b) + ".edr");
                if (std::filesystem::exists(edr)) last_edr = edr;
            }
            if (!last_edr.empty())
                diag.energy = compute_energy_diagnostics(gmx_, last_edr, dir / "diagnostics_scratch",
                                                          /*include_pressure=*/true, /*include_density=*/true);

            const auto monitor = dir / cfg_.metad_monitor_file;
            if (std::filesystem::exists(monitor)) {
                diag.structure = compute_structural_diagnostics({monitor}, cfg_.n_prot, dt_colvar_ps,
                                                                 cfg_.adaptive_pair_contact_threshold);
                diag.has_structure = true;
            }
            append_diagnostics_jsonl(cfg_.metadDir() / "diagnostics.jsonl", diag, label);
        } catch (const std::exception& e) {
            std::cerr << "Warning: " << label << " diagnostics failed (ignored, purely observational): "
                      << e.what() << '\n';
        }
    }
}

void MetadynamicsRunner::finalize() {
    const auto metad_dir = cfg_.metadDir();
    const auto dirs = sh_.dryRun() ? run_dirs() : existing_run_dirs();

    if (cfg_.metad_sum_hills) {
        std::vector<std::string> hills;
        for (const auto& dir : dirs) {
            const auto path = dir / kHillsFile;
            if (std::filesystem::exists(path) || sh_.dryRun()) hills.push_back(path.string());
        }

        if (hills.empty()) {
            std::cerr << "Skipping FES: no HILLS file found under " << metad_dir << '\n';
        } else {
            sh_.run({cfg_.plumed, "sum_hills", "--hills", join(hills, ","),
                     "--outfile", (metad_dir / "fes.dat").string(), "--mintozero"},
                    metad_dir);

            if (cfg_.metad_fes_stride > 0) {
                const auto conv_dir = metad_dir / "fes_convergence";
                if (!sh_.dryRun()) std::filesystem::create_directories(conv_dir);
                sh_.run({cfg_.plumed, "sum_hills", "--hills", join(hills, ","),
                         "--outfile", (conv_dir / "fes_").string(),
                         "--stride", std::to_string(cfg_.metad_fes_stride), "--mintozero"},
                        metad_dir);
            }
        }
    }

    if (cfg_.metad_center) {
        for (const auto& dir : dirs) {
            const auto name = last_batch_basename(dir);
            const auto tpr = dir / (name + ".tpr");
            const auto xtc = dir / (name + ".xtc");
            const auto gro = dir / (name + ".gro");
            if (std::filesystem::exists(xtc) || sh_.dryRun())
                gmx_.trjconv(tpr, xtc, dir / "md_center.xtc", "1\n0\n", {"-center", "-pbc mol", "-ur compact"});
            if (std::filesystem::exists(gro) || sh_.dryRun()) {
                gmx_.trjconv(tpr, gro, dir / "md_center.gro", "1\n0\n", {"-center", "-pbc mol"});
                sh_.run({cfg_.gmx, "editconf", "-f", (dir / "md_center.gro").string(),
                         "-o", (dir / "md_center.pdb").string()});
            }
        }
    }

    run_report(cfg_, sh_, metad_dir, "metadynamics");
}

} // namespace cg
