#include "MetadynamicsRunner.hpp"
#include "MetadynamicsInternal.hpp"
#include "ColvarTable.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "GromacsDriver.hpp"
#include "IndexBuilder.hpp"
#include "MdpWriter.hpp"
#include "PlumedWriter.hpp"
#include "Reporter.hpp"
#include "RunDiagnostics.hpp"
#include "RunLedger.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace cg {

MetadynamicsRunner::MetadynamicsRunner(Config& cfg, Shell& sh, GromacsDriver& gmx)
    : cfg_(cfg), sh_(sh), gmx_(gmx) {}

std::vector<std::filesystem::path> MetadynamicsRunner::run_dirs() const {
    if (cfg_.metad_walkers == 1) return {cfg_.metadDir()};
    std::vector<std::filesystem::path> out;
    for (const auto& d : walker_dirs(cfg_.metad_walkers)) out.push_back(cfg_.metadDir() / d);
    return out;
}

std::vector<std::filesystem::path> MetadynamicsRunner::existing_run_dirs() const {
    std::vector<std::filesystem::path> discovered;
    for (int i = 0;; ++i) {
        const auto dir = cfg_.metadDir() / ("walker" + std::to_string(i));
        if (!std::filesystem::is_directory(dir)) break;
        discovered.push_back(dir);
    }
    if (!discovered.empty()) return discovered;
    return {cfg_.metadDir()};
}

std::string MetadynamicsRunner::batch_name(int batch) const {
    return "md_" + zero_padded(batch);
}

bool MetadynamicsRunner::batch_completed(int batch) const {
    if (batch < 0) return true;
    const auto base = batch_name(batch);
    for (const auto& dir : run_dirs()) {
        const auto log = dir / (base + ".log");
        if (!std::filesystem::exists(log)) return false;
        if (read_text(log).find("Finished mdrun") == std::string::npos) return false;
        if (!std::filesystem::exists(dir / (base + ".gro"))) return false;
        if (!std::filesystem::exists(dir / (base + ".cpt"))) return false;
    }
    return true;
}

int MetadynamicsRunner::first_pending_batch() const {
    const int total = cfg_.metadNumChunks();
    for (int i = 0; i < total; ++i)
        if (!batch_completed(i)) return i;
    return total;
}

void MetadynamicsRunner::prepare_common() {
    const auto npt_gro = cfg_.resultDir() / "npt.gro";
    const auto npt_cpt = cfg_.resultDir() / "npt.cpt";
    const auto model_src = cfg_.metadModelPath();
    const auto metad_dir = cfg_.metadDir();

    require_file(npt_gro, "NPT coordinates", sh_.dryRun());
    require_file(npt_cpt, "NPT checkpoint", sh_.dryRun());
    if (!cfg_.pmf) require_file(model_src, "TorchScript CV model", sh_.dryRun());

    load_cv_parameters();

    if (!sh_.dryRun()) std::filesystem::create_directories(metad_dir);

    auto ledger = RunLedger::load(metad_dir / "run.ledger");
    ledger.require("batch_nsteps", std::to_string(cfg_.metadChunkNsteps()), "--metad-chunk-us/--md-dt-ps");
    ledger.require("dt_ps", to_string_fixed(cfg_.md_dt_ps, 6), "--md-dt-ps");
    ledger.require("walkers", std::to_string(cfg_.metad_walkers), "--metad-walkers");
    ledger.require("nodes", std::to_string(cfg_.metad_nodes), "--metad-nodes");
    ledger.require("n_prot", std::to_string(cfg_.n_prot), "--n-prot");
    if (!sh_.dryRun()) {
        ledger.require("sigma", cfg_.metad_sigma, "--metad-sigma / cv_params.pkl");
        ledger.require("grid_min", cfg_.metad_grid_min, "--metad-grid-min / cv_params.pkl");
        ledger.require("grid_max", cfg_.metad_grid_max, "--metad-grid-max / cv_params.pkl");
    }
    if (!sh_.dryRun()) ledger.save(metad_dir / "run.ledger");

    const auto index_src = cfg_.resultDir() / "index.ndx";
    generate_index(sh_, cfg_, npt_gro);

    for (const auto& dir : run_dirs()) {
        if (!sh_.dryRun()) std::filesystem::create_directories(dir);
        copy_if_distinct(index_src, dir / "index.ndx", sh_.dryRun());
        if (!cfg_.pmf) copy_if_distinct(model_src, dir / "CVs_torchscript.pt", sh_.dryRun());
    }
}

void MetadynamicsRunner::run() {
    const int total = cfg_.metadNumChunks();
    if (total <= 0) throw std::runtime_error("Metadynamics has nothing to run (zero steps requested).");

    prepare_common();

    const int start = first_pending_batch();
    if (start >= total) {
        std::cerr << "All " << total << " metadynamics batch(es) already complete; finalizing only.\n";
        analyze();
        finalize();
        return;
    }

    std::cerr << (start ? "Resuming" : "Starting") << " metadynamics from batch " << start << " / " << total
              << " (t = "
              << static_cast<double>(start) * static_cast<double>(cfg_.metadChunkNsteps()) * cfg_.md_dt_ps / 1000.0
              << " ns, bias " << (start ? "restored from HILLS" : "starting flat") << ")\n";

    for (int batch = start; batch < total; ++batch) {
        if (total > 1) {
            std::cerr << "\n============================================================\n"
                      << "METADYNAMICS BATCH " << batch << " / " << total << '\n'
                      << "============================================================\n";
        }
        prepare_batch(batch);
        const auto backups_before = plumed_backup_counts();
        run_batch(batch);
        if (batch == start) check_walker_sharing();
        check_bias_continuity(batch, backups_before);
    }

    analyze();
    finalize();
}

} // namespace cg
