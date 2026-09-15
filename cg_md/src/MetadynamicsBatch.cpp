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

std::vector<std::string> MetadynamicsRunner::plumed_outputs() const {
    return {kHillsFile, kColvarFile, cfg_.metad_monitor_file};
}

void MetadynamicsRunner::record_plumed_lengths(const std::filesystem::path& dir, int batch) const {
    const auto path = dir / "plumed_records.ledger";
    auto ledger = RunLedger::load(path);
    for (const auto& name : plumed_outputs())
        ledger.set("after_" + zero_padded(batch) + "_" + name,
                   std::to_string(count_plumed_records(dir / name)));
    ledger.save(path);
}

void MetadynamicsRunner::rewind_plumed_outputs(const std::filesystem::path& dir, int batch) const {
    const auto ledger = RunLedger::load(dir / "plumed_records.ledger");
    for (const auto& name : plumed_outputs()) {
        const auto key = "after_" + zero_padded(batch - 1) + "_" + name;
        const auto recorded = ledger.get(key);

        std::size_t keep = 0;
        bool have_count = false;
        if (!recorded.empty()) {
            try {
                keep = static_cast<std::size_t>(std::stoull(recorded));
                have_count = true;
            } catch (const std::exception&) {
                std::cerr << "Warning: unreadable record count '" << recorded << "' for " << key
                          << "; falling back to counting the file's clock restarts.\n";
            }
        }

        const std::size_t removed =
            have_count ? truncate_plumed_file_to(dir / name, keep)
                       : truncate_plumed_file_to_segment(dir / name, static_cast<std::size_t>(batch));

        if (removed)
            std::cerr << "Discarded " << removed << " record(s) from " << (dir / name)
                      << " left by an interrupted attempt at batch " << batch << ".\n";
    }
}

void MetadynamicsRunner::prepare_batch(int batch) {
    const auto name = batch_name(batch);
    const auto batch_steps = cfg_.metadChunkNsteps();
    const auto done_steps = static_cast<std::int64_t>(batch) * batch_steps;
    const double tinit_ps = static_cast<double>(done_steps) * cfg_.md_dt_ps;

    MetadPlumedOptions options;
    options.restart = batch > 0;
    options.hills_file = kHillsFile;

    const auto previous = batch_name(batch - 1);
    int walker = 0;
    for (const auto& dir : run_dirs()) {
        const auto mdp = write_metad_batch_mdp(cfg_, batch, batch_steps, tinit_ps, walker++);
        if (batch > 0 && !sh_.dryRun()) rewind_plumed_outputs(dir, batch);

        const auto input_gro = batch == 0 ? cfg_.resultDir() / "npt.gro" : dir / (previous + ".gro");
        const auto input_cpt = batch == 0 ? cfg_.resultDir() / "npt.cpt" : dir / (previous + ".cpt");
        require_file(input_gro, "metadynamics batch input coordinates", sh_.dryRun());
        require_file(input_cpt, "metadynamics batch input checkpoint", sh_.dryRun());

        gmx_.grompp(mdp, input_gro, dir / (name + ".tpr"), input_cpt, {}, dir / "index.ndx");
        write_metad_plumed_dat(cfg_, dir / "plumed.dat", "index.ndx", "CVs_torchscript.pt",
                               kColvarFile, options);
    }
}

void MetadynamicsRunner::run_batch(int batch) {
    const auto name = batch_name(batch);
    std::vector<std::string> cmd;

    if (cfg_.metad_walkers == 1) {
        cmd = {cfg_.gmx, "mdrun", "-plumed", "plumed.dat", "-v", "-deffnm", name};
    } else {
        cmd = {cfg_.mpirun, "-np", std::to_string(cfg_.metad_walkers),
               cfg_.gmx, "mdrun", "-plumed", "plumed.dat", "-multidir"};
        for (const auto& dir : walker_dirs(cfg_.metad_walkers)) cmd.push_back(dir);
        cmd.insert(cmd.end(), {"-v", "-deffnm", name});
    }

    if (cfg_.ntomp > 0) cmd.insert(cmd.end(), {"-ntomp", std::to_string(cfg_.ntomp)});
    append_pinning(cfg_, cmd);

    sh_.run(cmd, cfg_.metadDir());

    if (!sh_.dryRun())
        for (const auto& dir : run_dirs()) record_plumed_lengths(dir, batch);
}

} // namespace cg
