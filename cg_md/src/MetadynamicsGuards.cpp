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

void MetadynamicsRunner::check_walker_sharing() const {
    if (sh_.dryRun() || cfg_.metad_walkers <= 1) return;

    const auto log = run_dirs().front() / "PLUMED.OUT";
    if (!std::filesystem::exists(log)) {
        std::cerr << "Note: no " << log << ", so the walker count could not be checked.\n";
        return;
    }

    const int joined = plumed_walker_count(read_tail(log, 256 * 1024));
    if (joined == cfg_.metad_walkers) {
        std::cerr << "PLUMED reports " << joined << " walkers sharing hills.\n";
        return;
    }
    if (joined < 0) {
        std::cerr << "Note: " << log << " does not state how many walkers are active; "
                  << "check it says " << cfg_.metad_walkers << ", not 1.\n";
        return;
    }

    throw std::runtime_error(
        "PLUMED joined " + std::to_string(joined) + " walker(s) but --metad-walkers is " +
        std::to_string(cfg_.metad_walkers) + " (" + log.string() + "). The walkers are not "
        "sharing hills, so this is " + std::to_string(cfg_.metad_walkers) + " independent "
        "metadynamics runs and the summed FES would be wrong without anything failing.\n"
        "Usually GROMACS was built without `plumed patch`: the native PLUMED interface does not "
        "pass the multi-simulation communicator. See cgmd_build_lib/05_gromacs.sh.\n"
        "This batch ran without sharing, so delete " + cfg_.metadDir().string() +
        " before resuming - a completed batch is not re-run.");
}

std::vector<std::size_t> MetadynamicsRunner::plumed_backup_counts() const {
    std::vector<std::size_t> counts;
    for (const auto& dir : run_dirs()) {
        std::size_t n = 0;
        std::error_code ec;
        for (const auto& entry : std::filesystem::directory_iterator(dir, ec)) {
            if (ec) break;
            const auto name = entry.path().filename().string();
            for (const auto& output : plumed_outputs())
                if (is_plumed_backup(name, output)) { ++n; break; }
        }
        counts.push_back(n);
    }
    return counts;
}

void MetadynamicsRunner::check_bias_continuity(int batch,
                                                const std::vector<std::size_t>& before) const {
    if (sh_.dryRun() || batch <= 0) return;

    const auto dirs = run_dirs();
    const auto after = plumed_backup_counts();
    if (after.size() != dirs.size() || before.size() != dirs.size()) return;

    for (std::size_t i = 0; i < dirs.size(); ++i) {
        const auto ledger = RunLedger::load(dirs[i] / "plumed_records.ledger");
        for (const auto& name : plumed_outputs()) {
            const auto previous = ledger.get("after_" + zero_padded(batch - 1) + "_" + name);
            const auto current = ledger.get("after_" + zero_padded(batch) + "_" + name);
            if (previous.empty() || current.empty()) continue;
            std::size_t was = 0, now = 0;
            try {
                was = static_cast<std::size_t>(std::stoull(previous));
                now = static_cast<std::size_t>(std::stoull(current));
            } catch (const std::exception&) {
                continue;
            }
            if (now >= was) continue;

            throw std::runtime_error(
                (dirs[i] / name).string() + " shrank across batch " + std::to_string(batch) +
                ": " + std::to_string(was) + " records after batch " + std::to_string(batch - 1) +
                ", " + std::to_string(now) + " now. A PLUMED output only shrinks if it was "
                "reopened instead of appended to, i.e. RESTART did not take effect, so the bias "
                "from batches 0.." + std::to_string(batch - 1) + " was thrown away.\n"
                "No bck.* was left, so backups are off (PLUMED_MAXBACKUP=0) and those records are "
                "gone. Fix the RESTART line in " + (dirs[i] / "plumed.dat").string() +
                " and restart the run.");
        }

        if (after[i] <= before[i]) continue;

        throw std::runtime_error(
            "Batch " + std::to_string(batch) + " backed up its PLUMED output in " +
            dirs[i].string() + " instead of appending. PLUMED does that when RESTART is not in "
            "effect: it renamed HILLS to bck.0.HILLS and resumed from a flat bias, discarding the "
            "hills from batches 0.." + std::to_string(batch - 1) + ".\n"
            "Those hills are still in the bck.* file - do not delete it. Move it back over HILLS "
            "and re-run this batch; the PLUMED outputs are rewound automatically.\n"
            "cg_md writes RESTART into every batch after the first, so check the top of " +
            (dirs[i] / "plumed.dat").string() + ".");
    }
}

} // namespace cg
