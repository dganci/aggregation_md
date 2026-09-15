#include "AdaptiveSampler.hpp"
#include "Config.hpp"
#include "GromacsDriver.hpp"
#include "Reporter.hpp"
#include "Shell.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
namespace cg {

void AdaptiveSampler::center_chunks() const {
    try {
        const auto chunks = completed_chunks();
        if (chunks.empty()) return;

        std::vector<std::string> xtcs;
        for (const auto chunk : chunks) {
            const auto path = chunk_base(chunk).string() + ".xtc";
            if (std::filesystem::exists(path)) xtcs.push_back(path);
        }
        if (xtcs.empty()) return;

        const auto joined = cfg_.resultDir() / "md_all.xtc";
        const auto tpr = chunk_tpr(chunks.front());

        std::vector<std::string> cmd = {cfg_.gmx, "trjcat", "-f"};
        cmd.insert(cmd.end(), xtcs.begin(), xtcs.end());
        cmd.insert(cmd.end(), {"-o", joined.string(), "-cat"});
        sh_.run(cmd);

        const auto whole = cfg_.resultDir() / "md_all_whole.xtc";
        gmx_.trjconv(tpr, joined, whole, "0\n", {"-pbc", "whole"});
        gmx_.trjconv(tpr, whole, cfg_.resultDir() / "md_all_center.xtc",
                     "1\n0\n", {"-center", "-pbc", "mol", "-ur", "compact"});
        gmx_.trjconv(tpr, cfg_.resultDir() / "md_all_center.xtc",
                     cfg_.resultDir() / "md_all_center.gro", "0\n", {"-dump", "0"});

        if (!cfg_.keep_intermediate_trajectories) {
            std::error_code ec;
            for (const char* name : {"md_all.xtc", "md_all_whole.xtc"})
                std::filesystem::remove(cfg_.resultDir() / name, ec);
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: could not concatenate/centre the chunk trajectories (ignored): "
                  << e.what() << '\n';
    }
}

void AdaptiveSampler::report() const {
    center_chunks();
    run_report(cfg_, sh_, cfg_.resultDir(), "adaptive");
}

} // namespace cg
