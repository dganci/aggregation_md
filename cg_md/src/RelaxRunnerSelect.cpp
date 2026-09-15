#include "RelaxRunner.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "GromacsDriver.hpp"
#include "GroUtils.hpp"
#include "MdpWriter.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"
#include "TopologyEditor.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace cg {

void RelaxRunner::write_distinct_copies() {
    const double total_ps = cfg_.relax_us * 1'000'000.0;
    const double start_ps = cfg_.relax_discard_frac * total_ps;
    const double span_ps = total_ps - start_ps;
    if (cfg_.n_prot < 1 || span_ps <= 0.0) return;

    const double spacing_ps = cfg_.n_prot > 1 ? span_ps / static_cast<double>(cfg_.n_prot - 1)
                                              : span_ps;
    const auto multi = dir() / "copies.pdb";
    gmx_.trjconv(dir() / "md.tpr", dir() / "md.xtc", multi, "1\n",
                 {"-b", to_string_fixed(start_ps, 1), "-dt", to_string_fixed(spacing_ps, 1),
                  "-pbc", "whole"});

    if (sh_.dryRun()) return;

    std::vector<std::vector<std::string>> models;
    std::vector<std::string> cur;
    for (const auto& line : read_lines(multi)) {
        if (starts_with(line, "ENDMDL")) {
            if (!cur.empty()) { cur.push_back("END"); models.push_back(cur); cur.clear(); }
            continue;
        }
        if (starts_with(line, "MODEL") || starts_with(line, "TITLE") ||
            starts_with(line, "REMARK") || starts_with(line, "CRYST1") ||
            starts_with(line, "END")) continue;
        cur.push_back(line);
    }
    if (!cur.empty()) { cur.push_back("END"); models.push_back(cur); }

    if (static_cast<int>(models.size()) < cfg_.n_prot)
        throw std::runtime_error(
            "--distinct-copies needs " + std::to_string(cfg_.n_prot) + " frames from the "
            "relaxation, but trjconv produced " + std::to_string(models.size()) + " from " +
            (dir() / "md.xtc").string() + ".\n"
            "The relaxation trajectory is too short or too sparsely written for the number of "
            "copies requested. Raise --relax-us, or raise --relax-frames (which is what sets the "
            "relaxation's own write-out stride - --nst-xtc does not apply here), or drop "
            "--distinct-copies to replicate one representative instead.");

    for (int i = 0; i < cfg_.n_prot; ++i) write_lines(cfg_.cgRelaxedCopyPath(i), models[i]);

    std::cerr << "Relaxation done: " << cfg_.n_prot << " distinct starting conformations written to "
              << cfg_.cgRelaxedCopyPath(0).parent_path() << ",\n"
              << "one per copy, spaced " << to_string_fixed(spacing_ps / 1000.0, 1)
              << " ns apart over the equilibrated part of the trajectory.\n";
}

void RelaxRunner::select_representative() {
    const auto discard_ps = cfg_.relax_discard_frac * cfg_.relax_us * 1'000'000.0;

    sh_.runShell("printf '1\\n1\\n' | " + shell_quote(cfg_.gmx) + " cluster" +
                 " -s " + shell_quote((dir() / "md.tpr").string()) +
                 " -f " + shell_quote((dir() / "md_prot.xtc").string()) +
                 " -method gromos" +
                 " -o " + shell_quote((dir() / "rmsd-clust.xpm").string()) +
                 " -cutoff " + to_string_fixed(cfg_.relax_cluster_cutoff_nm, 3) +
                 " -b " + to_string_fixed(discard_ps, 1) +
                 " -cl " + shell_quote((dir() / "clusters.pdb").string()) +
                 " -g " + shell_quote((dir() / "cluster.log").string()) +
                 " -dist " + shell_quote((dir() / "rmsd_dist.xvg").string()) +
                 " -sz " + shell_quote((dir() / "cluster_sizes.xvg").string()));

    if (sh_.dryRun()) return;

    const auto clusters = dir() / "clusters.pdb";
    if (!std::filesystem::exists(clusters))
        throw std::runtime_error("Clustering produced no structures: " + clusters.string());

    std::vector<std::string> first;
    for (const auto& line : read_lines(clusters)) {
        if (starts_with(line, "ENDMDL")) break;
        if (starts_with(line, "MODEL") || starts_with(line, "TITLE") || starts_with(line, "REMARK")) continue;
        first.push_back(line);
    }
    if (first.empty()) throw std::runtime_error("Could not extract a representative from " + clusters.string());
    first.push_back("END");

    write_lines(cfg_.cgRelaxedPath(), first);

    std::cerr << "Relaxation done: " << cfg_.cgRelaxedPath() << " holds the centroid of the most\n"
              << "populated cluster, and Packmol will now replicate it instead of the raw structure.\n"
              << "Cluster populations: " << (dir() / "cluster_sizes.xvg") << '\n';
}

} // namespace cg
