#include "Workflow.hpp"
#include "FileUtils.hpp"
#include "PdbUtils.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
namespace cg {

void Workflow::ensure_dirs() {
    if (sh_.dryRun()) return;

    for (const auto& dir : {cfg_.pdbPath().parent_path(), cfg_.cgPath().parent_path(),
                            cfg_.topologyPath().parent_path(), cfg_.packRawPath().parent_path(),
                            cfg_.packCleanPath().parent_path(), cfg_.packInputPath().parent_path(),
                            cfg_.beadsPath().parent_path(), cfg_.systemDir(), cfg_.resultDir()}) {
        std::filesystem::create_directories(dir);
    }

    write_text(cfg_.resultDir() / "run_params.txt",
               "# Identity of this run. Any change here yields a different run\n"
               "# directory; see Config::identityString().\n"
               "run = " + cfg_.runName() + "\n\n" + cfg_.identityString());
    ensure_residuetypes();
}

void Workflow::ensure_residuetypes() {
    const std::vector<std::string> phospho = {"SEP", "TPO", "PTR"};
    const auto dst = cfg_.resultDir() / "residuetypes.dat";

    std::string text;
    if (std::filesystem::exists(dst)) {
        text = read_text(dst);
    } else {
        std::vector<std::filesystem::path> candidates;
        if (const char* gmxlib = std::getenv("GMXLIB")) candidates.emplace_back(std::filesystem::path(gmxlib) / "residuetypes.dat");
        if (const char* gmxdata = std::getenv("GMXDATA")) candidates.emplace_back(std::filesystem::path(gmxdata) / "top" / "residuetypes.dat");
        for (const char* prefix : {"/opt/gromacs-plumed", "/usr/local/gromacs", "/usr"})
            candidates.emplace_back(std::filesystem::path(prefix) / "share" / "gromacs" / "top" / "residuetypes.dat");
        for (const auto& c : candidates) {
            if (std::filesystem::exists(c)) { text = read_text(c); break; }
        }
        if (text.empty())
            std::cerr << "Note: no stock residuetypes.dat found; writing one with the "
                         "phosphorylated residues only.\n";
    }

    const auto lists = [&text](const std::string& res) {
        std::istringstream in(text);
        std::string line;
        while (std::getline(in, line)) {
            const auto fields = split_ws(line);
            if (!fields.empty() && fields[0] == res) return true;
        }
        return false;
    };

    std::string added;
    for (const auto& res : phospho) {
        if (!lists(res)) {
            if (!text.empty() && text.back() != '\n') text += '\n';
            text += res + "\tProtein\n";
            added += (added.empty() ? "" : ", ") + res;
        }
    }
    if (!added.empty() || !std::filesystem::exists(dst)) {
        write_text(dst, text);
        if (!added.empty()) std::cerr << "residuetypes.dat: added " << added << " as Protein -> " << dst << '\n';
    }
}

void Workflow::check_inputs() {
    if (sh_.dryRun()) return;

    if (cfg_.stage != "all" && cfg_.stage != "cg" && cfg_.stage != "prepare" &&
        cfg_.stage != "relax") return;

    if (!std::filesystem::exists(cfg_.pdbPath()))
        throw std::runtime_error("Input PDB not found: " + cfg_.pdbPath().string());

    const auto modified = detect_modified_residue(cfg_.pdbPath());
    if (!modified.empty() && !cfg_.is_phospho) {
        throw std::runtime_error(
            "Input PDB " + cfg_.pdbPath().string() + " contains modified residue '" + modified +
            "' but --phospho was not set. Without it the residue is written to martinize2 unchanged, "
            "is not recognized against the base force field, and its whole molecule is silently "
            "dropped. Re-run with --phospho, which renames it to its unmodified parent in a throwaway "
            "copy and supplies the martini3IDP force-field/mapping extensions this residue needs.");
    }

    if (cfg_.seq_length > 0) {
        const int actual = count_pdb_residues(cfg_.pdbPath());
        if (actual != cfg_.seq_length) {
            throw std::runtime_error(
                "Input PDB " + cfg_.pdbPath().string() + " contains " + std::to_string(actual) +
                " residue(s), but --seq-length is " + std::to_string(cfg_.seq_length) +
                ". Fix --seq-length (and --ss-string, if set, to the same length) to match the "
                "PDB's actual residue count - check whether a terminal residue is missing from the "
                "deposited structure (e.g. via the martinize2 'Missing atom ...:OXT' message, which "
                "names the residue actually treated as the C-terminus).");
        }
    }

    if (const auto chains = pdb_chain_ids(cfg_.pdbPath());
        chains.size() > 1 && cfg_.merge_chains.empty()) {
        throw std::runtime_error(
            "Input PDB " + cfg_.pdbPath().string() + " contains " + std::to_string(chains.size()) +
            " chains (" + join(chains, ",") + ") but --merge-chains was not given. This pipeline "
            "treats one 'protomer' as one CG molecule: pass --merge-chains " + join(chains, ",") +
            " to coarse-grain the chains into a single molecule (the right choice for a coiled-coil "
            "dimer, whose two chains are one physical unit), so that --n-prot then counts DIMERS and "
            "every inter-protomer descriptor becomes inter-dimer. Alternatively supply a "
            "single-chain PDB.");
    }

    if (cfg_.use_elastic && !trim(cfg_.elastic_units).empty()) {
        const auto [lo, hi] = pdb_residue_number_range(cfg_.pdbPath());
        std::vector<std::string> outside;
        std::istringstream units(trim(cfg_.elastic_units));
        for (std::string unit; std::getline(units, unit, ',');) {
            const auto colon = unit.find(':');
            if (colon == std::string::npos) continue;
            try {
                const int a = std::stoi(trim(unit.substr(0, colon)));
                const int b = std::stoi(trim(unit.substr(colon + 1)));
                if (a < lo || b > hi) outside.push_back(trim(unit));
            } catch (const std::exception&) {
                continue;
            }
        }
        if (!outside.empty()) {
            throw std::runtime_error(
                "--elastic-units " + trim(cfg_.elastic_units) + " names residue(s) outside " +
                cfg_.pdbPath().string() + ", whose residues run " + std::to_string(lo) + "-" +
                std::to_string(hi) + ": " + join(outside, ", ") + ".\n"
                "martinize2 reads -eunit as PDB residue numbers and quietly drops whatever falls "
                "outside the structure, so a range written in construct-local numbering does not "
                "fail - it builds a smaller elastic network over the wrong residues, leaving parts "
                "of a coiled-coil unrestrained while the topology still looks normal." +
                (lo > 1 ? " If these ranges are local positions, add " + std::to_string(lo - 1) +
                          " to each." : std::string{}));
        }
    }

    if (cfg_.martiniMainItp() != "martini_v3.0.0.itp") {
        const auto itp = cfg_.project_dir / "martini_v300" / cfg_.martiniMainItp();
        if (!std::filesystem::exists(itp))
            throw std::runtime_error(
                "--martini-lambda-pw " + to_string_fixed(cfg_.martini_lambda_pw, 2) +
                " richiede " + itp.string() + ", che non esiste. Generalo con:\n"
                "  python3 tools/make_lambda_itp.py " + to_string_fixed(cfg_.martini_lambda_pw, 2));
    }

    if (cfg_.solvation_mode == "gromacs" && cfg_.stage != "cg" &&
        !std::filesystem::exists(cfg_.martiniWaterPath())) {
        throw std::runtime_error(
            "--solvation-mode gromacs needs a Martini water box at " +
            cfg_.martiniWaterPath().string() + ", which does not exist.\n"
            "Either point --martini-water-gro at one, or pass --solvation-mode insane, which "
            "builds the solvent itself and needs no such file.");
    }

    if (cfg_.is_phospho) {
        const auto ff_dir = cfg_.martinizeFfDirPath();
        const auto map_dir = cfg_.martinizeMapDirPath();
        if (!std::filesystem::is_directory(ff_dir))
            throw std::runtime_error(
                "Missing Martinize2 extra force-field directory for phosphorylated input: " + ff_dir.string() +
                "\nSet --martinize-ff-dir PATH, or place the SEP/TPO/PTR force-field block(s) there "
                "(see README.md).");
        if (!std::filesystem::is_directory(map_dir))
            throw std::runtime_error(
                "Missing Martinize2 extra mapping directory for phosphorylated input: " + map_dir.string() +
                "\nSet --martinize-map-dir PATH, or place the SEP/TPO/PTR mapping file(s) there "
                "(see README.md).");

        for (const auto& residue : detect_modified_residues(cfg_.pdbPath())) {
            if (std::filesystem::exists(map_dir / (residue + ".mapping"))) continue;
            throw std::runtime_error(
                "Input PDB " + cfg_.pdbPath().string() + " contains residue '" + residue +
                "', for which there is no " + (map_dir / (residue + ".mapping")).string() + ".\n"
                "Only the modified residues with a modification/mapping pair are parametrized; this "
                "one would be renamed to its unmodified parent and its extra atoms silently dropped, "
                "so the run would finish having simulated the UNMODIFIED residue.\n"
                "Add the pair (force_fields/charmm/modification.ff + "
                "force_fields/martini3IDP/modification.ff + " + residue +
                ".mapping, following SEP's), or use a structure without this residue.");
        }

        if (!std::filesystem::is_directory(ff_dir / "charmm"))
            throw std::runtime_error(
                "--martinize-ff-dir (" + ff_dir.string() + ") has no 'charmm' subdirectory. "
                "vermouth needs -ff-dir to be the PARENT of one subdirectory per force-field name "
                "(a 'charmm/' subdirectory for input-residue recognition, plus a subdirectory named "
                "after your -ff, e.g. 'martini3IDP/', for the output CG parameters), as siblings. "
                "If your SEP/TPO/PTR .ff files currently live directly under --martinize-ff-dir, "
                "point --martinize-ff-dir at their parent directory instead. Verify with: "
                "martinize2 -ff-dir " + ff_dir.string() + " -list-blocks | grep -i sep");
    }
}

} // namespace cg
