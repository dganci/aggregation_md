#include "Preparer.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "GroUtils.hpp"
#include "PdbUtils.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cg {

void Preparer::verify_bead_count() const {
    if (sh_.dryRun()) return;

    const auto actual = count_pdb_atoms(cfg_.cgPath());
    if (actual == cfg_.atoms_per_prot) return;

    throw std::runtime_error(
        "Coarse-grained bead count mismatch: martinize2 produced " + std::to_string(actual) +
        " bead(s) for protomer '" + cfg_.protomer_name + "' (" + cfg_.cgPath().string() +
        "), but --atoms-per-prot is " + std::to_string(cfg_.atoms_per_prot) +
        ". Set --atoms-per-prot " + std::to_string(actual) + " and re-run.\n"
        "--atoms-per-prot is what slices the packed system into per-protomer index groups, so a "
        "wrong value does not merely fail: it silently mixes beads from neighbouring chains into "
        "each 'protomer' group, and every COM, contact number and collective variable built on "
        "them becomes meaningless. It is checked here, right after coarse-graining, rather than "
        "only in generate_index(): by then the run has already paid for the relaxation, Packmol, "
        "solvation and energy minimization.");
}

namespace {

std::set<std::string> itp_bead_types(const std::filesystem::path& itp) {
    std::set<std::string> beads;
    std::string section;
    for (const auto& line : read_lines(itp)) {
        const auto s = trim(line);
        if (starts_with(s, "[")) {
            section = trim(s.substr(1, s.find(']') == std::string::npos ? std::string::npos
                                                                       : s.find(']') - 1));
            continue;
        }
        if (section != "atoms" || s.empty() || starts_with(s, ";")) continue;
        std::istringstream in(s);
        std::string index, type;
        if (in >> index >> type) beads.insert(type);
    }
    return beads;
}

std::set<std::string> lambda_covered_types(const std::filesystem::path& itp) {
    static const std::string kMarker = "; bead proteici coperti:";
    std::ifstream in(itp);
    std::string line;
    for (int i = 0; i < 20 && std::getline(in, line); ++i) {
        if (!starts_with(line, kMarker)) continue;
        std::istringstream fields(line.substr(kMarker.size()));
        std::set<std::string> covered;
        for (std::string type; fields >> type;) covered.insert(type);
        return covered;
    }
    return {};
}

} // namespace

void Preparer::verify_lambda_coverage() const {
    if (sh_.dryRun()) return;
    if (cfg_.martiniMainItp() == "martini_v3.0.0.itp") return;

    const auto lambda_itp = cfg_.project_dir / "martini_v300" / cfg_.martiniMainItp();
    const auto protomer_itp = cfg_.prepDir() / (cfg_.protomer_name + "_0.itp");
    if (!std::filesystem::exists(protomer_itp)) return;

    const auto covered = lambda_covered_types(lambda_itp);
    if (covered.empty()) {
        throw std::runtime_error(
            lambda_itp.string() + " has no '; bead proteici coperti:' header, so which bead "
            "types it rescales cannot be verified. Regenerate it:\n"
            "  python3 tools/make_lambda_itp.py " + to_string_fixed(cfg_.martini_lambda_pw, 2));
    }

    std::vector<std::string> missing;
    for (const auto& bead : itp_bead_types(protomer_itp)) {
        if (covered.count(bead) == 0) missing.push_back(bead);
    }
    if (missing.empty()) return;

    throw std::runtime_error(
        "Incomplete lambda_PW rescaling: protomer '" + cfg_.protomer_name + "' uses bead type(s) " +
        join(missing, ", ") + ", which " + lambda_itp.string() + " does not rescale against water.\n"
        "This would not fail - the run would simply proceed with those bead types at their stock "
        "protein-water strength while every other type is at lambda " +
        to_string_fixed(cfg_.martini_lambda_pw, 2) + ", i.e. this construct would use a different "
        "force field from the others in the campaign, and the aggregation propensities would no "
        "longer be comparable. Regenerate the file over every construct:\n"
        "  python3 tools/make_lambda_itp.py " + to_string_fixed(cfg_.martini_lambda_pw, 2));
}

void Preparer::verify_periodic_margin() const {
    if (sh_.dryRun()) return;

    const auto packed = cfg_.packCleanPath();
    const auto solvated = std::filesystem::exists(cfg_.solvatedPath()) ? cfg_.solvatedPath()
                                                                      : cfg_.beadsPath();
    if (!std::filesystem::exists(packed) || !std::filesystem::exists(solvated)) return;

    const double dmax_nm = pdb_max_molecule_dmax(packed, cfg_.atoms_per_prot) / 10.0;
    if (dmax_nm <= 0.0) return;

    const double image = min_image_distance_nm(read_gro_box_vectors(solvated));
    const double needed = dmax_nm + 2.0 * cfg_.rcoulomb_nm;
    if (periodic_margin_ok(dmax_nm, image, cfg_.rcoulomb_nm)) return;

    throw std::runtime_error(
        "Periodic box too small for this protomer: its longest internal distance is " +
        to_string_fixed(dmax_nm, 1) + " nm and the nearest periodic image is " +
        to_string_fixed(image, 1) + " nm away, so parts of one chain sit within the " +
        to_string_fixed(cfg_.rcoulomb_nm, 2) + " nm cutoff of themselves (" +
        to_string_fixed(needed, 1) + " nm would leave one cutoff of headroom).\n"
        "A chain interacting with itself across the boundary produces contacts and an extended, "
        "pinned conformation that no log file reports as wrong, and that feed straight into the "
        "contact-number collective variable this study measures. Raise the margin:\n"
        "  --box-margin-nm " + to_string_fixed(cfg_.box_margin_nm + (needed - image) + 0.2, 1));
}


void Preparer::verify_pmf_wall() const {
    if (!cfg_.pmf || sh_.dryRun()) return;

    const auto solvated = std::filesystem::exists(cfg_.solvatedPath()) ? cfg_.solvatedPath()
                                                                      : cfg_.beadsPath();
    if (!std::filesystem::exists(solvated)) return;

    const double image = min_image_distance_nm(read_gro_box_vectors(solvated));
    if (image <= 0.0) return;
    const double limit = 0.5 * image - cfg_.rcoulomb_nm;
    if (cfg_.pmf_wall_nm <= limit) return;

    throw std::runtime_error(
        "PMF wall outside the usable range of this box: --pmf-wall-nm is " +
        to_string_fixed(cfg_.pmf_wall_nm, 2) + " nm but the nearest periodic image is " +
        to_string_fixed(image, 2) + " nm away, leaving " + to_string_fixed(limit, 2) +
        " nm once half the image distance and one " + to_string_fixed(cfg_.rcoulomb_nm, 2) +
        " nm cutoff are taken out.\n"
        "Past half the image distance the minimum-image convention folds a long COM separation "
        "onto a shorter one, so the sampled 'unbound' plateau is the lattice, not the dissociated "
        "dimer - and the standard-state correction read off that plateau is wrong by an amount "
        "nothing in the output reveals. Either lower the wall:\n"
        "  --pmf-wall-nm " + to_string_fixed(limit, 1) + "\n"
        "or enlarge the box:\n"
        "  --box-margin-nm " +
        to_string_fixed(cfg_.box_margin_nm + 2.0 * (cfg_.pmf_wall_nm - limit) + 0.2, 1));
}

} // namespace cg
