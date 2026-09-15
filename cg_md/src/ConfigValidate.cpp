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

void require_one_of(const std::string& name, const std::string& value, const std::vector<std::string>& allowed) {
    for (const auto& x : allowed) if (value == x) return;
    throw std::runtime_error("Invalid " + name + ": " + value);
}

} // namespace

void Config::validate() const {
    require_one_of("--solvation-mode", solvation_mode, {"gromacs", "insane"});
    require_one_of("--thermostat", thermostat_mode, {"system", "protein-solvent", "legacy"});
    require_one_of("--stage", stage, {"all", "cg", "prepare", "em", "nvt", "npt", "production", "adaptive",
                                      "metadynamics", "metad", "pmf", "fes", "center", "report",
                                      "relax", "contact-map"});

    if (adaptive && metadynamics) throw std::runtime_error("--adaptive and --metadynamics are mutually exclusive production modes");

    if (use_elastic && elastic_units == kDefaultElasticUnits)
        throw std::runtime_error(
            "--elastic-units is still the built-in default (" + kDefaultElasticUnits + "), which "
            "belongs to no construct in particular. An elastic network over the wrong residues "
            "does not fail - it restrains the wrong part of the chain and leaves the rest free.\n"
            "Derive the ranges with  python3 tools/ss_from_pdb.py <pdb>  and pass them, or pass "
            "--no-elastic for a disordered construct.");

    if (n_prot < 2 && stage != "relax" && stage != "cg")
        throw std::runtime_error("--n-prot must be >= 2: every descriptor and sampling metric in this "
                                 "pipeline is defined between protomer pairs, so a single copy has "
                                 "nothing to measure. Use a plain GROMACS run for a single-chain control.");
    if (seq_length <= 0) throw std::runtime_error("--seq-length must be > 0");
    if (!ss_string.empty() && ss_string.size() != static_cast<std::size_t>(seq_length)) {
        throw std::runtime_error("--ss-string length (" + std::to_string(ss_string.size()) +
                                 ") must equal --seq-length (" + std::to_string(seq_length) + ")");
    }
    if (atoms_per_prot <= 0) throw std::runtime_error("--atoms-per-prot must be > 0");
    if (box_margin_nm <= 0.0) throw std::runtime_error("--box-margin-nm must be > 0");
    if (martini_lambda_pw < 0.5 || martini_lambda_pw > 2.0)
        throw std::runtime_error("--martini-lambda-pw must be in [0.5, 2.0]");
    if (md_dt_ps <= 0.0 || nvt_dt_ps <= 0.0 || npt_dt_ps <= 0.0) throw std::runtime_error("time steps must be > 0");
    if (plumed_stride <= 0) throw std::runtime_error("--plumed-stride must be > 0");
    if (contact_nn <= 0 || contact_mm <= contact_nn)
        throw std::runtime_error("--contact-mm must be > --contact-nn > 0");
    if (epsilon_r <= 0.0) throw std::runtime_error("--epsilon-r must be > 0");
    if (rcoulomb_nm <= 0.0 || rvdw_nm <= 0.0) throw std::runtime_error("--rcoulomb-nm/--rvdw-nm must be > 0");
    if (nstlist <= 0) throw std::runtime_error("--nstlist must be > 0");
    if (verlet_buffer_tolerance != -1.0 && verlet_buffer_tolerance <= 0.0)
        throw std::runtime_error("--verlet-buffer-tolerance must be > 0, or exactly -1 to fix rlist");
    if (verlet_buffer_tolerance == -1.0 && rlist_nm < std::max(rcoulomb_nm, rvdw_nm))
        throw std::runtime_error("--rlist-nm must be >= the largest cut-off: with "
                                 "--verlet-buffer-tolerance -1 it IS the pair list radius, and a "
                                 "list shorter than the cut-off misses interactions outright");
    require_one_of("--pin", pin, {"auto", "on", "off"});
    if (pin != "on" && (pinoffset >= 0 || pinstride >= 0))
        throw std::runtime_error("--pinoffset/--pinstride only take effect with --pin on");
    if (nst_xtc <= 0 || nst_energy <= 0 || nst_log <= 0)
        throw std::runtime_error("--nst-xtc/--nst-energy/--nst-log must be > 0");
    if (xtc_precision <= 0) throw std::runtime_error("--xtc-precision must be > 0");

    if (adaptive) {
        if (adaptive_chunk_us <= 0.0) throw std::runtime_error("--adaptive-chunk-us must be > 0");
        if (adaptive_max_total_us < adaptive_chunk_us)
            throw std::runtime_error("--adaptive-max-total-us must be >= --adaptive-chunk-us");
        if (adaptive_min_total_us > adaptive_max_total_us)
            throw std::runtime_error("--adaptive-min-total-us must be <= --adaptive-max-total-us");
        if (adaptiveChunkNsteps() <= 0)
            throw std::runtime_error("--adaptive-chunk-us is shorter than one --md-dt-ps step");
        if (adaptive_max_pattern_growth < 1.0)
            throw std::runtime_error("--adaptive-max-pattern-growth is a ratio of a total to a part and cannot be < 1");
        if (adaptive_max_pattern_jsd < 0.0 || adaptive_max_pattern_jsd > 1.0)
            throw std::runtime_error("--adaptive-max-pattern-jsd is in bits and lies in [0, 1]");
        if (adaptive_min_effective_samples < 0 || adaptive_min_assembly_events < 0)
            throw std::runtime_error("--adaptive-min-effective-samples/--adaptive-min-assembly-events must be >= 0");
        if (adaptive_event_residence_ps < 0.0)
            throw std::runtime_error("--adaptive-event-residence-ps must be >= 0");
        if (adaptive_min_time_over_its < 0.0)
            throw std::runtime_error("--adaptive-min-time-over-its must be >= 0 (0 turns the gate off)");
    }

    if (metadynamics) {
        if (metad_walkers <= 0) throw std::runtime_error("--metad-walkers must be > 0");
        if (metad_nodes <= 0) throw std::runtime_error("--metad-nodes must be > 0");
        if (metad_pace <= 0) throw std::runtime_error("--metad-pace must be > 0");
        if (metad_print_stride <= 0) throw std::runtime_error("--metad-print-stride must be > 0");
        if (metad_biasfactor <= 1.0)
            throw std::runtime_error("--metad-biasfactor must be > 1 for well-tempered metadynamics");
        if (metad_height <= 0.0) throw std::runtime_error("--metad-height must be > 0");
        if (metad_grid_bin <= 0) throw std::runtime_error("--metad-grid-bin must be > 0");
        if (metad_wall_kt < 0.0) throw std::runtime_error("--metad-wall-kt must be >= 0 (0 disables the walls)");
        if (metad_wall_margin_frac <= 0.0 || metad_wall_margin_frac >= 0.5)
            throw std::runtime_error("--metad-wall-margin-frac must be in (0, 0.5): the walls sit that "
                                     "fraction of the grid span inside each bound, so at 0.5 they would "
                                     "meet in the middle of the grid.");
        if (metad_rct_ustride <= 0) throw std::runtime_error("--metad-rct-ustride must be > 0");
        if (metadTotalUs() <= 0.0)
            throw std::runtime_error("Metadynamics needs a positive --metad-total-us (or --md-total-us)");
        if (metad_chunk_us < 0.0) throw std::runtime_error("--metad-chunk-us must be >= 0");
        if (metad_chunk_us > 0.0 && metadChunkNsteps() <= 0)
            throw std::runtime_error("--metad-chunk-us is shorter than one --md-dt-ps step");
        if (!metad_grid && metad_calc_rct)
            throw std::runtime_error("--no-metad-grid is incompatible with c(t) reweighting; "
                                     "add --no-metad-rct if you really want a gridless run "
                                     "(you will not be able to reweight it)");
    }

    if (distinct_copies && !relax)
        throw std::runtime_error(
            "--distinct-copies needs --relax: the per-copy conformations are frames of the "
            "relaxation trajectory, and without one every copy would silently be the same "
            "structure. Add --relax (with --relax-us), or drop --distinct-copies.");

    if (relax) {
        if (relax_us <= 0.0) throw std::runtime_error("--relax-us must be > 0");
        if (relax_cluster_cutoff_nm <= 0.0) throw std::runtime_error("--relax-cluster-cutoff-nm must be > 0");
        if (relax_discard_frac < 0.0 || relax_discard_frac >= 1.0)
            throw std::runtime_error("--relax-discard-frac must be in [0, 1)");
        if (relax_frames < 10)
            throw std::runtime_error("--relax-frames must be >= 10: fewer frames than that cannot "
                                     "be clustered into anything meaningful");
    }

    if (pmf) {
        if (n_prot != 2)
            throw std::runtime_error("--pmf requires exactly --n-prot 2: the potential of mean force "
                                     "is defined between one pair of chains, and its standard-state "
                                     "correction assumes a two-body system. Use the metadynamics "
                                     "stage for many-chain aggregation.");
        if (pmf_wall_nm <= 0.0) throw std::runtime_error("--pmf-wall-nm must be > 0");
        if (pmf_wall_kappa <= 0.0) throw std::runtime_error("--pmf-wall-kappa must be > 0");
        if (pmf_sigma_nm <= 0.0) throw std::runtime_error("--pmf-sigma-nm must be > 0");
        if (pmf_bound_cutoff_nm < 0.0) throw std::runtime_error("--pmf-bound-cutoff-nm must be >= 0");
        if (pmf_plateau_frac <= 0.0 || pmf_plateau_frac >= 0.95)
            throw std::runtime_error("--pmf-plateau-frac must be in (0, 0.95)");
        const double estimated_image_nm = box_margin_nm + 2.0 * packmol_cluster_radius_A / 10.0;
        if (pmf_wall_nm > 0.5 * estimated_image_nm)
            std::cerr << "Warning: --pmf-wall-nm (" << pmf_wall_nm << " nm) looks larger than half the "
                         "expected minimum-image distance (" << to_string_fixed(estimated_image_nm, 1)
                      << " nm). Beyond that a COM separation is aliased by the minimum-image "
                         "convention and the PMF's unbound plateau is meaningless.\n";
    }
}

} // namespace cg
