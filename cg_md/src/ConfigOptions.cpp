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

[[noreturn]] void help_and_exit();

namespace {

std::filesystem::path absolute_clean(const std::filesystem::path& p) {
    return std::filesystem::weakly_canonical(std::filesystem::absolute(p));
}


using Setter = std::function<void(Config&, const std::string&)>;

std::unordered_map<std::string, Setter> option_setters() {
    return {
        {"--project-dir", [](Config& c, const std::string& v) { c.project_dir = absolute_clean(v); }},
        {"--gmx", [](Config& c, const std::string& v) { c.gmx = v; }},
        {"--ss-string", [](Config& c, const std::string& v) { c.ss_string = v; }},
        {"--martinize", [](Config& c, const std::string& v) { c.martinize = v; }},
        {"--insane", [](Config& c, const std::string& v) { c.insane = v; }},
        {"--packmol", [](Config& c, const std::string& v) { c.packmol = v; }},
        {"--plumed", [](Config& c, const std::string& v) { c.plumed = v; }},
        {"--mpirun", [](Config& c, const std::string& v) { c.mpirun = v; }},
        {"--python", [](Config& c, const std::string& v) { c.python = v; }},
        {"--n-prot", [](Config& c, const std::string& v) { c.n_prot = std::stoi(v); }},
        {"--protomer", [](Config& c, const std::string& v) { c.protomer_name = v; }},
        {"--run-tag", [](Config& c, const std::string& v) { c.run_tag = v; }},
        {"--seq-length", [](Config& c, const std::string& v) { c.seq_length = std::stoi(v); }},
        {"--atoms-per-prot", [](Config& c, const std::string& v) { c.atoms_per_prot = std::stoi(v); }},
        {"--md-dt-ps", [](Config& c, const std::string& v) { c.md_dt_ps = std::stod(v); }},
        {"--md-total-us", [](Config& c, const std::string& v) { c.md_total_us = std::stod(v); }},
        {"--plumed-stride", [](Config& c, const std::string& v) { c.plumed_stride = std::stoi(v); }},
        {"--stage", [](Config& c, const std::string& v) { c.stage = v; }},
        {"--ntomp", [](Config& c, const std::string& v) { c.ntomp = std::stoi(v); }},
        {"--pin", [](Config& c, const std::string& v) { c.pin = v; }},
        {"--pinoffset", [](Config& c, const std::string& v) { c.pinoffset = std::stoi(v); }},
        {"--pinstride", [](Config& c, const std::string& v) { c.pinstride = std::stoi(v); }},
        {"--thermostat", [](Config& c, const std::string& v) { c.thermostat_mode = v; }},
        {"--elastic-units", [](Config& c, const std::string& v) { c.elastic_units = v; }},
        {"--packmol-box-A", [](Config& c, const std::string& v) { c.packmol_box_A = std::stod(v); c.packmol_cluster_radius_A = 0.5 * c.packmol_box_A; }},
        {"--packmol-cluster-radius-A", [](Config& c, const std::string& v) { c.packmol_cluster_radius_A = std::stod(v); c.packmol_box_A = 2.0 * c.packmol_cluster_radius_A; }},
        {"--packmol-tolerance-A", [](Config& c, const std::string& v) { c.packmol_tolerance_A = std::stod(v); }},
        {"--packmol-radius-A", [](Config& c, const std::string& v) { c.packmol_radius_A = std::stod(v); }},
        {"--solvation-mode", [](Config& c, const std::string& v) { c.solvation_mode = v; }},
        {"--box-type", [](Config& c, const std::string& v) { c.box_type = v; }},
        {"--box-margin-nm", [](Config& c, const std::string& v) { c.box_margin_nm = std::stod(v); }},
        {"--solvate-radius-nm", [](Config& c, const std::string& v) { c.solvate_radius_nm = std::stod(v); }},
        {"--martini-lambda-pw", [](Config& c, const std::string& v) { c.martini_lambda_pw = std::stod(v); }},
        {"--martini-water-gro", [](Config& c, const std::string& v) { c.martini_water_gro = v; }},
        {"--martinize-ff-dir", [](Config& c, const std::string& v) { c.martinize_ff_dir = v; }},
        {"--martinize-map-dir", [](Config& c, const std::string& v) { c.martinize_map_dir = v; }},
        {"--temperature-K", [](Config& c, const std::string& v) { c.temperature_K = std::stod(v); }},
        {"--salt-M", [](Config& c, const std::string& v) { c.salt_M = std::stod(v); }},
        {"--seed", [](Config& c, const std::string& v) { c.seed = std::stoi(v); }},
        {"--nvt-dt-ps", [](Config& c, const std::string& v) { c.nvt_dt_ps = std::stod(v); }},
        {"--nvt-nsteps", [](Config& c, const std::string& v) { c.nvt_nsteps = std::stoll(v); }},
        {"--npt-dt-ps", [](Config& c, const std::string& v) { c.npt_dt_ps = std::stod(v); }},
        {"--npt-nsteps", [](Config& c, const std::string& v) { c.npt_nsteps = std::stoll(v); }},
        {"--contact-r0-nm", [](Config& c, const std::string& v) { c.contact_r0_nm = std::stod(v); }},
        {"--epsilon-r", [](Config& c, const std::string& v) { c.epsilon_r = std::stod(v); }},
        {"--rcoulomb-nm", [](Config& c, const std::string& v) { c.rcoulomb_nm = std::stod(v); }},
        {"--rvdw-nm", [](Config& c, const std::string& v) { c.rvdw_nm = std::stod(v); }},
        {"--nstlist", [](Config& c, const std::string& v) { c.nstlist = std::stoi(v); }},
        {"--rlist-nm", [](Config& c, const std::string& v) { c.rlist_nm = std::stod(v); }},
        {"--verlet-buffer-tolerance", [](Config& c, const std::string& v) { c.verlet_buffer_tolerance = std::stod(v); }},
        {"--contact-nn", [](Config& c, const std::string& v) { c.contact_nn = std::stoi(v); }},
        {"--contact-mm", [](Config& c, const std::string& v) { c.contact_mm = std::stoi(v); }},
        {"--nst-xtc", [](Config& c, const std::string& v) { c.nst_xtc = std::stoi(v); }},
        {"--nst-energy", [](Config& c, const std::string& v) { c.nst_energy = std::stoi(v); }},
        {"--nst-log", [](Config& c, const std::string& v) { c.nst_log = std::stoi(v); }},
        {"--xtc-precision", [](Config& c, const std::string& v) { c.xtc_precision = std::stoi(v); }},
        {"--contact-map-cutoff-nm", [](Config& c, const std::string& v) { c.contact_map_cutoff_nm = std::stod(v); }},
        {"--contact-map-stride", [](Config& c, const std::string& v) { c.contact_map_stride = std::stoi(v); }},
        {"--report-script", [](Config& c, const std::string& v) { c.report_script = v; }},
        {"--metad-total-us", [](Config& c, const std::string& v) { c.metad_total_us = std::stod(v); }},
        {"--metad-chunk-us", [](Config& c, const std::string& v) { c.metad_chunk_us = std::stod(v); }},
        {"--metad-feature-cols", [](Config& c, const std::string& v) { c.metad_feature_cols = v; }},
        {"--metad-grid-min", [](Config& c, const std::string& v) { c.metad_grid_min = v; }},
        {"--metad-grid-max", [](Config& c, const std::string& v) { c.metad_grid_max = v; }},
        {"--metad-grid-bin", [](Config& c, const std::string& v) { c.metad_grid_bin = std::stoi(v); }},
        {"--metad-wall-kt", [](Config& c, const std::string& v) { c.metad_wall_kt = std::stod(v); }},
        {"--metad-wall-margin-frac", [](Config& c, const std::string& v) { c.metad_wall_margin_frac = std::stod(v); }},
        {"--metad-grid-wstride", [](Config& c, const std::string& v) { c.metad_grid_wstride = std::stoi(v); }},
        {"--metad-rct-ustride", [](Config& c, const std::string& v) { c.metad_rct_ustride = std::stoi(v); }},
        {"--metad-fes-stride", [](Config& c, const std::string& v) { c.metad_fes_stride = std::stoi(v); }},
        {"--relax-us", [](Config& c, const std::string& v) { c.relax_us = std::stod(v); }},
        {"--relax-cluster-cutoff-nm", [](Config& c, const std::string& v) { c.relax_cluster_cutoff_nm = std::stod(v); }},
        {"--relax-discard-frac", [](Config& c, const std::string& v) { c.relax_discard_frac = std::stod(v); }},
        {"--relax-frames", [](Config& c, const std::string& v) { c.relax_frames = std::stoi(v); }},
        {"--pmf-wall-nm", [](Config& c, const std::string& v) { c.pmf_wall_nm = std::stod(v); }},
        {"--pmf-wall-kappa", [](Config& c, const std::string& v) { c.pmf_wall_kappa = std::stod(v); }},
        {"--pmf-sigma-nm", [](Config& c, const std::string& v) { c.pmf_sigma_nm = std::stod(v); }},
        {"--pmf-bound-cutoff-nm", [](Config& c, const std::string& v) { c.pmf_bound_cutoff_nm = std::stod(v); }},
        {"--pmf-plateau-frac", [](Config& c, const std::string& v) { c.pmf_plateau_frac = std::stod(v); }},
        {"--metad-walkers", [](Config& c, const std::string& v) { c.metad_walkers = std::stoi(v); }},
        {"--metad-nodes", [](Config& c, const std::string& v) { c.metad_nodes = std::stoi(v); }},
        {"--metad-pace", [](Config& c, const std::string& v) { c.metad_pace = std::stoi(v); }},
        {"--metad-print-stride", [](Config& c, const std::string& v) { c.metad_print_stride = std::stoi(v); }},
        {"--metad-height", [](Config& c, const std::string& v) { c.metad_height = std::stod(v); }},
        {"--metad-biasfactor", [](Config& c, const std::string& v) { c.metad_biasfactor = std::stod(v); }},
        {"--metad-sigma", [](Config& c, const std::string& v) { c.metad_sigma = v; }},
        {"--metad-model", [](Config& c, const std::string& v) { c.metad_model = v; }},
        {"--metad-cv-params", [](Config& c, const std::string& v) { c.metad_cv_params = v; }},
        {"--adaptive-chunk-us", [](Config& c, const std::string& v) { c.adaptive_chunk_us = std::stod(v); }},
        {"--adaptive-max-total-us", [](Config& c, const std::string& v) { c.adaptive_max_total_us = std::stod(v); }},
        {"--adaptive-min-total-us", [](Config& c, const std::string& v) { c.adaptive_min_total_us = std::stod(v); }},
        {"--adaptive-min-cn-transitions", [](Config& c, const std::string& v) { c.adaptive_min_cn_state_transitions = std::stoi(v); }},
        {"--adaptive-min-bidirectional-events", [](Config& c, const std::string& v) { c.adaptive_min_bidirectional_events = std::stoi(v); }},
        {"--adaptive-min-cn-range", [](Config& c, const std::string& v) { c.adaptive_min_cn_range = std::stod(v); }},
        {"--adaptive-min-rg-range", [](Config& c, const std::string& v) { c.adaptive_min_rg_global_range = std::stod(v); }},
        {"--adaptive-min-contact-patterns", [](Config& c, const std::string& v) { c.adaptive_min_unique_contact_patterns = std::stoi(v); }},
        {"--adaptive-min-lcc-unique", [](Config& c, const std::string& v) { c.adaptive_min_largest_cluster_unique = std::stoi(v); }},
        {"--adaptive-pair-contact-threshold", [](Config& c, const std::string& v) { c.adaptive_pair_contact_threshold = std::stod(v); }},
        {"--adaptive-max-pattern-growth", [](Config& c, const std::string& v) { c.adaptive_max_pattern_growth = std::stod(v); }},
        {"--adaptive-max-pattern-jsd", [](Config& c, const std::string& v) { c.adaptive_max_pattern_jsd = std::stod(v); }},
        {"--adaptive-min-effective-samples", [](Config& c, const std::string& v) { c.adaptive_min_effective_samples = std::stoi(v); }},
        {"--adaptive-min-assembly-events", [](Config& c, const std::string& v) { c.adaptive_min_assembly_events = std::stoi(v); }},
        {"--adaptive-event-residence-ps", [](Config& c, const std::string& v) { c.adaptive_event_residence_ps = std::stod(v); }},
        {"--adaptive-min-time-over-its", [](Config& c, const std::string& v) { c.adaptive_min_time_over_its = std::stod(v); }},
        {"--adaptive-smooth-window", [](Config& c, const std::string& v) { c.adaptive_smooth_window_frames = std::stoi(v); }},
        {"--adaptive-min-residence", [](Config& c, const std::string& v) { c.adaptive_min_residence_frames = std::stoi(v); }},
        {"--adaptive-pattern-downsample", [](Config& c, const std::string& v) { c.adaptive_pattern_downsample = std::stoi(v); }}
    };
}

} // namespace

Config Config::fromArgs(int argc, char** argv) {
    Config cfg;
    const auto setters = option_setters();

    const auto value = [&](int& i) -> std::string {
        if (i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + argv[i]);
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") help_and_exit();
        if (arg == "--dry-run") { cfg.dry_run = true; continue; }
        if (arg == "--phospho") { cfg.is_phospho = true; continue; }
        if (arg == "--modify") { cfg.modify = value(i); continue; }
        if (arg == "--merge-chains") { cfg.merge_chains = value(i); continue; }
        if (arg == "--dssp") { cfg.use_dssp = true; continue; }
        if (arg == "--relax") { cfg.relax = true; continue; }
        if (arg == "--distinct-copies") { cfg.distinct_copies = true; continue; }
        if (arg == "--no-elastic") { cfg.use_elastic = false; continue; }
        if (arg == "--no-ions") { cfg.add_ions = false; continue; }
        if (arg == "--adaptive") { cfg.adaptive = true; continue; }
        if (arg == "--metadynamics") { cfg.metadynamics = true; continue; }
        if (arg == "--no-metad-sum-hills") { cfg.metad_sum_hills = false; continue; }
        if (arg == "--no-metad-center") { cfg.metad_center = false; continue; }
        if (arg == "--no-metad-grid") { cfg.metad_grid = false; continue; }
        if (arg == "--keep-intermediate-trajectories") { cfg.keep_intermediate_trajectories = true; continue; }
        if (arg == "--no-metad-rct") { cfg.metad_calc_rct = false; continue; }
        if (arg == "--no-metad-work") { cfg.metad_calc_work = false; continue; }
        if (arg == "--no-report") { cfg.report = false; continue; }
        if (arg == "--pmf") { cfg.pmf = true; continue; }

        const auto it = setters.find(arg);
        if (it == setters.end()) throw std::runtime_error("Unknown argument: " + arg);
        const std::string raw = value(i);
        try {
            it->second(cfg, raw);
        } catch (const std::invalid_argument&) {
            throw std::runtime_error("Not a number for " + arg + ": \"" + raw + "\"");
        } catch (const std::out_of_range&) {
            throw std::runtime_error("Value out of range for " + arg + ": \"" + raw + "\"");
        }
    }

    cfg.project_dir = absolute_clean(cfg.project_dir);
    if (cfg.stage == "metad" || cfg.stage == "metadynamics") cfg.metadynamics = true;
    if (cfg.stage == "adaptive") cfg.adaptive = true;
    if (cfg.stage == "pmf") cfg.pmf = true;
    if (cfg.stage == "relax") cfg.relax = true;

    if (cfg.pmf) {
        cfg.metadynamics = true;
        cfg.metad_nodes = 1;
        if (cfg.metad_sigma.empty()) cfg.metad_sigma = to_string_fixed(cfg.pmf_sigma_nm, 4);
        if (cfg.metad_grid_min.empty()) cfg.metad_grid_min = "0";
        if (cfg.metad_grid_max.empty())
            cfg.metad_grid_max = to_string_fixed(cfg.pmf_wall_nm + 1.0, 4);
        cfg.metad_feature_cols = "d_1_2";
    }

    cfg.validate();
    cfg.runId();
    return cfg;
}

} // namespace cg
