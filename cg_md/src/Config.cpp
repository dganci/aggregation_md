#include "Config.hpp"

#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cg {
namespace {

constexpr double ps_per_us = 1'000'000.0;

std::filesystem::path absolute_clean(const std::filesystem::path& p) {
    return std::filesystem::weakly_canonical(std::filesystem::absolute(p));
}

template <class T>
std::int64_t steps_from_us(T us, double dt_ps) {
    return static_cast<std::int64_t>(static_cast<double>(us) * ps_per_us / dt_ps + 0.5);
}

[[noreturn]] void help_and_exit() {
    std::cout
        << "Usage: cg_md [options]\n"
        << "\nCore:\n"
        << "  --project-dir PATH              --dry-run\n"
        << "  --stage all|prepare|em|nvt|npt|production|adaptive|metadynamics|center\n"
        << "  --gmx CMD                       --mkdssp CMD\n"
        << "  --martinize CMD                 --insane CMD                 --packmol CMD\n"
        << "  --plumed CMD                    --mpirun CMD                --python CMD\n"
        << "\nSystem:\n"
        << "  --n-prot N                      --protomer NAME\n"
        << "  --seq-length N                  --atoms-per-prot N            --phospho\n"
        << "  --temperature-K T               --salt-M M                    --ntomp N\n"
        << "\nPreparation:\n"
        << "  --no-elastic                    --elastic-units A:B\n"
        << "  --packmol-box-A A               --packmol-cluster-radius-A A\n"
        << "  --packmol-tolerance-A A         --packmol-radius-A A\n"
        << "  --solvation-mode gromacs|insane --box-type TYPE              --box-margin-nm NM\n"
        << "  --solvate-radius-nm NM          --martini-water-gro PATH      --no-ions\n"
        << "\nSimulation:\n"
        << "  --nvt-dt-ps PS                  --nvt-nsteps N\n"
        << "  --npt-dt-ps PS                  --npt-nsteps N\n"
        << "  --md-dt-ps PS                   --md-total-us US\n"
        << "  --plumed-stride N               --contact-r0-nm NM\n"
        << "  --thermostat system|protein-solvent|legacy\n"
        << "\nMetadynamics:\n"
        << "  --metadynamics                  --metad-walkers N             --metad-model PATH\n"
        << "  --metad-cv-params PATH          --metad-sigma CSV\n"
        << "  --metad-nodes N                 --metad-pace N\n"
        << "  --metad-height X                --metad-biasfactor X          --metad-print-stride N\n"
        << "  --no-metad-sum-hills            --no-metad-center\n"
        << "\nAdaptive:\n"
        << "  --adaptive                      --adaptive-chunk-us US         --adaptive-max-total-us US\n"
        << "  --adaptive-min-total-us US      --adaptive-min-cn-transitions N\n"
        << "  --adaptive-min-bidirectional-events N\n"
        << "  --adaptive-min-cn-range X       --adaptive-min-rg-range X\n"
        << "  --adaptive-min-contact-patterns N --adaptive-min-lcc-unique N\n"
        << "  --adaptive-pair-contact-threshold X --adaptive-smooth-window N\n"
        << "  --adaptive-min-residence N      --adaptive-pattern-downsample N\n";
    std::exit(0);
}

using Setter = std::function<void(Config&, const std::string&)>;

std::unordered_map<std::string, Setter> option_setters() {
    return {
        {"--project-dir", [](Config& c, const std::string& v) { c.project_dir = absolute_clean(v); }},
        {"--gmx", [](Config& c, const std::string& v) { c.gmx = v; }},
        {"--mkdssp", [](Config& c, const std::string& v) { c.mkdssp = v; }},
        {"--martinize", [](Config& c, const std::string& v) { c.martinize = v; }},
        {"--insane", [](Config& c, const std::string& v) { c.insane = v; }},
        {"--packmol", [](Config& c, const std::string& v) { c.packmol = v; }},
        {"--plumed", [](Config& c, const std::string& v) { c.plumed = v; }},
        {"--mpirun", [](Config& c, const std::string& v) { c.mpirun = v; }},
        {"--python", [](Config& c, const std::string& v) { c.python = v; }},
        {"--n-prot", [](Config& c, const std::string& v) { c.n_prot = std::stoi(v); }},
        {"--protomer", [](Config& c, const std::string& v) { c.protomer_name = v; }},
        {"--seq-length", [](Config& c, const std::string& v) { c.seq_length = std::stoi(v); }},
        {"--atoms-per-prot", [](Config& c, const std::string& v) { c.atoms_per_prot = std::stoi(v); }},
        {"--md-dt-ps", [](Config& c, const std::string& v) { c.md_dt_ps = std::stod(v); }},
        {"--md-total-us", [](Config& c, const std::string& v) { c.md_total_us = std::stod(v); }},
        {"--plumed-stride", [](Config& c, const std::string& v) { c.plumed_stride = std::stoi(v); }},
        {"--stage", [](Config& c, const std::string& v) { c.stage = v; }},
        {"--ntomp", [](Config& c, const std::string& v) { c.ntomp = std::stoi(v); }},
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
        {"--martini-water-gro", [](Config& c, const std::string& v) { c.martini_water_gro = v; }},
        {"--martinize-ff-dir", [](Config& c, const std::string& v) { c.martinize_ff_dir = v; }},
        {"--martinize-map-dir", [](Config& c, const std::string& v) { c.martinize_map_dir = v; }},
        {"--temperature-K", [](Config& c, const std::string& v) { c.temperature_K = std::stod(v); }},
        {"--salt-M", [](Config& c, const std::string& v) { c.salt_M = std::stod(v); }},
        {"--nvt-dt-ps", [](Config& c, const std::string& v) { c.nvt_dt_ps = std::stod(v); }},
        {"--nvt-nsteps", [](Config& c, const std::string& v) { c.nvt_nsteps = std::stoll(v); }},
        {"--npt-dt-ps", [](Config& c, const std::string& v) { c.npt_dt_ps = std::stod(v); }},
        {"--npt-nsteps", [](Config& c, const std::string& v) { c.npt_nsteps = std::stoll(v); }},
        {"--contact-r0-nm", [](Config& c, const std::string& v) { c.contact_r0_nm = std::stod(v); }},
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
        {"--adaptive-smooth-window", [](Config& c, const std::string& v) { c.adaptive_smooth_window_frames = std::stoi(v); }},
        {"--adaptive-min-residence", [](Config& c, const std::string& v) { c.adaptive_min_residence_frames = std::stoi(v); }},
        {"--adaptive-pattern-downsample", [](Config& c, const std::string& v) { c.adaptive_pattern_downsample = std::stoi(v); }}
    };
}

void require_one_of(const std::string& name, const std::string& value, const std::vector<std::string>& allowed) {
    for (const auto& x : allowed) if (value == x) return;
    throw std::runtime_error("Invalid " + name + ": " + value);
}

} // namespace

std::string Config::systemName() const { return std::to_string(n_prot) + "x" + protomer_name; }
std::filesystem::path Config::pdbPath() const { return project_dir / "pdbs" / "AA" / (protomer_name + ".pdb"); }
std::filesystem::path Config::cgPath() const { return project_dir / "pdbs" / "CG" / (protomer_name + "_cg.pdb"); }
std::filesystem::path Config::topologyPath() const { return project_dir / "system" / "topologies" / (protomer_name + "_topo.top"); }
std::filesystem::path Config::packRawPath() const { return project_dir / "pdbs" / "PACKMOL" / "raw" / (systemName() + ".pdb"); }
std::filesystem::path Config::packCleanPath() const { return project_dir / "pdbs" / "PACKMOL" / "clean" / (systemName() + ".pdb"); }
std::filesystem::path Config::packInputPath() const { return project_dir / "pdbs" / "PACKMOL" / "inps" / (systemName() + ".inp"); }
std::filesystem::path Config::boxedPath() const { return project_dir / "system" / "beads" / (systemName() + "_boxed.gro"); }
std::filesystem::path Config::solvatedPath() const { return project_dir / "system" / "beads" / (systemName() + "_solvated.gro"); }
std::filesystem::path Config::ionsTprPath() const { return resultDir() / "ions.tpr"; }
std::filesystem::path Config::martiniWaterPath() const { return martini_water_gro.is_absolute() ? martini_water_gro : project_dir / martini_water_gro; }
std::filesystem::path Config::martinizeFfDirPath() const { return martinize_ff_dir.is_absolute() ? martinize_ff_dir : project_dir / martinize_ff_dir; }
std::filesystem::path Config::martinizeMapDirPath() const { return martinize_map_dir.is_absolute() ? martinize_map_dir : project_dir / martinize_map_dir; }
std::filesystem::path Config::beadsPath() const { return project_dir / "system" / "beads" / (systemName() + "_beads.gro"); }
std::filesystem::path Config::systemDir() const { return project_dir / "system"; }
std::filesystem::path Config::resultDir() const { return project_dir / "runs" / systemName(); }
std::filesystem::path Config::metadDir() const { return resultDir() / "metadynamics"; }
std::filesystem::path Config::metadModelPath() const {
    if (!metad_model.empty()) return metad_model.is_absolute() ? metad_model : project_dir / metad_model;
    return metadDir() / "CVs" / protomer_name / "CVs_torchscript.pt";
}

std::filesystem::path Config::metadCvParamsPath() const {
    if (!metad_cv_params.empty()) return metad_cv_params.is_absolute() ? metad_cv_params : project_dir / metad_cv_params;
    return metadDir() / "CVs" / protomer_name / "cv_params.pkl";
}

std::int64_t Config::mdNsteps() const { return steps_from_us(md_total_us, md_dt_ps); }
std::int64_t Config::adaptiveChunkNsteps() const { return steps_from_us(adaptive_chunk_us, md_dt_ps); }
int Config::adaptiveMaxChunks() const { return adaptive_chunk_us > 0.0 ? static_cast<int>(adaptive_max_total_us / adaptive_chunk_us + 0.999999) : 0; }

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
        if (arg == "--no-elastic") { cfg.use_elastic = false; continue; }
        if (arg == "--no-ions") { cfg.add_ions = false; continue; }
        if (arg == "--adaptive") { cfg.adaptive = true; continue; }
        if (arg == "--metadynamics") { cfg.metadynamics = true; continue; }
        if (arg == "--no-metad-sum-hills") { cfg.metad_sum_hills = false; continue; }
        if (arg == "--no-metad-center") { cfg.metad_center = false; continue; }

        const auto it = setters.find(arg);
        if (it == setters.end()) throw std::runtime_error("Unknown argument: " + arg);
        it->second(cfg, value(i));
    }

    cfg.project_dir = absolute_clean(cfg.project_dir);
    if (cfg.stage == "metad" || cfg.stage == "metadynamics") cfg.metadynamics = true;
    cfg.validate();
    return cfg;
}

void Config::validate() const {
    require_one_of("--solvation-mode", solvation_mode, {"gromacs", "insane"});
    require_one_of("--thermostat", thermostat_mode, {"system", "protein-solvent", "legacy"});
    require_one_of("--stage", stage, {"all", "prepare", "em", "nvt", "npt", "production", "adaptive", "metadynamics", "metad", "fes", "center"});

    if (adaptive && metadynamics) throw std::runtime_error("--adaptive and --metadynamics are mutually exclusive production modes");

    if (n_prot <= 0) throw std::runtime_error("--n-prot must be > 0");
    if (seq_length <= 0) throw std::runtime_error("--seq-length must be > 0");
    if (atoms_per_prot <= 0) throw std::runtime_error("--atoms-per-prot must be > 0");
    if (box_margin_nm <= 0.0) throw std::runtime_error("--box-margin-nm must be > 0");
    if (md_dt_ps <= 0.0 || nvt_dt_ps <= 0.0 || npt_dt_ps <= 0.0) throw std::runtime_error("time steps must be > 0");
    if (plumed_stride <= 0) throw std::runtime_error("--plumed-stride must be > 0");
    if (metadynamics) {
        if (metad_walkers <= 0) throw std::runtime_error("--metad-walkers must be > 0");
        if (metad_nodes <= 0) throw std::runtime_error("--metad-nodes must be > 0");
        if (metad_pace <= 0) throw std::runtime_error("--metad-pace must be > 0");
        if (metad_print_stride <= 0) throw std::runtime_error("--metad-print-stride must be > 0");
    }
}

void Config::print() const {
    std::cerr
        << "Configuration\n"
        << "  project_dir       = " << project_dir << '\n'
        << "  system            = " << systemName() << '\n'
        << "  n_prot            = " << n_prot << '\n'
        << "  protomer_name     = " << protomer_name << '\n'
        << "  seq_length        = " << seq_length << '\n'
        << "  atoms_per_prot    = " << atoms_per_prot << '\n'
        << "  md_dt_ps          = " << md_dt_ps << '\n'
        << "  md_total_us       = " << md_total_us << '\n'
        << "  md_nsteps         = " << mdNsteps() << '\n'
        << "  plumed_stride     = " << plumed_stride << '\n'
        << "  thermostat_mode   = " << thermostat_mode << '\n'
        << "  solvation_mode    = " << solvation_mode << '\n'
        << "  packmol_R_A       = " << packmol_cluster_radius_A << '\n'
        << "  box_type          = " << box_type << '\n'
        << "  box_margin_nm     = " << box_margin_nm << '\n'
        << "  solvate_radius_nm = " << solvate_radius_nm << '\n'
        << "  martini_water_gro = " << martiniWaterPath() << '\n'
        << "  martinize_ff_dir  = " << martinizeFfDirPath() << '\n'
        << "  martinize_map_dir = " << martinizeMapDirPath() << '\n'
        << "  add_ions          = " << (add_ions ? "true" : "false") << '\n'
        << "  pdb_path          = " << pdbPath() << '\n'
        << "  topology_path     = " << topologyPath() << '\n'
        << "  result_dir        = " << resultDir() << '\n'
        << "  adaptive          = " << (adaptive ? "true" : "false") << '\n'
        << "  adaptive_chunk_us = " << adaptive_chunk_us << '\n'
        << "  adaptive_max_us   = " << adaptive_max_total_us << '\n'
        << "  metadynamics      = " << (metadynamics ? "true" : "false") << '\n';
    if (metadynamics) {
        std::cerr
            << "  metad_dir         = " << metadDir() << '\n'
            << "  metad_model       = " << metadModelPath() << '\n'
            << "  metad_cv_params   = " << metadCvParamsPath() << '\n'
            << "  metad_walkers     = " << metad_walkers << '\n'
            << "  metad_sigma       = " << (metad_sigma.empty() ? "<from cv_params.pkl>" : metad_sigma) << '\n';
    }
}

} // namespace cg
