#include "Reporter.hpp"
#include "Config.hpp"
#include "Shell.hpp"

#include <exception>
#include <iostream>
#include <vector>

namespace cg {
namespace {

std::filesystem::path executable_relative_script() {
#ifdef __linux__
    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec) {
        const auto candidate = exe.parent_path().parent_path() / "scripts" / "analyze_run.py";
        if (std::filesystem::exists(candidate)) return candidate;
    }
#endif
    return {};
}

} // namespace

std::filesystem::path resolve_sibling_script(const Config& cfg, const char* name) {
    const auto in_project = cfg.project_dir / "scripts" / name;
    if (std::filesystem::exists(in_project)) return in_project;
#ifdef __linux__
    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec) {
        const auto candidate = exe.parent_path().parent_path() / "scripts" / name;
        if (std::filesystem::exists(candidate)) return candidate;
    }
#endif
    return {};
}

namespace {
std::filesystem::path sibling_script(const Config& cfg, const char* name) {
    return resolve_sibling_script(cfg, name);
}
} // namespace

void run_contact_map(const Config& cfg,
                     const Shell& sh,
                     const std::filesystem::path& traj,
                     const std::filesystem::path& top,
                     const std::filesystem::path& out_dir) {
    if (!cfg.report) return;
    try {
        const auto script = sibling_script(cfg, "contact_map.py");
        if (script.empty()) return;
        if (!sh.dryRun() && (!std::filesystem::exists(traj) || !std::filesystem::exists(top))) return;

        const int rc = sh.run({cfg.python, script.string(),
                               "--traj", traj.string(), "--top", top.string(),
                               "--n-prot", std::to_string(cfg.n_prot),
                               "--cutoff-nm", std::to_string(cfg.contact_map_cutoff_nm),
                               "--stride", std::to_string(cfg.contact_map_stride),
                               "--out-dir", out_dir.string(),
                               "--label", cfg.systemName()},
                              {}, /*check=*/false);
        if (rc != 0)
            std::cerr << "Warning: contact-map analysis exited with code " << rc
                      << " (ignored, purely observational).\n";
    } catch (const std::exception& e) {
        std::cerr << "Warning: contact-map analysis failed (ignored): " << e.what() << '\n';
    }
}

std::filesystem::path resolve_report_script(const Config& cfg) {
    if (!cfg.report_script.empty()) {
        const auto explicit_path = cfg.report_script.is_absolute() ? cfg.report_script
                                                                   : cfg.project_dir / cfg.report_script;
        return std::filesystem::exists(explicit_path) ? explicit_path : std::filesystem::path{};
    }

    const auto in_project = cfg.project_dir / "scripts" / "analyze_run.py";
    if (std::filesystem::exists(in_project)) return in_project;

    return executable_relative_script();
}

void run_report(const Config& cfg,
                const Shell& sh,
                const std::filesystem::path& run_dir,
                const std::string& mode_in) {
    if (!cfg.report) return;
    const std::string mode = (cfg.pmf && mode_in == "metadynamics") ? "pmf" : mode_in;

    try {
        const auto script = resolve_report_script(cfg);
        if (script.empty()) {
            std::cerr << "Note: skipping report (scripts/analyze_run.py not found; "
                         "pass --report-script PATH or --no-report to silence this).\n";
            return;
        }

        std::vector<std::string> cmd = {
            cfg.python, script.string(),
            "--run-dir", run_dir.string(),
            "--mode", mode,
            "--n-prot", std::to_string(cfg.n_prot),
            "--dt-ps", std::to_string(cfg.md_dt_ps),
            "--plumed-stride", std::to_string(cfg.plumed_stride),
            "--temperature-K", std::to_string(cfg.temperature_K),
            "--contact-threshold", std::to_string(cfg.adaptive_pair_contact_threshold),
            "--biasfactor", std::to_string(cfg.metad_biasfactor)
        };

        if (cfg.pmf) {
            cmd.insert(cmd.end(), {
                "--pmf-wall-nm", std::to_string(cfg.pmf_wall_nm),
                "--pmf-bound-cutoff-nm", std::to_string(cfg.pmf_bound_cutoff_nm),
                "--pmf-plateau-frac", std::to_string(cfg.pmf_plateau_frac)
            });
        }

        const int rc = sh.run(cmd, {}, /*check=*/false);
        if (rc != 0)
            std::cerr << "Warning: report generation exited with code " << rc
                      << " (ignored, purely observational).\n";
    } catch (const std::exception& e) {
        std::cerr << "Warning: report generation failed (ignored, purely observational): "
                  << e.what() << '\n';
    }
}

} // namespace cg
