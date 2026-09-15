#include "AdaptiveSampler.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "GromacsDriver.hpp"
#include "IndexBuilder.hpp"
#include "MdpWriter.hpp"
#include "PlumedWriter.hpp"
#include "Reporter.hpp"
#include "SamplingMonitor.hpp"
#include "RunDiagnostics.hpp"
#include "RunLedger.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
namespace cg {

AdaptiveSampler::AdaptiveSampler(Config& cfg, Shell& sh, GromacsDriver& gmx) : cfg_(cfg), sh_(sh), gmx_(gmx) {}

SamplingRules AdaptiveSampler::rules() const {
    SamplingRules r;
    r.min_total_us = cfg_.adaptive_min_total_us;
    r.min_cn_state_transitions = cfg_.adaptive_min_cn_state_transitions;
    r.min_bidirectional_events = cfg_.adaptive_min_bidirectional_events;
    r.min_cn_range = cfg_.adaptive_min_cn_range;
    r.min_rg_global_range = cfg_.adaptive_min_rg_global_range;
    r.min_unique_contact_patterns = cfg_.adaptive_min_unique_contact_patterns;
    r.min_largest_cluster_unique = cfg_.adaptive_min_largest_cluster_unique;
    r.pair_contact_threshold = cfg_.adaptive_pair_contact_threshold;
    r.max_pattern_growth_ratio = cfg_.adaptive_max_pattern_growth;
    r.max_pattern_jsd_bits = cfg_.adaptive_max_pattern_jsd;
    r.min_effective_samples = cfg_.adaptive_min_effective_samples;
    r.min_assembly_events = cfg_.adaptive_min_assembly_events;
    r.event_residence_ps = cfg_.adaptive_event_residence_ps;
    r.min_time_over_its = cfg_.adaptive_min_time_over_its;
    r.smooth_window_frames = cfg_.adaptive_smooth_window_frames;
    r.min_residence_frames = cfg_.adaptive_min_residence_frames;
    r.pattern_downsample = cfg_.adaptive_pattern_downsample;
    return r;
}

std::filesystem::path AdaptiveSampler::chunk_base(int chunk) const {
    return cfg_.resultDir() / ("md_chunk_" + zero_padded(chunk));
}

std::filesystem::path AdaptiveSampler::chunk_tpr(int chunk) const { return chunk_base(chunk).string() + ".tpr"; }
std::filesystem::path AdaptiveSampler::chunk_colvar(int chunk) const { return cfg_.resultDir() / ("COLVAR_chunk_" + zero_padded(chunk) + ".dat"); }
std::filesystem::path AdaptiveSampler::chunk_plumed(int chunk) const { return cfg_.resultDir() / ("plumed_chunk_" + zero_padded(chunk) + ".dat"); }

bool AdaptiveSampler::chunk_completed(int chunk) const {
    const auto base = chunk_base(chunk).string();
    const auto log_path = base + ".log";
    if (!std::filesystem::exists(log_path)) return false;
    if (read_tail(log_path, 64 * 1024).find("Finished mdrun") == std::string::npos) return false;
    return std::filesystem::exists(base + ".gro") && std::filesystem::exists(base + ".cpt");
}

std::vector<int> AdaptiveSampler::completed_chunks() const {
    std::vector<int> out;
    for (int i = 0; i < cfg_.adaptiveMaxChunks(); ++i) {
        if (!chunk_completed(i)) break;
        out.push_back(i);
    }
    return out;
}

std::vector<std::filesystem::path> AdaptiveSampler::completed_colvars() const {
    std::vector<std::filesystem::path> paths;
    for (const auto chunk : completed_chunks()) {
        const auto path = chunk_colvar(chunk);
        if (std::filesystem::exists(path)) paths.push_back(path);
    }
    return paths;
}

void AdaptiveSampler::run() {
    const auto npt_gro = cfg_.resultDir() / "npt.gro";
    const auto npt_cpt = cfg_.resultDir() / "npt.cpt";
    const bool require_inputs = !sh_.dryRun();
    if (require_inputs && !std::filesystem::exists(npt_gro)) throw std::runtime_error("Adaptive production requires existing npt.gro: " + npt_gro.string());
    if (require_inputs && !std::filesystem::exists(npt_cpt)) throw std::runtime_error("Adaptive production requires existing npt.cpt: " + npt_cpt.string());

    const auto chunk_steps = cfg_.adaptiveChunkNsteps();
    const auto max_chunks = cfg_.adaptiveMaxChunks();
    if (chunk_steps <= 0 || max_chunks <= 0) throw std::runtime_error("Invalid adaptive chunk settings.");

    auto ledger = RunLedger::load(cfg_.resultDir() / "adaptive.ledger");
    ledger.require("chunk_nsteps", std::to_string(chunk_steps), "--adaptive-chunk-us/--md-dt-ps");
    ledger.require("dt_ps", to_string_fixed(cfg_.md_dt_ps, 6), "--md-dt-ps");
    ledger.require("n_prot", std::to_string(cfg_.n_prot), "--n-prot");
    ledger.require("plumed_stride", std::to_string(cfg_.plumed_stride), "--plumed-stride");
    if (!sh_.dryRun()) ledger.save(cfg_.resultDir() / "adaptive.ledger");

    const auto sampling_rules = rules();
    const double dt_colvar_ps = cfg_.md_dt_ps * static_cast<double>(cfg_.plumed_stride);
    const auto done = completed_chunks();
    const int start = done.empty() ? 0 : done.back() + 1;

    std::cerr << (start ? "Resuming" : "Starting") << " adaptive production from chunk " << start
              << " (t = " << static_cast<double>(start) * static_cast<double>(chunk_steps) * cfg_.md_dt_ps / 1000.0
              << " ns)\n";

    const auto index_path = cfg_.resultDir() / "index.ndx";
    generate_index(sh_, cfg_, npt_gro);

    for (int chunk = start; chunk < max_chunks; ++chunk) {
        std::cerr << "\n============================================================\n"
                  << "ADAPTIVE CHUNK " << chunk << " / " << max_chunks << '\n'
                  << "============================================================\n";

        const auto input_gro = chunk == 0 ? npt_gro : std::filesystem::path(chunk_base(chunk - 1).string() + ".gro");
        const auto input_cpt = chunk == 0 ? npt_cpt : std::filesystem::path(chunk_base(chunk - 1).string() + ".cpt");
        if (require_inputs && !std::filesystem::exists(input_gro)) throw std::runtime_error("Missing input GRO for adaptive chunk: " + input_gro.string());
        if (require_inputs && !std::filesystem::exists(input_cpt)) throw std::runtime_error("Missing input CPT for adaptive chunk: " + input_cpt.string());

        const double tinit_ps = static_cast<double>(chunk) * static_cast<double>(chunk_steps) * cfg_.md_dt_ps;
        const auto mdp = write_md_chunk_mdp(cfg_, chunk, chunk_steps, tinit_ps);
        write_plumed_dat(cfg_, chunk_plumed(chunk), chunk_colvar(chunk));
        gmx_.grompp(mdp, input_gro, chunk_tpr(chunk), input_cpt, {}, index_path);
        gmx_.mdrun(chunk_base(chunk), chunk_tpr(chunk), chunk_plumed(chunk));

        if (sh_.dryRun()) continue;

        auto colvars = completed_colvars();
        const auto current = chunk_colvar(chunk);
        if (std::find(colvars.begin(), colvars.end(), current) == colvars.end()) colvars.push_back(current);

        auto metrics = evaluate_sampling(colvars, cfg_.n_prot, dt_colvar_ps, sampling_rules);
        if (sampling_rules.min_time_over_its > 0.0)
            apply_cv_readiness(metrics, sampling_rules, cv_readiness(colvars));
        append_metrics_jsonl(cfg_.resultDir() / "adaptive_sampling_metrics.jsonl", metrics, chunk);
        std::cerr << "\nSampling metrics:\n" << metrics.to_json(chunk) << '\n';

        try {
            RunDiagnostics diag;
            diag.energy = compute_energy_diagnostics(gmx_, chunk_base(chunk).string() + ".edr",
                                                      cfg_.resultDir() / "diagnostics_scratch",
                                                      /*include_pressure=*/true, /*include_density=*/true);
            diag.structure = compute_structural_diagnostics({current}, cfg_.n_prot, dt_colvar_ps,
                                                             sampling_rules.pair_contact_threshold);
            diag.has_structure = true;
            append_diagnostics_jsonl(cfg_.resultDir() / "diagnostics.jsonl", diag, "chunk_" + zero_padded(chunk));
        } catch (const std::exception& e) {
            std::cerr << "Warning: chunk " << chunk << " diagnostics failed (ignored, purely observational): "
                      << e.what() << '\n';
        }

        if (metrics.stop) {
            std::cerr << "\nSTOP CONDITION REACHED. Total sampled time: " << metrics.total_time_us << " us\n";
            report();
            return;
        }
    }

    std::cerr << "\nAdaptive maximum total time reached without satisfying the stop condition.\n";
    report();
}

CvReadiness AdaptiveSampler::cv_readiness(const std::vector<std::filesystem::path>& colvars) const {
    CvReadiness cv;
    const auto script = resolve_sibling_script(cfg_, "cv_readiness.py");
    if (script.empty()) {
        cv.note = "scripts/cv_readiness.py not found";
        std::cerr << "Warning: " << cv.note << " - the CV-trainability gate cannot pass\n";
        return cv;
    }
    const auto out = cfg_.resultDir() / "cv_readiness.txt";
    std::vector<std::string> cmd = {cfg_.python, script.string(), "--n-prot", std::to_string(cfg_.n_prot),
                                    "--out", out.string(), "--colvar"};
    for (const auto& path : colvars) cmd.push_back(path.string());
    const int rc = sh_.run(cmd, {}, /*check=*/false);
    cv = read_cv_readiness(out);
    if (rc != 0 || !cv.ok)
        std::cerr << "Warning: cv_readiness.py " << (rc ? "exited " + std::to_string(rc) : "gave no usable result")
                  << (cv.note.empty() ? "" : " (" + cv.note + ")") << " - the CV-trainability gate cannot pass\n";
    return cv;
}

} // namespace cg
