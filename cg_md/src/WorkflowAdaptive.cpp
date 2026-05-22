#include "Workflow.hpp"
#include "FileUtils.hpp"
#include "IndexBuilder.hpp"
#include "MdpWriter.hpp"
#include "PlumedWriter.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace cg {

SamplingRules Workflow::adaptive_rules() const {
    SamplingRules r;
    r.min_total_us = cfg_.adaptive_min_total_us;
    r.min_cn_state_transitions = cfg_.adaptive_min_cn_state_transitions;
    r.min_bidirectional_events = cfg_.adaptive_min_bidirectional_events;
    r.min_cn_range = cfg_.adaptive_min_cn_range;
    r.min_rg_global_range = cfg_.adaptive_min_rg_global_range;
    r.min_unique_contact_patterns = cfg_.adaptive_min_unique_contact_patterns;
    r.min_largest_cluster_unique = cfg_.adaptive_min_largest_cluster_unique;
    r.pair_contact_threshold = cfg_.adaptive_pair_contact_threshold;
    r.smooth_window_frames = cfg_.adaptive_smooth_window_frames;
    r.min_residence_frames = cfg_.adaptive_min_residence_frames;
    r.pattern_downsample = cfg_.adaptive_pattern_downsample;
    return r;
}

std::filesystem::path Workflow::chunk_base(int chunk) const {
    return cfg_.resultDir() / ("md_chunk_" + zero_padded(chunk));
}

std::filesystem::path Workflow::chunk_tpr(int chunk) const { return chunk_base(chunk).string() + ".tpr"; }
std::filesystem::path Workflow::chunk_colvar(int chunk) const { return cfg_.resultDir() / ("COLVAR_chunk_" + zero_padded(chunk) + ".dat"); }
std::filesystem::path Workflow::chunk_plumed(int chunk) const { return cfg_.resultDir() / ("plumed_chunk_" + zero_padded(chunk) + ".dat"); }

bool Workflow::chunk_completed(int chunk) const {
    const auto log_path = chunk_base(chunk).string() + ".log";
    return std::filesystem::exists(log_path) && read_text(log_path).find("Finished mdrun") != std::string::npos;
}

std::vector<int> Workflow::completed_chunks() const {
    std::vector<int> out;
    for (int i = 0; i < cfg_.adaptiveMaxChunks(); ++i)
        if (chunk_completed(i)) out.push_back(i);
    return out;
}

std::vector<std::filesystem::path> Workflow::completed_colvars() const {
    std::vector<std::filesystem::path> paths;
    for (const auto chunk : completed_chunks()) {
        const auto path = chunk_colvar(chunk);
        if (std::filesystem::exists(path)) paths.push_back(path);
    }
    return paths;
}

void Workflow::run_adaptive_production() {
    const auto npt_gro = cfg_.resultDir() / "npt.gro";
    const auto npt_cpt = cfg_.resultDir() / "npt.cpt";
    if (!std::filesystem::exists(npt_gro)) throw std::runtime_error("Adaptive production requires existing npt.gro: " + npt_gro.string());
    if (!std::filesystem::exists(npt_cpt)) throw std::runtime_error("Adaptive production requires existing npt.cpt: " + npt_cpt.string());

    const auto chunk_steps = cfg_.adaptiveChunkNsteps();
    const auto max_chunks = cfg_.adaptiveMaxChunks();
    if (chunk_steps <= 0 || max_chunks <= 0) throw std::runtime_error("Invalid adaptive chunk settings.");

    write_md_chunk_mdp(cfg_, chunk_steps);
    const auto rules = adaptive_rules();
    const double dt_colvar_ps = cfg_.md_dt_ps * static_cast<double>(cfg_.plumed_stride);
    const auto done = completed_chunks();
    const int start = done.empty() ? 0 : done.back() + 1;

    std::cerr << (start ? "Resuming" : "Starting") << " adaptive production from chunk " << start << '\n';

    for (int chunk = start; chunk < max_chunks; ++chunk) {
        std::cerr << "\n============================================================\n"
                  << "ADAPTIVE CHUNK " << chunk << " / " << max_chunks << '\n'
                  << "============================================================\n";

        const auto input_gro = chunk == 0 ? npt_gro : std::filesystem::path(chunk_base(chunk - 1).string() + ".gro");
        const auto input_cpt = chunk == 0 ? npt_cpt : std::filesystem::path(chunk_base(chunk - 1).string() + ".cpt");
        if (!std::filesystem::exists(input_gro)) throw std::runtime_error("Missing input GRO for adaptive chunk: " + input_gro.string());
        if (!std::filesystem::exists(input_cpt)) throw std::runtime_error("Missing input CPT for adaptive chunk: " + input_cpt.string());

        generate_index(sh_, cfg_, input_gro);
        write_plumed_dat(cfg_, chunk_plumed(chunk), chunk_colvar(chunk));
        grompp(cfg_.systemDir() / "md_chunk.mdp", input_gro, chunk_tpr(chunk), input_cpt, {}, cfg_.resultDir() / "index.ndx");
        mdrun(chunk_base(chunk), chunk_tpr(chunk), chunk_plumed(chunk));

        auto colvars = completed_colvars();
        const auto current = chunk_colvar(chunk);
        if (std::find(colvars.begin(), colvars.end(), current) == colvars.end()) colvars.push_back(current);

        const auto metrics = evaluate_sampling(colvars, cfg_.n_prot, dt_colvar_ps, rules);
        append_metrics_jsonl(cfg_.resultDir() / "adaptive_sampling_metrics.jsonl", metrics, chunk);
        std::cerr << "\nSampling metrics:\n" << metrics.to_json(chunk) << '\n';

        if (metrics.stop) {
            std::cerr << "\nSTOP CONDITION REACHED. Total sampled time: " << metrics.total_time_us << " us\n";
            return;
        }
    }

    std::cerr << "\nAdaptive maximum total time reached without satisfying the stop condition.\n";
}

} // namespace cg
