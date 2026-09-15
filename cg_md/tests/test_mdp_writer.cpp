#include "test_framework.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "MdpWriter.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

using namespace cg;

namespace {

Config make_config() {
    Config cfg;
    cfg.project_dir = std::filesystem::temp_directory_path() / "cg_md_test_mdp";
    std::filesystem::create_directories(cfg.systemDir());
    cfg.protomer_name = "test";
    cfg.md_dt_ps = 0.01;
    return cfg;
}

bool contains(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}

} // namespace

CG_TEST(every_mdp_carries_martinis_uniform_dielectric) {
    const auto cfg = make_config();
    write_em_mdp(cfg);
    write_nvt_mdp(cfg);
    write_npt_mdp(cfg);
    write_md_mdp(cfg);
    const auto chunk = write_md_chunk_mdp(cfg, 0, 1000, 0.0);
    const auto batch = write_metad_batch_mdp(cfg, 0, 1000, 0.0);

    for (const auto& file : {cfg.systemDir() / "emin.mdp", cfg.systemDir() / "nvt.mdp",
                             cfg.systemDir() / "npt.mdp", cfg.systemDir() / "md.mdp",
                             chunk, batch}) {
        const auto text = read_text(file);
        CG_CHECK(contains(text, "epsilon_r = 15.0"));
        CG_CHECK(contains(text, "epsilon_rf = 0"));
        CG_CHECK(contains(text, "coulombtype = Reaction-Field"));
        CG_CHECK(contains(text, "rcoulomb = 1.10"));
        CG_CHECK(contains(text, "rvdw = 1.10"));
    }
}

CG_TEST(chunk_mdps_are_numbered_and_advance_the_simulation_clock) {
    const auto cfg = make_config();

    const auto first = write_md_chunk_mdp(cfg, 0, 50000, 0.0);
    const auto second = write_md_chunk_mdp(cfg, 1, 50000, 500.0);

    CG_CHECK(first != second);
    CG_CHECK(contains(first.filename().string(), "md_chunk_000"));
    CG_CHECK(contains(second.filename().string(), "md_chunk_001"));

    const auto text_first = read_text(first);
    const auto text_second = read_text(second);

    CG_CHECK(contains(text_first, "tinit = 0.000000"));
    CG_CHECK(contains(text_second, "tinit = 500.000000"));
    CG_CHECK(contains(text_second, "nsteps = 50000"));
    CG_CHECK(contains(text_second, "continuation = yes"));
    CG_CHECK(contains(text_second, "gen_vel = no"));

    CG_CHECK(std::filesystem::exists(first));
    CG_CHECK(std::filesystem::exists(second));
}

CG_TEST(metad_batch_mdps_are_numbered_separately_from_adaptive_chunks) {
    const auto cfg = make_config();
    const auto batch = write_metad_batch_mdp(cfg, 2, 25000, 500.0);

    CG_CHECK(contains(batch.filename().string(), "md_metad_002"));
    const auto text = read_text(batch);
    CG_CHECK(contains(text, "tinit = 500.000000"));
    CG_CHECK(contains(text, "nsteps = 25000"));
}

CG_TEST(metad_batch_arithmetic_covers_the_requested_total_exactly_once) {
    Config cfg = make_config();
    cfg.md_dt_ps = 0.01;
    cfg.metad_total_us = 0.001;
    cfg.metad_chunk_us = 0.0003;

    CG_CHECK_EQ(cfg.metadNsteps(), static_cast<std::int64_t>(100000));
    CG_CHECK_EQ(cfg.metadChunkNsteps(), static_cast<std::int64_t>(30000));
    CG_CHECK_EQ(cfg.metadNumChunks(), 4);
    CG_CHECK_EQ(cfg.effectiveMetadNsteps(), static_cast<std::int64_t>(120000));

    std::int64_t covered = 0;
    for (int b = 0; b < cfg.metadNumChunks(); ++b) covered += cfg.metadChunkNsteps();
    CG_CHECK_EQ(covered, cfg.effectiveMetadNsteps());
}

CG_TEST(extending_a_batched_run_keeps_the_clock_aligned_with_the_state) {
    Config first = make_config();
    first.md_dt_ps = 0.01;
    first.metad_total_us = 0.001;
    first.metad_chunk_us = 0.0003;

    Config extended = first;
    extended.metad_total_us = 0.0012;

    CG_CHECK_EQ(first.metadChunkNsteps(), extended.metadChunkNsteps());
    CG_CHECK_EQ(first.metadNumChunks(), 4);
    CG_CHECK_EQ(extended.metadNumChunks(), 4);
    CG_CHECK_EQ(first.effectiveMetadNsteps(), extended.effectiveMetadNsteps());
    for (int b = 0; b < extended.metadNumChunks(); ++b) {
        const auto t_first = static_cast<std::int64_t>(b) * first.metadChunkNsteps();
        const auto t_extended = static_cast<std::int64_t>(b) * extended.metadChunkNsteps();
        CG_CHECK_EQ(t_first, t_extended);
    }
}
