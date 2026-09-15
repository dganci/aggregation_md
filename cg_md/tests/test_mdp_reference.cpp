#include "test_framework.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "MdpSpec.hpp"
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

CG_TEST(each_segment_gets_its_own_ld_seed_when_a_seed_is_fixed) {
    Config cfg = make_config();
    cfg.seed = 42;

    write_nvt_mdp(cfg);
    write_npt_mdp(cfg);
    write_md_chunk_mdp(cfg, 0, 1000, 0.0);
    write_md_chunk_mdp(cfg, 1, 1000, 10.0);
    write_metad_batch_mdp(cfg, 0, 1000, 0.0);
    write_metad_batch_mdp(cfg, 1, 1000, 10.0);

    const auto seed_of = [&](const std::string& file) {
        const auto text = read_text(cfg.systemDir() / file);
        const auto pos = text.find("ld-seed = ");
        CG_CHECK(pos != std::string::npos);
        return text.substr(pos + 10, text.find('\n', pos) - pos - 10);
    };

    const std::vector<std::string> seeds = {
        seed_of("nvt.mdp"), seed_of("npt.mdp"),
        seed_of("md_chunk_000.mdp"), seed_of("md_chunk_001.mdp"),
        seed_of("md_metad_000.mdp"), seed_of("md_metad_001.mdp")};

    for (std::size_t i = 0; i < seeds.size(); ++i)
        for (std::size_t j = i + 1; j < seeds.size(); ++j)
            CG_CHECK(seeds[i] != seeds[j]);

    write_nvt_mdp(cfg);
    CG_CHECK(seed_of("nvt.mdp") == seeds[0]);

    write_metad_batch_mdp(cfg, 1, 1000, 10.0, 1);
    write_metad_batch_mdp(cfg, 1, 1000, 10.0, 2);
    const auto w0 = seed_of("md_metad_001.mdp");
    const auto w1 = seed_of("md_metad_001_w1.mdp");
    const auto w2 = seed_of("md_metad_001_w2.mdp");
    CG_CHECK(w0 != w1 && w1 != w2 && w0 != w2);
    CG_CHECK(std::stoi(w1) - std::stoi(w0) == kSeedWalkerStride);
}

CG_TEST(an_unset_seed_stays_at_the_gromacs_random_convention) {
    Config cfg = make_config();
    write_md_chunk_mdp(cfg, 3, 1000, 30.0);
    CG_CHECK(contains(read_text(cfg.systemDir() / "md_chunk_003.mdp"), "ld-seed = -1"));
}

CG_TEST(metad_falls_back_to_md_total_us_and_a_single_batch_by_default) {
    Config cfg = make_config();
    cfg.md_total_us = 2.0;
    CG_CHECK_EQ(cfg.metadNsteps(), cfg.mdNsteps());
    CG_CHECK_EQ(cfg.metadNumChunks(), 1);
    CG_CHECK_EQ(cfg.metadChunkNsteps(), cfg.metadNsteps());
}

CG_TEST(every_dynamics_mdp_sets_the_trajectory_precision) {
    auto cfg = make_config();
    write_nvt_mdp(cfg);
    write_npt_mdp(cfg);
    write_md_mdp(cfg);

    for (const char* name : {"nvt.mdp", "npt.mdp", "md.mdp"})
        CG_CHECK(contains(read_text(cfg.systemDir() / name), "compressed-x-precision = 100"));

    cfg.xtc_precision = 1000;
    write_md_mdp(cfg);
    CG_CHECK(contains(read_text(cfg.systemDir() / "md.mdp"), "compressed-x-precision = 1000"));

    std::filesystem::remove_all(cfg.project_dir);
}

CG_TEST(dynamics_mdps_match_the_martini_reference) {
    auto cfg = make_config();
    cfg.temperature_K = 310.0;
    write_npt_mdp(cfg);
    write_md_mdp(cfg);

    for (const char* name : {"npt.mdp", "md.mdp"}) {
        const auto text = read_text(cfg.systemDir() / name);
        CG_CHECK(contains(text, "nsttcouple = 20"));
        CG_CHECK(contains(text, "nstpcouple = 20"));
        CG_CHECK(contains(text, "tau-p = 4.000000"));
        CG_CHECK(contains(text, "refcoord-scaling = all"));
        CG_CHECK(contains(text, "nstcomm = 100"));
        CG_CHECK(contains(text, "verlet-buffer-tolerance = -1.0000"));
        CG_CHECK(contains(text, "rlist = 1.35"));
        CG_CHECK(contains(text, "compressed-x-precision = 100"));
        CG_CHECK(!contains(text, "energygrps"));
    }

    write_nvt_mdp(cfg);
    const auto nvt = read_text(cfg.systemDir() / "nvt.mdp");
    CG_CHECK(contains(nvt, "pcoupl = no"));
    CG_CHECK(contains(nvt, "nsttcouple = 20"));
    CG_CHECK(!contains(nvt, "nstpcouple"));

    write_em_mdp(cfg);
    const auto em = read_text(cfg.systemDir() / "emin.mdp");
    CG_CHECK(!contains(em, "nsttcouple"));
    CG_CHECK(!contains(em, "tcoupl"));
    CG_CHECK(contains(em, "lincs-warnangle = 90"));
    CG_CHECK(!contains(read_text(cfg.systemDir() / "md.mdp"), "lincs-warnangle"));

    std::filesystem::remove_all(cfg.project_dir);
}

CG_TEST(the_aligned_thermostat_keeps_ions_with_the_solvent) {
    auto cfg = make_config();
    cfg.thermostat_mode = "protein-solvent";
    cfg.temperature_K = 310.0;
    write_md_mdp(cfg);

    const auto text = read_text(cfg.systemDir() / "md.mdp");
    CG_CHECK(contains(text, "tc-grps = Protein Solvent_and_ions"));
    CG_CHECK(contains(text, "tau-t = 1.0 1.0"));
    CG_CHECK(!contains(text, "ION"));

    std::filesystem::remove_all(cfg.project_dir);
}
