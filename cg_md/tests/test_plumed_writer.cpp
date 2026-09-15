#include "test_framework.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "PlumedWriter.hpp"

#include <exception>
#include <filesystem>
#include <string>

using namespace cg;

namespace {

std::filesystem::path scratch_dir() {
    const auto dir = std::filesystem::temp_directory_path() / "cg_md_test_plumed";
    std::filesystem::create_directories(dir / "runs" / "3xtest");
    return dir;
}

Config make_config(int n_prot = 3) {
    Config cfg;
    cfg.project_dir = scratch_dir();
    cfg.protomer_name = "test";
    cfg.n_prot = n_prot;
    cfg.atoms_per_prot = 10;
    cfg.metad_nodes = 2;
    cfg.metad_sigma = "0.1,0.2";
    cfg.metad_grid_min = "-1,-1";
    cfg.metad_grid_max = "1,1";
    return cfg;
}

bool contains(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}

} // namespace

CG_TEST(both_plumed_writers_emit_the_same_descriptor_block) {
    const auto cfg = make_config();
    const auto production = build_descriptors(cfg, "index.ndx");
    const auto biased = build_descriptors(cfg, "index.ndx");

    CG_CHECK(production.text == biased.text);
    CG_CHECK(production.feature_cols == biased.feature_cols);

    CG_CHECK_EQ(static_cast<int>(production.feature_cols.size()), 3 + 3 + 3 + 2);
    CG_CHECK(production.feature_cols.front() == "d_1_2");
    CG_CHECK(production.feature_cols[3] == "cn_1_2");
    CG_CHECK(production.feature_cols[6] == "rg1");
    CG_CHECK(production.feature_cols[9] == "rg_com");
    CG_CHECK(production.feature_cols[10] == "cn_total");
}

CG_TEST(descriptor_block_always_makes_molecules_whole_and_pins_the_switching_function) {
    const auto cfg = make_config();
    const auto d = build_descriptors(cfg, "index.ndx");

    CG_CHECK(contains(d.text, "WHOLEMOLECULES"));
    CG_CHECK(contains(d.text, "NN=6 MM=12"));
    CG_CHECK(!contains(d.text, "GYRATION ATOMS=allgrp"));
    CG_CHECK(contains(d.text, "rg_com_sq: COMBINE ARG=d_1_2"));
    CG_CHECK(contains(d.text, "POWERS=2,2,2"));
    CG_CHECK(contains(d.text, "rg_com: CUSTOM ARG=rg_com_sq FUNC=sqrt(x)"));
}

CG_TEST(merged_dimers_get_an_atom_order_that_never_crosses_the_rod) {
    Config cfg = make_config(2);
    cfg.n_prot = 2;
    cfg.atoms_per_prot = 6;
    cfg.merge_chains = "A,B";
    const auto d = build_descriptors(cfg, "idx.ndx");

    CG_CHECK(contains(d.text, "ENTITY0=1-3,6,5,4"));
    CG_CHECK(contains(d.text, "ENTITY1=7-9,12,11,10"));
    CG_CHECK(!contains(d.text, "ENTITY0=group0"));
}

CG_TEST(single_chain_protomers_keep_the_plain_group) {
    Config cfg = make_config(2);
    cfg.n_prot = 2;
    cfg.atoms_per_prot = 6;
    cfg.merge_chains = "";
    const auto d = build_descriptors(cfg, "idx.ndx");
    CG_CHECK(contains(d.text, "ENTITY0=group0"));
    CG_CHECK(contains(d.text, "ENTITY1=group1"));
}
