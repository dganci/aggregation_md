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

CG_TEST(metad_plumed_dat_uses_the_trained_feature_list_and_enables_reweighting) {
    auto cfg = make_config();
    const auto out = scratch_dir() / "plumed_metad.dat";

    cfg.metad_feature_cols = "d_1_2,d_1_3,d_2_3,cn_total,rg_com";

    MetadPlumedOptions options;
    options.restart = true;
    const auto files = write_metad_plumed_dat(cfg, out, "index.ndx", "CVs_torchscript.pt",
                                              "COLVAR", options);

    const auto text = read_text(out);
    CG_CHECK(contains(text, "ARG=d_1_2,d_1_3,d_2_3,cn_total,rg_com"));
    CG_CHECK_EQ(static_cast<int>(files.feature_cols.size()), 5);

    CG_CHECK(contains(text, "RESTART"));
    CG_CHECK(contains(text, "FILE=HILLS"));
    CG_CHECK(contains(text, "GRID_MIN=-1,-1"));
    CG_CHECK(contains(text, "GRID_BIN=200,200"));
    CG_CHECK(contains(text, "CALC_RCT"));
    CG_CHECK(contains(text, "metad.rbias"));
    CG_CHECK(contains(text, "COLVAR_monitor"));
    CG_CHECK(contains(text, "WHOLEMOLECULES"));

    std::filesystem::remove(out);
}

CG_TEST(metad_plumed_dat_refuses_a_feature_the_descriptors_do_not_define) {
    auto cfg = make_config();
    cfg.metad_feature_cols = "d_1_2,rg_of_something_else";

    bool threw = false;
    try {
        write_metad_plumed_dat(cfg, scratch_dir() / "bad.dat", "index.ndx",
                               "CVs_torchscript.pt", "COLVAR");
    } catch (const std::exception&) {
        threw = true;
    }
    CG_CHECK(threw);
}

CG_TEST(pmf_mode_biases_the_pair_distance_directly_and_needs_no_model) {
    auto cfg = make_config(/*n_prot=*/2);
    cfg.pmf = true;
    cfg.metad_nodes = 1;
    cfg.metad_sigma = "0.05";
    cfg.metad_grid_min = "0";
    cfg.metad_grid_max = "7";
    cfg.pmf_wall_nm = 6.0;
    cfg.pmf_wall_kappa = 2000.0;
    cfg.metad_feature_cols = "d_1_2";

    const auto out = scratch_dir() / "plumed_pmf.dat";
    write_metad_plumed_dat(cfg, out, "index.ndx", "unused.pt", "COLVAR");
    const auto text = read_text(out);

    CG_CHECK(contains(text, "UPPER_WALLS ARG=d_1_2"));
    CG_CHECK(contains(text, "AT=6.0000"));
    CG_CHECK(contains(text, "KAPPA=2000.0"));
    CG_CHECK(contains(text, "\nARG=d_1_2\n"));
    CG_CHECK(!contains(text, "PYTORCH_MODEL"));
    CG_CHECK(!contains(text, "mycv"));
    CG_CHECK(contains(text, "CALC_RCT"));
    CG_CHECK(contains(text, "wall.bias"));
    CG_CHECK(contains(text, "WHOLEMOLECULES"));

    std::filesystem::remove(out);
}

CG_TEST(metad_plumed_dat_without_restart_omits_the_restart_directive) {
    auto cfg = make_config();
    const auto out = scratch_dir() / "plumed_first.dat";
    write_metad_plumed_dat(cfg, out, "index.ndx", "CVs_torchscript.pt", "COLVAR");

    const auto text = read_text(out);
    CG_CHECK(!contains(text, "\nRESTART\n"));
    std::filesystem::remove(out);
}


CG_TEST(metadynamics_walls_the_grid_boundary) {
    auto cfg = make_config(/*n_prot=*/2);
    cfg.metadynamics = true;
    cfg.metad_nodes = 2;
    cfg.metad_sigma = "0.1,0.1";
    cfg.metad_grid_min = "-2,0";
    cfg.metad_grid_max = "2,10";
    cfg.metad_feature_cols = "d_1_2,cn_1_2";
    cfg.temperature_K = 310.0;
    cfg.metad_wall_kt = 50.0;
    cfg.metad_wall_margin_frac = 0.05;

    const auto out = scratch_dir() / "plumed_walls.dat";
    write_metad_plumed_dat(cfg, out, "index.ndx", "CVs_torchscript.pt", "COLVAR");
    const auto text = read_text(out);

    CG_CHECK(contains(text, "lwall: LOWER_WALLS ARG=mycv.node-0,mycv.node-1 AT=-1.8000,0.5000"));
    CG_CHECK(contains(text, "uwall: UPPER_WALLS ARG=mycv.node-0,mycv.node-1 AT=1.8000,9.5000"));
    CG_CHECK(contains(text, "KAPPA=3221.8541,515.4967"));
    CG_CHECK(contains(text, "GRID_MIN=-2,0"));
    CG_CHECK(contains(text, "GRID_MAX=2,10"));
    CG_CHECK(contains(text, "lwall.bias"));
    CG_CHECK(contains(text, "uwall.bias"));

    std::filesystem::remove(out);
}

CG_TEST(metadynamics_walls_can_be_disabled) {
    auto cfg = make_config(/*n_prot=*/2);
    cfg.metadynamics = true;
    cfg.metad_nodes = 1;
    cfg.metad_sigma = "0.1";
    cfg.metad_grid_min = "-2";
    cfg.metad_grid_max = "2";
    cfg.metad_feature_cols = "d_1_2";
    cfg.metad_wall_kt = 0.0;

    const auto out = scratch_dir() / "plumed_nowalls.dat";
    write_metad_plumed_dat(cfg, out, "index.ndx", "CVs_torchscript.pt", "COLVAR");
    const auto text = read_text(out);

    CG_CHECK(!contains(text, "LOWER_WALLS"));
    CG_CHECK(!contains(text, "UPPER_WALLS"));
    CG_CHECK(!contains(text, "lwall"));
    CG_CHECK(!contains(text, "uwall"));

    std::filesystem::remove(out);
}

CG_TEST(non_numeric_grid_bounds_skip_the_walls_instead_of_throwing) {
    auto cfg = make_config(/*n_prot=*/2);
    cfg.metadynamics = true;
    cfg.metad_nodes = 1;
    cfg.metad_sigma = "SIGMA_FROM_CV_PARAMS";
    cfg.metad_grid_min = "GRID_MIN_FROM_CV_PARAMS";
    cfg.metad_grid_max = "GRID_MAX_FROM_CV_PARAMS";
    cfg.metad_feature_cols = "d_1_2";

    const auto out = scratch_dir() / "plumed_placeholder.dat";
    write_metad_plumed_dat(cfg, out, "index.ndx", "CVs_torchscript.pt", "COLVAR");
    const auto text = read_text(out);

    CG_CHECK(!contains(text, "LOWER_WALLS"));
    CG_CHECK(!contains(text, "lwall.bias"));
    CG_CHECK(contains(text, "METAD"));

    std::filesystem::remove(out);
}
