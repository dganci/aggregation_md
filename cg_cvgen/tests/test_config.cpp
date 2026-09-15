#include "test_framework.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"

#include <exception>
#include <filesystem>
#include <regex>
#include <vector>

using namespace cgcv;

namespace {

std::filesystem::path temp_file(const std::string& name) {
    const auto path = std::filesystem::temp_directory_path() / ("cg_cvgen_test_" + name);
    write_text(path, "placeholder\n");
    return path;
}

Config make_valid_config() {
    Config cfg;
    cfg.input_colvar = temp_file("colvar.dat");
    cfg.backend = temp_file("backend.py");
    return cfg;
}

} // namespace

CG_TEST(config_validate_accepts_sane_defaults) {
    const auto cfg = make_valid_config();
    cfg.validate();
    std::filesystem::remove(cfg.input_colvar);
    std::filesystem::remove(cfg.backend);
}

CG_TEST(config_validate_rejects_missing_input_colvar) {
    auto cfg = make_valid_config();
    std::filesystem::remove(cfg.input_colvar);
    bool threw = false;
    try { cfg.validate(); } catch (const std::exception&) { threw = true; }
    CG_CHECK(threw);
    std::filesystem::remove(cfg.backend);
}

CG_TEST(config_validate_rejects_out_of_range_split_ratio) {
    auto cfg = make_valid_config();
    cfg.split_ratio = 1.5;
    bool threw = false;
    try { cfg.validate(); } catch (const std::exception&) { threw = true; }
    CG_CHECK(threw);
    std::filesystem::remove(cfg.input_colvar);
    std::filesystem::remove(cfg.backend);
}

CG_TEST(config_validate_rejects_negative_margin) {
    auto cfg = make_valid_config();
    cfg.margin = -0.1;
    bool threw = false;
    try { cfg.validate(); } catch (const std::exception&) { threw = true; }
    CG_CHECK(threw);
    std::filesystem::remove(cfg.input_colvar);
    std::filesystem::remove(cfg.backend);
}

CG_TEST(config_to_json_serializes_core_fields) {
    auto cfg = make_valid_config();
    cfg.n_cvs = 4;
    cfg.lags = {5, 9, 13};
    cfg.feature_regex = "^d_[0-9]+_[0-9]+$";
    const auto json = cfg.to_json();
    CG_CHECK(json.find("\"n_cvs\": 4") != std::string::npos);
    CG_CHECK(json.find("[5, 9, 13]") != std::string::npos);
    CG_CHECK(json.find("d_[0-9]+_[0-9]+") != std::string::npos);
    std::filesystem::remove(cfg.input_colvar);
    std::filesystem::remove(cfg.backend);
}

CG_TEST(config_path_is_output_dir_slash_cvgen_config_json) {
    Config cfg;
    cfg.output_dir = "some/output/dir";
    CG_CHECK_EQ(cfg.config_path(), std::filesystem::path("some/output/dir/cvgen_config.json"));
}

CG_TEST(default_feature_regex_matches_actual_cg_md_colvar_columns) {
    const std::regex rx(Config{}.feature_regex);
    for (const auto& col : {"d_1_2", "d_3_10", "cn_1_2", "cn_total", "rg1", "rg10", "rg_com"}) {
        CG_CHECK(std::regex_search(std::string(col), rx));
    }
    for (const auto& col : {"time", "walker"}) {
        CG_CHECK(!std::regex_search(std::string(col), rx));
    }
}
