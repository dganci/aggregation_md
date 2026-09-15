#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cgcv {

/// Every option accepted by the cg_cvgen CLI, plus JSON serialization
/// (to_json()) so the values can be handed off unchanged to the Python DeepTICA
/// backend (scripts/cvgen_backend.py, see Runner::run()).
struct Config {
    std::filesystem::path input_colvar;
    std::filesystem::path output_dir = "CVs";
    std::filesystem::path backend;
    std::string python = "python3";

    /// Lag times to score, in FRAMES.
    std::vector<int> lags = {50, 100, 250, 500};
    /// Scan the implied timescales of the loaded data and take the lags from
    /// the plateau, overriding `lags`.
    bool auto_lags = true;
    std::vector<int> hidden_layers = {32, 16};
    int selected_lag = 0;
    int n_cvs = 3;
    int max_epochs = 200;
    int patience = 50;
    int probe_size = 1024;
    int seed = 42;

    double equilibration_time_ps = 50.0;
    double split_ratio = 0.8;
    double margin = 0.5;

    /// Selects which COLVAR columns are fed to DeepTICA as input features.
    std::string feature_regex = R"(^d_[0-9]+_[0-9]+$|^cn_[0-9]+_[0-9]+$|^cn_total$|^rg[0-9]+$|^rg_com$)";
    /// Sort each permutable descriptor block before training, removing the
    /// chain-labelling degeneracy of identical protomers (Samantray et al.;
    /// TICAgg).
    bool permutation_invariant = true;
    bool save_embeddings = false;
    bool dry_run = false;
    bool verbose = false;

    /// Throws std::runtime_error for any internally-inconsistent or
    /// missing-file option.
    void validate() const;
    /// output_dir/"cvgen_config.json" — the file to_json() is written to and
    /// that cvgen_backend.py reads back via --config.
    std::filesystem::path config_path() const;
    /// Serializes every option the Python backend needs (everything except
    /// `backend` and `python`, which only matter for building the command line
    /// in Runner::run()).
    std::string to_json() const;
    std::string summary() const;
};

Config parse_args(int argc, char** argv);
void print_help(const char* exe);

} // namespace cgcv
