#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cgcv {

struct Config {
    std::filesystem::path input_colvar;
    std::filesystem::path output_dir = "CVs";
    std::filesystem::path backend;
    std::string python = "python3";

    std::vector<int> lags = {5, 7, 10, 12};
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

    std::string feature_regex = R"(^d[0-9]+$|^rg$|^cn_total$)";
    bool save_embeddings = false;
    bool dry_run = false;
    bool verbose = false;

    void validate() const;
    std::filesystem::path config_path() const;
    std::string to_json() const;
    std::string summary() const;
};

Config parse_args(int argc, char** argv);
void print_help(const char* exe);

} // namespace cgcv
