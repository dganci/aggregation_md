#include "Config.hpp"
#include "FileUtils.hpp"
#include "TextUtils.hpp"

#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace cgcv {

namespace {

std::string take(int& i, int argc, char** argv, const std::string& opt) {
    if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
    return argv[++i];
}

double to_double(const std::string& s, const std::string& opt) {
    try { return std::stod(s); }
    catch (...) { throw std::runtime_error("Invalid value for " + opt + ": " + s); }
}

int to_int(const std::string& s, const std::string& opt) {
    try { return std::stoi(s); }
    catch (...) { throw std::runtime_error("Invalid value for " + opt + ": " + s); }
}

} // namespace

void Config::validate() const {
    require_file(input_colvar, "input COLVAR file");
    require_file(backend, "Python backend");
    if (output_dir.empty()) throw std::runtime_error("--output-dir cannot be empty");
    if (python.empty()) throw std::runtime_error("--python cannot be empty");
    if (lags.empty()) throw std::runtime_error("--lags cannot be empty");
    if (hidden_layers.empty()) throw std::runtime_error("--hidden-layers cannot be empty");
    if (n_cvs <= 0) throw std::runtime_error("--n-cvs must be > 0");
    if (max_epochs <= 0) throw std::runtime_error("--max-epochs must be > 0");
    if (patience <= 0) throw std::runtime_error("--patience must be > 0");
    if (probe_size <= 0) throw std::runtime_error("--probe-size must be > 0");
    if (split_ratio <= 0.0 || split_ratio > 1.0) throw std::runtime_error("--split-ratio must be in (0, 1]");
    if (margin < 0.0) throw std::runtime_error("--margin must be >= 0");
    if (feature_regex.empty()) throw std::runtime_error("--feature-regex cannot be empty");
}

std::filesystem::path Config::config_path() const { return output_dir / "cvgen_config.json"; }

std::string Config::to_json() const {
    std::ostringstream out;
    out << "{\n"
        << "  \"input_colvar\": \"" << json_escape(std::filesystem::absolute(input_colvar).string()) << "\",\n"
        << "  \"output_dir\": \"" << json_escape(std::filesystem::absolute(output_dir).string()) << "\",\n"
        << "  \"lags\": " << to_json_array(lags) << ",\n"
        << "  \"selected_lag\": " << selected_lag << ",\n"
        << "  \"n_cvs\": " << n_cvs << ",\n"
        << "  \"hidden_layers\": " << to_json_array(hidden_layers) << ",\n"
        << "  \"max_epochs\": " << max_epochs << ",\n"
        << "  \"patience\": " << patience << ",\n"
        << "  \"probe_size\": " << probe_size << ",\n"
        << "  \"seed\": " << seed << ",\n"
        << "  \"equilibration_time_ps\": " << equilibration_time_ps << ",\n"
        << "  \"split_ratio\": " << split_ratio << ",\n"
        << "  \"margin\": " << margin << ",\n"
        << "  \"feature_regex\": \"" << json_escape(feature_regex) << "\",\n"
        << "  \"save_embeddings\": " << (save_embeddings ? "true" : "false") << ",\n"
        << "  \"verbose\": " << (verbose ? "true" : "false") << "\n"
        << "}\n";
    return out.str();
}

std::string Config::summary() const {
    std::ostringstream out;
    out << "cg_cvgen\n"
        << "  input_colvar      = " << input_colvar << '\n'
        << "  output_dir        = " << output_dir << '\n'
        << "  backend           = " << backend << '\n'
        << "  lags              = " << to_json_array(lags) << '\n'
        << "  selected_lag      = " << (selected_lag > 0 ? std::to_string(selected_lag) : "auto") << '\n'
        << "  n_cvs             = " << n_cvs << '\n'
        << "  hidden_layers     = " << to_json_array(hidden_layers) << '\n'
        << "  max_epochs        = " << max_epochs << '\n'
        << "  equilibration_ps  = " << equilibration_time_ps << '\n'
        << "  feature_regex     = " << feature_regex << '\n';
    return out.str();
}

void print_help(const char* exe) {
    std::cout
        << "Usage:\n  " << exe << " --input-colvar COLVAR --output-dir CVs [options]\n\n"
        << "Core options:\n"
        << "  --input-colvar PATH          Unbiased COLVAR file produced by PLUMED\n"
        << "  --output-dir PATH            Directory receiving CVs_torchscript.pt and parameters\n"
        << "  --lags 5,7,10,12             Lag times to score\n"
        << "  --selected-lag N             Force final lag; default chooses best score\n"
        << "  --n-cvs N                    Number of DeepTICA CVs\n"
        << "  --hidden-layers 32,16        Neural-network hidden layers\n"
        << "  --equilibration-time-ps X    Discard frames with time <= X\n"
        << "  --max-epochs N               Training epochs\n"
        << "  --feature-regex REGEX        Columns used as model input\n"
        << "  --python CMD                 Python executable\n"
        << "  --backend PATH               cvgen_backend.py path\n"
        << "  --save-embeddings            Store final embeddings\n"
        << "  --dry-run                    Write config and print backend command only\n";
}

Config parse_args(int argc, char** argv) {
    Config cfg;
    cfg.backend = default_backend_path(argc > 0 ? argv[0] : nullptr);

    const std::unordered_map<std::string, std::function<void(int&)>> opts = {
        {"--input-colvar", [&](int& i){ cfg.input_colvar = take(i, argc, argv, argv[i]); }},
        {"--output-dir", [&](int& i){ cfg.output_dir = take(i, argc, argv, argv[i]); }},
        {"--backend", [&](int& i){ cfg.backend = take(i, argc, argv, argv[i]); }},
        {"--python", [&](int& i){ cfg.python = take(i, argc, argv, argv[i]); }},
        {"--lags", [&](int& i){ cfg.lags = parse_int_list(take(i, argc, argv, argv[i])); }},
        {"--lag-list", [&](int& i){ cfg.lags = parse_int_list(take(i, argc, argv, argv[i])); }},
        {"--selected-lag", [&](int& i){ cfg.selected_lag = to_int(take(i, argc, argv, argv[i]), "--selected-lag"); }},
        {"--n-cvs", [&](int& i){ cfg.n_cvs = to_int(take(i, argc, argv, argv[i]), "--n-cvs"); }},
        {"--hidden-layers", [&](int& i){ cfg.hidden_layers = parse_layers(take(i, argc, argv, argv[i])); }},
        {"--max-epochs", [&](int& i){ cfg.max_epochs = to_int(take(i, argc, argv, argv[i]), "--max-epochs"); }},
        {"--patience", [&](int& i){ cfg.patience = to_int(take(i, argc, argv, argv[i]), "--patience"); }},
        {"--probe-size", [&](int& i){ cfg.probe_size = to_int(take(i, argc, argv, argv[i]), "--probe-size"); }},
        {"--seed", [&](int& i){ cfg.seed = to_int(take(i, argc, argv, argv[i]), "--seed"); }},
        {"--equilibration-time-ps", [&](int& i){ cfg.equilibration_time_ps = to_double(take(i, argc, argv, argv[i]), "--equilibration-time-ps"); }},
        {"--split-ratio", [&](int& i){ cfg.split_ratio = to_double(take(i, argc, argv, argv[i]), "--split-ratio"); }},
        {"--margin", [&](int& i){ cfg.margin = to_double(take(i, argc, argv, argv[i]), "--margin"); }},
        {"--feature-regex", [&](int& i){ cfg.feature_regex = take(i, argc, argv, argv[i]); }}
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { print_help(argv[0]); std::exit(0); }
        if (arg == "--save-embeddings") { cfg.save_embeddings = true; continue; }
        if (arg == "--dry-run") { cfg.dry_run = true; continue; }
        if (arg == "--verbose") { cfg.verbose = true; continue; }
        const auto it = opts.find(arg);
        if (it == opts.end()) throw std::runtime_error("Unknown option: " + arg);
        it->second(i);
    }

    if (cfg.input_colvar.empty()) throw std::runtime_error("Missing required option: --input-colvar");
    cfg.validate();
    return cfg;
}

} // namespace cgcv
