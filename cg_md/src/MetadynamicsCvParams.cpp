#include "MetadynamicsRunner.hpp"
#include "Config.hpp"
#include "Shell.hpp"
#include "StringUtils.hpp"

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cg {
namespace {

constexpr const char* kCvParamsProbe =
    "import pickle, sys\n"
    "p = pickle.load(open(sys.argv[1], 'rb'))\n"
    "def as_list(v):\n"
    "    if v is None: return None\n"
    "    v = v.tolist() if hasattr(v, 'tolist') else list(v)\n"
    "    return [x if isinstance(x, str) else format(float(x), '.12g') for x in v]\n"
    "def emit(key, value):\n"
    "    if value is None: return\n"
    "    print(key + '\\t' + (','.join(value) if isinstance(value, list) else str(value)))\n"
    "emit('sigma', as_list(p.get('sigma')))\n"
    "emit('grid_min', as_list(p.get('grid_min')))\n"
    "emit('grid_max', as_list(p.get('grid_max')))\n"
    "fc = p.get('feature_cols')\n"
    "emit('feature_cols', list(fc) if fc is not None else None)\n"
    "emit('n_cvs', p.get('n_cvs'))\n"
    "pi = p.get('permutation_invariant')\n"
    "emit('permutation_invariant', None if pi is None else ('1' if pi else '0'))\n";


std::map<std::string, std::string> parse_key_value_lines(const std::string& text) {
    std::map<std::string, std::string> out;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        const auto tab = line.find('\t');
        if (tab == std::string::npos) continue;
        out[trim(line.substr(0, tab))] = trim(line.substr(tab + 1));
    }
    return out;
}


void require_cv_params(const std::filesystem::path& path) {
    if (std::filesystem::exists(path)) return;
    throw std::runtime_error(
        "Missing CV parameter pickle (pass --metad-sigma / --metad-grid-min / --metad-grid-max / "
        "--metad-feature-cols to bypass it): " + path.string());
}

} // namespace

void MetadynamicsRunner::load_cv_parameters() {
    if (cfg_.pmf) return;

    const bool need_sigma = cfg_.metad_sigma.empty();
    const bool need_grid = cfg_.metad_grid && (cfg_.metad_grid_min.empty() || cfg_.metad_grid_max.empty());
    const bool need_features = cfg_.metad_feature_cols.empty();

    if (sh_.dryRun()) {
        if (need_sigma) cfg_.metad_sigma = "SIGMA_FROM_CV_PARAMS";
        if (need_grid) {
            cfg_.metad_grid_min = "GRID_MIN_FROM_CV_PARAMS";
            cfg_.metad_grid_max = "GRID_MAX_FROM_CV_PARAMS";
        }
        return;
    }

    if (!need_sigma && !need_grid && !need_features) return;

    const auto params = cfg_.metadCvParamsPath();
    require_cv_params(params);

    const auto probe = parse_key_value_lines(
        sh_.runCapture({cfg_.python, "-c", kCvParamsProbe, params.string()}).output);

    const auto take = [&](const char* key) -> std::string {
        const auto it = probe.find(key);
        return it == probe.end() ? std::string{} : it->second;
    };

    if (need_sigma) cfg_.metad_sigma = take("sigma");
    if (need_grid) {
        cfg_.metad_grid_min = take("grid_min");
        cfg_.metad_grid_max = take("grid_max");
    }
    if (need_features) cfg_.metad_feature_cols = take("feature_cols");

    if (const auto pi = take("permutation_invariant"); !pi.empty()) {
        const bool trained_invariant = pi == "1";
        if (trained_invariant != cfg_.metad_permutation_invariant) {
            std::cerr << "Note: the trained CV is "
                      << (trained_invariant ? "permutation-invariant (sorted features)"
                                            : "index-ordered (raw features)")
                      << "; matching it.\n";
            cfg_.metad_permutation_invariant = trained_invariant;
        }
    }

    if (cfg_.metad_sigma.empty())
        throw std::runtime_error("Cannot extract sigma from " + params.string() + "; pass --metad-sigma");

    const auto sigma_values = split_csv(cfg_.metad_sigma);
    if (sigma_values.size() != static_cast<std::size_t>(cfg_.metad_nodes))
        throw std::runtime_error("--metad-sigma has " + std::to_string(sigma_values.size()) +
                                 " value(s) but --metad-nodes is " + std::to_string(cfg_.metad_nodes) +
                                 " - the CV was trained with a different number of components.");

    const auto n_cvs = take("n_cvs");
    if (!n_cvs.empty() && n_cvs != std::to_string(cfg_.metad_nodes))
        throw std::runtime_error("Trained CV has n_cvs=" + n_cvs + " but --metad-nodes is " +
                                 std::to_string(cfg_.metad_nodes));

    if (cfg_.metad_grid) {
        if (cfg_.metad_grid_min.empty() || cfg_.metad_grid_max.empty())
            throw std::runtime_error("Grid bounds not found in " + params.string() +
                                     "; pass --metad-grid-min/--metad-grid-max, or --no-metad-grid "
                                     "(note a gridless run also disables c(t) reweighting).");
        if (split_csv(cfg_.metad_grid_min).size() != static_cast<std::size_t>(cfg_.metad_nodes) ||
            split_csv(cfg_.metad_grid_max).size() != static_cast<std::size_t>(cfg_.metad_nodes))
            throw std::runtime_error("--metad-grid-min/--metad-grid-max must each hold exactly "
                                     "--metad-nodes values");
    }
}

} // namespace cg
