#pragma once

#include "Config.hpp"

namespace cgcv {

/// Writes cfg's JSON config (Config::to_json()) to Config::config_path(), then
/// invokes the Python backend (Config::backend) with `--config <that path>`.
class Runner {
public:
    explicit Runner(Config cfg) : cfg_(std::move(cfg)) {}
    /// Returns 0 on success; throws std::runtime_error if the backend exits
    /// non-zero.
    int run() const;

private:
    Config cfg_;
};

} // namespace cgcv
