#pragma once

#include "Config.hpp"

namespace cgcv {

class Runner {
public:
    explicit Runner(Config cfg) : cfg_(std::move(cfg)) {}
    int run() const;

private:
    Config cfg_;
};

} // namespace cgcv
