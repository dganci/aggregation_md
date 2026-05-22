#include "Runner.hpp"
#include "FileUtils.hpp"
#include "TextUtils.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cgcv {

int Runner::run() const {
    ensure_dir(cfg_.output_dir);
    write_text(cfg_.config_path(), cfg_.to_json());

    const auto cmd = cfg_.python + " " + shell_quote(cfg_.backend.string()) +
                     " --config " + shell_quote(cfg_.config_path().string());

    std::cout << cfg_.summary() << "\n";
    if (cfg_.dry_run) {
        std::cout << cmd << '\n';
        return 0;
    }

    const int code = std::system(cmd.c_str());
    if (code != 0) throw std::runtime_error("CV backend failed with exit code " + std::to_string(code));
    return 0;
}

} // namespace cgcv
