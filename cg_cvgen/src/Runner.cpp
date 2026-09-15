#include "Runner.hpp"
#include "FileUtils.hpp"
#include "TextUtils.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cgcv {
namespace {

int normalize_exit_code(int status) {
#ifdef _WIN32
    return status;
#else
    if (status == 0) return 0;
    return (status & 0xFF) == 0 ? status >> 8 : status;
#endif
}

} // namespace

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

    const int code = normalize_exit_code(std::system(cmd.c_str()));
    if (code != 0) throw std::runtime_error("CV backend failed with exit code " + std::to_string(code));
    return 0;
}

} // namespace cgcv
