#include "Shell.hpp"
#include "StringUtils.hpp"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace cg {
namespace {

FILE* open_pipe(const std::string& command) {
#ifdef _WIN32
    return _popen(command.c_str(), "r");
#else
    return popen(command.c_str(), "r");
#endif
}

int close_pipe(FILE* pipe) {
#ifdef _WIN32
    return _pclose(pipe);
#else
    return pclose(pipe);
#endif
}

} // namespace

Shell::Shell(bool dry_run, bool verbose) : dry_run_(dry_run), verbose_(verbose) {}

void Shell::setDryRun(bool value) { dry_run_ = value; }
bool Shell::dryRun() const { return dry_run_; }

int Shell::normalizeExit(int status) {
#ifdef _WIN32
    return status;
#else
    if (status == 0) return 0;
    return (status & 0xFF) == 0 ? status >> 8 : status;
#endif
}

std::string Shell::commandString(const std::vector<std::string>& args) {
    std::vector<std::string> quoted;
    quoted.reserve(args.size());
    for (const auto& arg : args) quoted.push_back(shell_quote(arg));
    return join(quoted, " ");
}

std::string Shell::withCwd(const std::string& command, const std::filesystem::path& cwd) {
    if (cwd.empty()) return command;
#ifdef _WIN32
    return "cd /d " + shell_quote(cwd.string()) + " && " + command;
#else
    return "cd " + shell_quote(cwd.string()) + " && " + command;
#endif
}

int Shell::execute(const std::string& command, const std::filesystem::path& cwd, bool check) const {
    const auto cmd = withCwd(command, cwd);
    if (verbose_) std::cerr << "\n$ " << cmd << '\n';
    if (dry_run_) return 0;

    const int raw = std::system(cmd.c_str());
    const int rc = normalizeExit(raw);
    if (check && rc != 0) {
        throw std::runtime_error("Command failed with code " + std::to_string(rc) +
                                 " (raw status " + std::to_string(raw) + "): " + cmd);
    }
    return rc;
}

int Shell::run(const std::vector<std::string>& args, const std::filesystem::path& cwd, bool check) const {
    return execute(commandString(args), cwd, check);
}

int Shell::runShell(const std::string& command, const std::filesystem::path& cwd, bool check) const {
    return execute(command, cwd, check);
}

CommandResult Shell::runCapture(const std::vector<std::string>& args,
                                const std::filesystem::path& cwd,
                                bool check) const {
    const auto cmd = withCwd(commandString(args) + " 2>&1", cwd);
    if (verbose_) std::cerr << "\n$ " << cmd << '\n';
    if (dry_run_) return {0, {}};

    auto* pipe = open_pipe(cmd);
    if (!pipe) throw std::runtime_error("Cannot open pipe for: " + cmd);

    std::array<char, 4096> buffer{};
    std::string output;
    while (std::fgets(buffer.data(), static_cast<int>(buffer.size()), pipe)) output += buffer.data();

    const int raw = close_pipe(pipe);
    const int rc = normalizeExit(raw);
    if (check && rc != 0) {
        throw std::runtime_error("Command failed with code " + std::to_string(rc) +
                                 " (raw status " + std::to_string(raw) + "): " + cmd +
                                 '\n' + output);
    }
    return {rc, output};
}

} // namespace cg
