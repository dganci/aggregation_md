#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

struct CommandResult {
    int exit_code = 0;
    std::string output;
};

/// Runs external programs (gmx, martinize2, packmol, insane, plumed, python)
/// through the platform shell via std::system/popen.
class Shell {
public:
    explicit Shell(bool dry_run = false, bool verbose = true);

    void setDryRun(bool value);
    bool dryRun() const;

    /// Quotes and joins `args`, then runs them as one command.
    int run(const std::vector<std::string>& args,
            const std::filesystem::path& cwd = {},
            bool check = true) const;

    /// Like run(), but takes a raw shell command string (for pipes/redirects
    /// that cannot be expressed as a plain argument vector).
    int runShell(const std::string& command,
                 const std::filesystem::path& cwd = {},
                 bool check = true) const;

    /// Like run(), but captures combined stdout+stderr instead of letting it
    /// pass through to the terminal.
    CommandResult runCapture(const std::vector<std::string>& args,
                             const std::filesystem::path& cwd = {},
                             bool check = true) const;

private:
    bool dry_run_;
    bool verbose_;

    static int normalizeExit(int status);
    static std::string commandString(const std::vector<std::string>& args);
    static std::string withCwd(const std::string& command, const std::filesystem::path& cwd);
    int execute(const std::string& command, const std::filesystem::path& cwd, bool check) const;
};

} // namespace cg
