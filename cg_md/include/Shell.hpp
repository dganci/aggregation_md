#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

struct CommandResult {
    int exit_code = 0;
    std::string output;
};

class Shell {
public:
    explicit Shell(bool dry_run = false, bool verbose = true);

    void setDryRun(bool value);
    bool dryRun() const;

    int run(const std::vector<std::string>& args,
            const std::filesystem::path& cwd = {},
            bool check = true) const;

    int runShell(const std::string& command,
                 const std::filesystem::path& cwd = {},
                 bool check = true) const;

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
