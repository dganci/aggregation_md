#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

struct Config;
class Shell;

/// Thin wrappers around the `gmx` subcommands shared by every pipeline stage
/// (grompp/mdrun/energy/trjconv).
void append_pinning(const Config& cfg, std::vector<std::string>& cmd);

class GromacsDriver {
public:
    GromacsDriver(Config& cfg, Shell& sh);

    void grompp(const std::filesystem::path& mdp,
                const std::filesystem::path& coordinates,
                const std::filesystem::path& output,
                const std::filesystem::path& checkpoint = {},
                const std::filesystem::path& restraints = {},
                const std::filesystem::path& index = {}) const;

    void mdrun(const std::filesystem::path& deffnm,
               const std::filesystem::path& tpr = {},
               const std::filesystem::path& plumed = {},
               int ntomp_override = 0) const;

    /// Extracts `observables` (e.g. {"Temperature", "Pressure"}) from `edr`
    /// into a single `xvg` file, in one `gmx energy` invocation (each
    /// observable becomes its own column, in the order given).
    void energy(const std::filesystem::path& edr,
                const std::filesystem::path& xvg,
                const std::vector<std::string>& observables) const;

    /// Convenience overload for the common single-observable case.
    void energy(const std::filesystem::path& edr,
                const std::filesystem::path& xvg,
                const std::string& observable) const;

    void trjconv(const std::filesystem::path& tpr,
                 const std::filesystem::path& in,
                 const std::filesystem::path& out,
                 const std::string& groups,
                 const std::vector<std::string>& options) const;

private:
    Config& cfg_;
    Shell& sh_;
};

} // namespace cg
