#pragma once

#include <filesystem>
#include <string>

namespace cg {

struct Config;
class Shell;

/// Locates one of the pipeline's Python scripts by name: <project-dir>/scripts/
/// first, then the scripts/ directory next to the executable's source tree.
std::filesystem::path resolve_sibling_script(const Config& cfg, const char* name);

/// Locates scripts/analyze_run.py: Config::report_script when set, otherwise
/// the copy shipped next to the executable's source tree
/// (<project-dir>/scripts/analyze_run.py, then <exe dir>/../scripts/...).
std::filesystem::path resolve_report_script(const Config& cfg);

/// Runs scripts/contact_map.py over a centred trajectory, producing the
/// inter-chain residue-residue contact map and the NPMI shape descriptor.
void run_contact_map(const Config& cfg,
                     const Shell& sh,
                     const std::filesystem::path& traj,
                     const std::filesystem::path& top,
                     const std::filesystem::path& out_dir);

/// Runs the offline analysis/plotting backend over `run_dir` and writes its
/// figures and summary into `run_dir/report/`.
void run_report(const Config& cfg,
                const Shell& sh,
                const std::filesystem::path& run_dir,
                const std::string& mode);

} // namespace cg
