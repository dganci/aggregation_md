#pragma once

#include <string>

namespace cg {

/// Number of walkers PLUMED reports in its own log (PLUMED.OUT), or -1 if the
/// log does not say.
int plumed_walker_count(const std::string& log_text);

/// True if `filename` is a PLUMED backup of `target`, i.e. bck.<N>.<target> or
/// bck.last.<target>.
bool is_plumed_backup(const std::string& filename, const std::string& target);

} // namespace cg
