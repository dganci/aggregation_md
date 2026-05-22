#pragma once

#include <string>

namespace cg {

class Shell;
struct Config;

std::string run_mkdssp(const Shell& sh, const Config& cfg);
std::string dssp_to_martinize_ss(const std::string& dssp_output, int seq_length);

} // namespace cg
