#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace cg {

struct Config;

std::string thermostat_block(const Config& cfg);
void write_em_mdp(const Config& cfg);
void write_nvt_mdp(const Config& cfg);
void write_npt_mdp(const Config& cfg);
void write_md_mdp(const Config& cfg);
void write_metad_mdp(const Config& cfg);
void write_md_chunk_mdp(const Config& cfg, std::int64_t chunk_nsteps);

} // namespace cg
