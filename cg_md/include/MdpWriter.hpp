#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace cg {

struct Config;

/// Renders the `tcoupl`/`tc-grps`/`tau-t`/`ref-t` block selected by
/// Config::thermostat_mode ("system" | "protein-solvent" | "legacy").
std::string thermostat_block(const Config& cfg);

/// Renders the non-bonded block shared by EM and every dynamics stage: Verlet
/// cut-off scheme, reaction-field electrostatics, the cut-offs, and -
/// critically - Martini's uniform relative dielectric `epsilon_r`
/// (Config::epsilon_r, 15 by default).
std::string nonbonded_block(const Config& cfg);

/// Each of these writes one GROMACS .mdp file under Config::systemDir().
void write_em_mdp(const Config& cfg);
void write_nvt_mdp(const Config& cfg);
void write_npt_mdp(const Config& cfg);
void write_md_mdp(const Config& cfg);

/// Like write_md_mdp(), but for one adaptive-sampling chunk / metadynamics
/// batch of `nsteps` steps whose simulation clock starts at `tinit_ps`.
std::filesystem::path write_md_chunk_mdp(const Config& cfg, int chunk,
                                         std::int64_t chunk_nsteps, double tinit_ps);
/// One .mdp per walker: they share everything but the ld-seed.
std::filesystem::path write_metad_batch_mdp(const Config& cfg, int batch,
                                            std::int64_t batch_nsteps, double tinit_ps,
                                            int walker = 0);

/// Production .mdp for the single-protomer relaxation run (--stage relax).
std::filesystem::path write_relax_mdp(const Config& cfg, std::int64_t nsteps);

} // namespace cg
