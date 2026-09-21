#pragma once

#include <cstdint>
#include <filesystem>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cg {

struct Config;

using Lines = std::vector<std::pair<std::string, std::string>>;

void emit(std::ostringstream& out, const Lines& lines);
void blank(std::ostringstream& out);
std::string fp(double x);

struct MdSpec {
    std::string title;
    std::filesystem::path filename;
    double dt_ps = 0.0;
    std::int64_t nsteps = 0;
    bool continuation = true;
    bool pressure = true;
    bool gen_vel = false;
    double tau_p = 12.0;
    double tinit_ps = 0.0;
    /// Distinguishes consecutive segments of one continued run.
    int segment_index = 0;
    /// Overrides Config::nst_xtc when > 0.
    int nst_xtc = 0;
};

constexpr int kSeedNvt = 0;
constexpr int kSeedNpt = 1;
constexpr int kSeedProduction = 2;
constexpr int kSeedAdaptiveBase = 100;
constexpr int kSeedMetadBase = 100000;
/// Walker w of a metadynamics batch adds w times this to the batch's own
/// segment index, so walkers started from one checkpoint with a fixed --seed
/// still draw different thermostat noise.
constexpr int kSeedWalkerStride = 10000;
constexpr int kSeedRelax = 3;

// C-rescale requires tau-p >= 25 * nstpcouple * dt, and grompp makes falling
// short a WARNING - which, with no -maxwarn, is fatal. At kNstCouple = 20 and
// dt = 10 fs that floor is 5.0 ps, so the previous 4.0 failed every NPT grompp
// the moment nstpcouple started being written out explicitly: before that
// GROMACS picked its own default of 10, the floor was 2.5, and 4.0 slipped
// under it. Raising tau-p rather than lowering nstpcouple or passing -maxwarn:
// nstpcouple = nstlist is what makes the coupling cheap, and the warning is a
// statement about integration quality, not noise to silence.
//
// 12.0 ps is the value the MdpSpec default already carried and a standard
// Martini choice; at kNstCouple = 20 it clears the floor by a factor of 2.4,
// so the criterion survives a change of nstlist or dt without silently
// breaking again.
constexpr double kTauP = 12.0;

/// The full .mdp text of one dynamics segment.
std::string production_mdp(const Config& cfg, const MdSpec& spec);

/// Writes production_mdp() to cfg.systemDir() / spec.filename (not in
/// --dry-run).
void write_md_spec(const Config& cfg, const MdSpec& spec);

/// ld-seed of one segment: Config::seed plus the segment's own offset.
int segment_ld_seed(const Config& cfg, int segment_index);

} // namespace cg
