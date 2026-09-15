#include "MdpSpec.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "MdpWriter.hpp"
#include "StringUtils.hpp"

#include <sstream>
#include <string>

namespace cg {

void emit(std::ostringstream& out, const Lines& lines) {
    for (const auto& [key, value] : lines) out << key << " = " << value << '\n';
}

void blank(std::ostringstream& out) { out << '\n'; }

std::string fp(double x) { return std::to_string(x); }

constexpr int kNstCouple = 20;

int segment_ld_seed(const Config& cfg, int segment_index) {
    if (cfg.seed < 0) return -1;
    return cfg.seed + segment_index;
}

std::string production_mdp(const Config& cfg, const MdSpec& spec) {
    std::ostringstream out;
    out << "; " << spec.title << '\n';
    emit(out, {
        {"integrator", "md"},
        {"dt", fp(spec.dt_ps)},
        {"nsteps", std::to_string(spec.nsteps)},
        {"tinit", fp(spec.tinit_ps)},
        {"continuation", spec.continuation ? "yes" : "no"},
        {"comm-mode", "Linear"},
        {"nstcomm", "100"}
    });
    blank(out);
    emit(out, {
        {"nstxout", "0"},
        {"nstvout", "0"},
        {"nstenergy", std::to_string(cfg.nst_energy)},
        {"nstlog", std::to_string(cfg.nst_log)},
        {"nstxout-compressed", std::to_string(spec.nst_xtc > 0 ? spec.nst_xtc : cfg.nst_xtc)},
        {"compressed-x-precision", std::to_string(cfg.xtc_precision)}
    });
    blank(out);
    out << thermostat_block(cfg);
    emit(out, {{"ld-seed", std::to_string(segment_ld_seed(cfg, spec.segment_index))},
               {"nsttcouple", std::to_string(kNstCouple)}});
    blank(out);

    if (spec.pressure) {
        emit(out, {
            {"pcoupl", "C-rescale"},
            {"pcoupltype", "isotropic"},
            {"tau-p", fp(spec.tau_p)},
            {"nstpcouple", std::to_string(kNstCouple)},
            {"ref-p", "1.0"},
            {"compressibility", "3e-4"},
            {"refcoord-scaling", "all"}
        });
    } else {
        emit(out, {{"pcoupl", "no"}});
    }

    blank(out);
    emit(out, {
        {"constraints", "none"},
        {"constraint-algorithm", "Lincs"},
        {"lincs-order", "8"},
        {"lincs-iter", "2"},
        {"gen_vel", spec.gen_vel ? "yes" : "no"}
    });
    if (spec.gen_vel)
        emit(out, {{"gen_temp", fp(cfg.temperature_K)},
                   {"gen_seed", std::to_string(segment_ld_seed(cfg, spec.segment_index))}});
    blank(out);
    out << nonbonded_block(cfg);
    return out.str();
}

void write_md_spec(const Config& cfg, const MdSpec& spec) {
    if (cfg.dry_run) return;
    write_text(cfg.systemDir() / spec.filename, production_mdp(cfg, spec));
}

} // namespace cg
