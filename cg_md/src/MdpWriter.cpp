#include "MdpWriter.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"

#include <sstream>
#include <utility>
#include <vector>

namespace cg {
namespace {

using Lines = std::vector<std::pair<std::string, std::string>>;

void emit(std::ostringstream& out, const Lines& lines) {
    for (const auto& [key, value] : lines) out << key << " = " << value << '\n';
}

void blank(std::ostringstream& out) { out << '\n'; }

std::string fp(double x) { return std::to_string(x); }

struct MdSpec {
    std::string title;
    std::filesystem::path filename;
    double dt_ps = 0.0;
    std::int64_t nsteps = 0;
    bool continuation = true;
    bool pressure = true;
    bool gen_vel = false;
    double tau_p = 12.0;
    std::string coulomb_modifier = "Potential-shift-verlet";
};

std::string production_mdp(const Config& cfg, const MdSpec& spec) {
    std::ostringstream out;
    out << "; " << spec.title << '\n';
    emit(out, {
        {"integrator", "md"},
        {"dt", fp(spec.dt_ps)},
        {"nsteps", std::to_string(spec.nsteps)},
        {"continuation", spec.continuation ? "yes" : "no"},
        {"comm-mode", "Linear"}
    });
    blank(out);
    emit(out, {
        {"nstxout", "0"},
        {"nstvout", "0"},
        {"nstenergy", "1000"},
        {"nstlog", "1000"},
        {"nstxout-compressed", "1000"}
    });
    blank(out);
    out << thermostat_block(cfg);
    blank(out);

    if (spec.pressure) {
        emit(out, {
            {"pcoupl", "C-rescale"},
            {"pcoupltype", "isotropic"},
            {"tau-p", fp(spec.tau_p)},
            {"ref-p", "1.0"},
            {"compressibility", "3e-4"}
        });
    } else {
        emit(out, {{"pcoupl", "no"}});
    }

    blank(out);
    emit(out, {
        {"constraints", "none"},
        {"gen_vel", spec.gen_vel ? "yes" : "no"}
    });
    if (spec.gen_vel) emit(out, {{"gen_temp", fp(cfg.temperature_K)}, {"gen_seed", "-1"}});
    blank(out);
    emit(out, {
        {"cutoff-scheme", "Verlet"},
        {"coulombtype", "Reaction-Field"},
        {"rvdw", "1.1"},
        {"rcoulomb", "1.1"},
        {"coulomb-modifier", spec.coulomb_modifier}
    });
    return out.str();
}

void write_md_spec(const Config& cfg, const MdSpec& spec) {
    write_text(cfg.systemDir() / spec.filename, production_mdp(cfg, spec));
}

} // namespace

std::string thermostat_block(const Config& cfg) {
    std::ostringstream out;
    emit(out, {{"tcoupl", "V-rescale"}});

    if (cfg.thermostat_mode == "legacy") {
        emit(out, {{"tc-grps", "Protein W ION"}, {"tau-t", "1.0 1.0 1.0"},
                   {"ref-t", fp(cfg.temperature_K) + " " + fp(cfg.temperature_K) + " " + fp(cfg.temperature_K)}});
    } else if (cfg.thermostat_mode == "protein-solvent") {
        emit(out, {{"tc-grps", "Protein Solvent_and_ions"}, {"tau-t", "1.0 1.0"},
                   {"ref-t", fp(cfg.temperature_K) + " " + fp(cfg.temperature_K)}});
    } else {
        emit(out, {{"tc-grps", "System"}, {"tau-t", "1.0"}, {"ref-t", fp(cfg.temperature_K)}});
    }
    return out.str();
}

void write_em_mdp(const Config& cfg) {
    std::ostringstream out;
    emit(out, {
        {"integrator", "steep"},
        {"emtol", "100"},
        {"emstep", fp(cfg.em_emstep)},
        {"nsteps", std::to_string(cfg.em_nsteps)}
    });
    blank(out);
    emit(out, {
        {"cutoff-scheme", "Verlet"},
        {"nstlist", "10"},
        {"vdwtype", "Cut-off"},
        {"rvdw", "1.1"},
        {"coulombtype", "reaction-field"},
        {"rcoulomb", "1.1"}
    });
    blank(out);
    emit(out, {{"constraints", "none"}, {"nstenergy", "10"}});
    write_text(cfg.systemDir() / "emin.mdp", out.str());
}

void write_nvt_mdp(const Config& cfg) {
    write_md_spec(cfg, {"NVT equilibration for Martini coarse-grained system", "nvt.mdp",
                        cfg.nvt_dt_ps, cfg.nvt_nsteps, false, false, true});
}

void write_npt_mdp(const Config& cfg) {
    write_md_spec(cfg, {"NPT equilibration for Martini coarse-grained system", "npt.mdp",
                        cfg.npt_dt_ps, cfg.npt_nsteps, true, true, false, 4.0, "Potential-shift"});
}

void write_md_mdp(const Config& cfg) {
    write_md_spec(cfg, {"Production MD for Martini coarse-grained system; total time = " + std::to_string(cfg.md_total_us) + " us",
                        "md.mdp", cfg.md_dt_ps, cfg.mdNsteps()});
}

void write_metad_mdp(const Config& cfg) {
    write_md_spec(cfg, {"Metadynamics MD for Martini coarse-grained system",
                        "md_metad.mdp", cfg.md_dt_ps, cfg.mdNsteps()});
}

void write_md_chunk_mdp(const Config& cfg, std::int64_t chunk_nsteps) {
    write_md_spec(cfg, {"Adaptive chunk MD for Martini coarse-grained system; chunk steps = " + std::to_string(chunk_nsteps),
                        "md_chunk.mdp", cfg.md_dt_ps, chunk_nsteps});
}

} // namespace cg
