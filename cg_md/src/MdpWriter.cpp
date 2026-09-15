#include "MdpWriter.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "MdpSpec.hpp"
#include "StringUtils.hpp"

#include <cstdint>
#include <sstream>
#include <string>

namespace cg {

std::string nonbonded_block(const Config& cfg) {
    std::ostringstream out;
    emit(out, {
        {"cutoff-scheme", "Verlet"},
        {"nstlist", std::to_string(cfg.nstlist)},
        {"verlet-buffer-tolerance", to_string_fixed(cfg.verlet_buffer_tolerance, 4)},
        {"rlist", to_string_fixed(cfg.rlist_nm, 2)},
        {"coulombtype", "Reaction-Field"},
        {"rcoulomb", to_string_fixed(cfg.rcoulomb_nm, 2)},
        {"epsilon_r", to_string_fixed(cfg.epsilon_r, 1)},
        {"epsilon_rf", "0"},
        {"coulomb-modifier", "Potential-shift"},
        {"vdwtype", "Cut-off"},
        {"vdw-modifier", "Potential-shift"},
        {"rvdw", to_string_fixed(cfg.rvdw_nm, 2)}
    });
    return out.str();
}

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
    out << nonbonded_block(cfg);
    blank(out);
    emit(out, {
        {"constraints", "none"},
        {"constraint-algorithm", "Lincs"}, {"lincs-order", "8"}, {"lincs-iter", "2"},
        {"lincs-warnangle", "90"},
        {"nstenergy", "10"}});
    if (cfg.dry_run) return;
    write_text(cfg.systemDir() / "emin.mdp", out.str());
}

void write_nvt_mdp(const Config& cfg) {
    write_md_spec(cfg, {"NVT equilibration for Martini coarse-grained system", "nvt.mdp",
                        cfg.nvt_dt_ps, cfg.nvt_nsteps, false, false, true,
                        kTauP, 0.0, kSeedNvt});
}

void write_npt_mdp(const Config& cfg) {
    write_md_spec(cfg, {"NPT equilibration for Martini coarse-grained system", "npt.mdp",
                        cfg.npt_dt_ps, cfg.npt_nsteps, true, true, false, kTauP,
                        0.0, kSeedNpt});
}

void write_md_mdp(const Config& cfg) {
    write_md_spec(cfg, {"Production MD for Martini coarse-grained system; total time = " + std::to_string(cfg.md_total_us) + " us",
                        "md.mdp", cfg.md_dt_ps, cfg.mdNsteps(),
                        true, true, false, kTauP, 0.0, kSeedProduction});
}

std::filesystem::path write_md_chunk_mdp(const Config& cfg, int chunk,
                                         std::int64_t chunk_nsteps, double tinit_ps) {
    const std::filesystem::path name = "md_chunk_" + zero_padded(chunk) + ".mdp";
    write_md_spec(cfg, {"Adaptive chunk " + std::to_string(chunk) + " MD; steps = " +
                            std::to_string(chunk_nsteps) + ", tinit = " + fp(tinit_ps) + " ps",
                        name, cfg.md_dt_ps, chunk_nsteps,
                        /*continuation=*/true, /*pressure=*/true, /*gen_vel=*/false,
                        /*tau_p=*/kTauP,
                        /*tinit_ps=*/tinit_ps, /*segment_index=*/kSeedAdaptiveBase + chunk});
    return cfg.systemDir() / name;
}

std::filesystem::path write_relax_mdp(const Config& cfg, std::int64_t nsteps) {
    const auto stride = static_cast<int>(std::max<std::int64_t>(1, nsteps / std::max(1, cfg.relax_frames)));
    const std::filesystem::path name = "md_relax.mdp";
    write_md_spec(cfg, {"Single-protomer relaxation before packing; steps = " + std::to_string(nsteps) +
                            ", ~" + std::to_string(cfg.relax_frames) + " frames",
                        name, cfg.md_dt_ps, nsteps,
                        /*continuation=*/true, /*pressure=*/true, /*gen_vel=*/false,
                        /*tau_p=*/kTauP, /*tinit_ps=*/0.0, /*segment_index=*/kSeedRelax,
                        /*nst_xtc=*/stride});
    return cfg.systemDir() / name;
}

std::filesystem::path write_metad_batch_mdp(const Config& cfg, int batch,
                                            std::int64_t batch_nsteps, double tinit_ps,
                                            int walker) {
    const std::filesystem::path name =
        "md_metad_" + zero_padded(batch) + (walker > 0 ? "_w" + std::to_string(walker) : "") + ".mdp";
    write_md_spec(cfg, {"Metadynamics batch " + std::to_string(batch) + " MD; steps = " +
                            std::to_string(batch_nsteps) + ", tinit = " + fp(tinit_ps) + " ps",
                        name, cfg.md_dt_ps, batch_nsteps,
                        /*continuation=*/true, /*pressure=*/true, /*gen_vel=*/false,
                        /*tau_p=*/kTauP,
                        /*tinit_ps=*/tinit_ps,
                        /*segment_index=*/kSeedMetadBase + batch + walker * kSeedWalkerStride});
    return cfg.systemDir() / name;
}

} // namespace cg
