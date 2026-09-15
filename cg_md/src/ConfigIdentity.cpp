#include "Config.hpp"
#include "StringUtils.hpp"

#include <cstdint>
#include <sstream>
#include <string>

namespace cg {

std::string Config::identityString() const {
    std::ostringstream o;
    const auto b = [](bool v) { return v ? "true" : "false"; };
    o << "n_prot = " << n_prot << '\n'
      << "protomer = " << protomer_name << '\n'
      << "seq_length = " << seq_length << '\n'
      << "atoms_per_prot = " << atoms_per_prot << '\n'
      << "phospho = " << b(is_phospho) << '\n'
      << "neutralize_if_phospho = " << b(neutralize_if_phospho) << '\n'
      << "ss_string = " << ss_string << '\n'
      << "use_dssp = " << b(use_dssp) << '\n'
      << "merge_chains = " << merge_chains << '\n'
      << "modify = " << modify << '\n'
      << "use_elastic = " << b(use_elastic) << '\n'
      << "elastic_units = " << (use_elastic ? elastic_units : std::string{}) << '\n'
      << "relax = " << b(relax) << '\n'
      << "relax_us = " << to_string_fixed(relax ? relax_us : 0.0, 6) << '\n'
      << "relax_frames = " << (relax ? relax_frames : 0) << '\n'
      << "relax_cluster_cutoff_nm = " << to_string_fixed(relax ? relax_cluster_cutoff_nm : 0.0, 4) << '\n'
      << "relax_discard_frac = " << to_string_fixed(relax ? relax_discard_frac : 0.0, 4) << '\n'
      << "packmol_box_A = " << to_string_fixed(packmol_box_A, 3) << '\n'
      << "packmol_cluster_radius_A = " << to_string_fixed(packmol_cluster_radius_A, 3) << '\n'
      << "packmol_tolerance_A = " << to_string_fixed(packmol_tolerance_A, 3) << '\n'
      << "packmol_radius_A = " << to_string_fixed(packmol_radius_A, 3) << '\n'
      << "solvation_mode = " << solvation_mode << '\n'
      << "box_type = " << box_type << '\n'
      << "box_margin_nm = " << to_string_fixed(box_margin_nm, 4) << '\n'
      << "solvate_radius_nm = " << to_string_fixed(solvate_radius_nm, 4) << '\n'
      << "add_ions = " << b(add_ions) << '\n'
      << "temperature_K = " << to_string_fixed(temperature_K, 3) << '\n'
      << "salt_M = " << to_string_fixed(salt_M, 4) << '\n'
      << "martini_lambda_pw = " << to_string_fixed(martini_lambda_pw, 4) << '\n'
      << "distinct_copies = " << b(distinct_copies) << '\n'
      << "epsilon_r = " << to_string_fixed(epsilon_r, 4) << '\n'
      << "rcoulomb_nm = " << to_string_fixed(rcoulomb_nm, 4) << '\n'
      << "rvdw_nm = " << to_string_fixed(rvdw_nm, 4) << '\n'
      << "nstlist = " << nstlist << '\n'
      << "seed = " << seed << '\n'
      << "em_emstep = " << to_string_fixed(em_emstep, 6) << '\n'
      << "em_nsteps = " << em_nsteps << '\n'
      << "nvt_dt_ps = " << to_string_fixed(nvt_dt_ps, 6) << '\n'
      << "nvt_nsteps = " << nvt_nsteps << '\n'
      << "npt_dt_ps = " << to_string_fixed(npt_dt_ps, 6) << '\n'
      << "npt_nsteps = " << npt_nsteps << '\n'
      << "md_dt_ps = " << to_string_fixed(md_dt_ps, 6) << '\n'
      << "thermostat_mode = " << thermostat_mode << '\n'
      << "plumed_stride = " << plumed_stride << '\n'
      << "contact_r0_nm = " << to_string_fixed(contact_r0_nm, 4) << '\n'
      << "contact_nn = " << contact_nn << '\n'
      << "contact_mm = " << contact_mm << '\n'
      << "nst_xtc = " << nst_xtc << '\n'
      << "nst_energy = " << nst_energy << '\n'
      << "nst_log = " << nst_log << '\n'
      << "adaptive = " << b(adaptive) << '\n'
      << "adaptive_chunk_us = " << to_string_fixed(adaptive ? adaptive_chunk_us : 0.0, 8) << '\n'
      << "metadynamics = " << b(metadynamics) << '\n'
      << "pmf = " << b(pmf) << '\n';
    if (metadynamics || pmf) {
        o << "metad_chunk_us = " << to_string_fixed(metad_chunk_us, 8) << '\n'
          << "metad_walkers = " << metad_walkers << '\n'
          << "metad_nodes = " << metad_nodes << '\n'
          << "metad_pace = " << metad_pace << '\n'
          << "metad_height = " << to_string_fixed(metad_height, 4) << '\n'
          << "metad_biasfactor = " << to_string_fixed(metad_biasfactor, 4) << '\n'
          << "metad_sigma = " << metad_sigma << '\n'
          << "metad_feature_cols = " << metad_feature_cols << '\n'
          << "metad_permutation_invariant = " << b(metad_permutation_invariant) << '\n'
          << "metad_grid = " << b(metad_grid) << '\n'
          << "metad_grid_min = " << metad_grid_min << '\n'
          << "metad_grid_max = " << metad_grid_max << '\n'
          << "metad_grid_bin = " << metad_grid_bin << '\n'
          << "metad_wall_kt = " << to_string_fixed(metad_wall_kt, 4) << '\n'
          << "metad_wall_margin_frac = " << to_string_fixed(metad_wall_margin_frac, 4) << '\n';
    }
    if (pmf) {
        o << "pmf_wall_nm = " << to_string_fixed(pmf_wall_nm, 4) << '\n'
          << "pmf_wall_kappa = " << to_string_fixed(pmf_wall_kappa, 4) << '\n'
          << "pmf_sigma_nm = " << to_string_fixed(pmf_sigma_nm, 4) << '\n';
    }
    return o.str();
}

std::string Config::runId() const {
    if (!run_tag.empty()) return run_tag;
    if (!frozen_run_id.empty()) return frozen_run_id;
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : identityString()) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    std::ostringstream o;
    o << std::hex << std::setfill('0') << std::setw(8) << static_cast<std::uint32_t>(h >> 32);
    frozen_run_id = o.str();
    return frozen_run_id;
}

std::string Config::runName() const { return systemName() + "_" + runId(); }

std::filesystem::path Config::runDir() const { return project_dir / "runs" / runName(); }
std::filesystem::path Config::prepDir() const { return runDir() / "prep"; }

std::filesystem::path Config::pdbPath() const { return project_dir / "pdbs" / "AA" / (protomer_name + ".pdb"); }
std::filesystem::path Config::cgPath() const { return prepDir() / (protomer_name + "_cg.pdb"); }
std::filesystem::path Config::cgRelaxedPath() const { return prepDir() / (protomer_name + "_cg_relaxed.pdb"); }

} // namespace cg
