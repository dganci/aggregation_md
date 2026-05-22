#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace cg {

struct Config {
    std::filesystem::path project_dir = ".";

    std::string gmx = "gmx_mpi";
    std::string mkdssp = "mkdssp";
    std::string martinize = "martinize2";
    std::string insane = "insane";
    std::string packmol = "packmol";
    std::string plumed = "plumed";
    std::string mpirun = "mpirun";
    std::string python = "python3";

    int n_prot = 10;
    std::string protomer_name = "desmin_head";
    int seq_length = 108;
    bool is_phospho = false;
    int atoms_per_prot = 237;

    double packmol_box_A = 145.0;
    double packmol_cluster_radius_A = 120.0;
    double packmol_tolerance_A = 2.5;
    double packmol_radius_A = 5.0;

    std::string solvation_mode = "gromacs";
    std::string box_type = "dodecahedron";
    double box_margin_nm = 1.5;
    double solvate_radius_nm = 0.21;
    std::filesystem::path martini_water_gro = "system/water.gro";
    bool add_ions = true;

    double temperature_K = 300.0;
    double salt_M = 0.15;

    double em_emstep = 0.01;
    int em_nsteps = 50000;

    double nvt_dt_ps = 0.01;
    std::int64_t nvt_nsteps = 100000;

    double npt_dt_ps = 0.01;
    std::int64_t npt_nsteps = 1000000;

    double md_dt_ps = 0.01;
    double md_total_us = 10.0;

    int plumed_stride = 100;
    double contact_r0_nm = 0.55;
    int ntomp = 0;

    bool metadynamics = false;
    int metad_walkers = 1;
    int metad_nodes = 3;
    int metad_pace = 500;
    int metad_print_stride = 100;
    double metad_height = 0.8;
    double metad_biasfactor = 10.0;
    std::string metad_sigma;
    std::filesystem::path metad_model;
    std::filesystem::path metad_cv_params;
    bool metad_sum_hills = true;
    bool metad_center = true;

    bool adaptive = false;
    double adaptive_chunk_us = 0.25;
    double adaptive_max_total_us = 5.0;
    double adaptive_min_total_us = 0.5;
    int adaptive_min_cn_state_transitions = 10;
    int adaptive_min_bidirectional_events = 1;
    double adaptive_min_cn_range = 20.0;
    double adaptive_min_rg_global_range = 0.5;
    int adaptive_min_unique_contact_patterns = 10;
    int adaptive_min_largest_cluster_unique = 3;
    double adaptive_pair_contact_threshold = 1.0;
    int adaptive_smooth_window_frames = 50;
    int adaptive_min_residence_frames = 5;
    int adaptive_pattern_downsample = 10;

    std::string thermostat_mode = "system";
    std::string elastic_units = "90:102";
    bool use_elastic = true;
    bool use_dssp_ss_string = true;
    bool neutralize_if_phospho = true;

    bool dry_run = false;
    std::string stage = "all";

    std::string systemName() const;
    std::int64_t mdNsteps() const;
    std::int64_t adaptiveChunkNsteps() const;
    int adaptiveMaxChunks() const;

    std::filesystem::path pdbPath() const;
    std::filesystem::path cgPath() const;
    std::filesystem::path topologyPath() const;
    std::filesystem::path packRawPath() const;
    std::filesystem::path packCleanPath() const;
    std::filesystem::path packInputPath() const;
    std::filesystem::path boxedPath() const;
    std::filesystem::path solvatedPath() const;
    std::filesystem::path ionsTprPath() const;
    std::filesystem::path martiniWaterPath() const;
    std::filesystem::path beadsPath() const;
    std::filesystem::path systemDir() const;
    std::filesystem::path resultDir() const;
    std::filesystem::path metadDir() const;
    std::filesystem::path metadModelPath() const;
    std::filesystem::path metadCvParamsPath() const;

    static Config fromArgs(int argc, char** argv);
    void validate() const;
    void print() const;
};

} // namespace cg
