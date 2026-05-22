#pragma once

#include "Config.hpp"
#include "SamplingMonitor.hpp"
#include "Shell.hpp"

#include <filesystem>
#include <initializer_list>
#include <string>
#include <vector>

namespace cg {

class Workflow {
public:
    explicit Workflow(Config cfg);
    void run();

private:
    using Step = void (Workflow::*)();

    Config cfg_;
    Shell sh_;

    void run_all();
    void run_steps(std::initializer_list<Step> steps);

    void ensure_dirs();
    void check_inputs();
    void coarse_grain();
    void move_generated_itps();
    void write_packmol_input();
    void run_packmol();
    void clean_packmol_pdb();
    void run_insane();
    void make_compact_box();
    void solvate_compact_box();
    void add_ions_to_solvated_system();
    void prepare_solvated_system();
    void neutralize_if_needed();
    void patch_topology();
    void write_all_mdps();
    void normalize_ions();

    std::string make_ndx_and_find_water_group(const std::filesystem::path& gro_file,
                                              const std::filesystem::path& ndx_path,
                                              const std::filesystem::path& log_path) const;
    void run_genion(const std::filesystem::path& tpr_path,
                    const std::filesystem::path& ndx_path,
                    const std::string& water_group,
                    bool include_salt);

    void prepare_em();
    void run_em();
    void analyze_em();
    void prepare_nvt();
    void run_nvt();
    void prepare_npt();
    void run_npt();
    void analyze_npt();
    void prepare_production();
    void run_production();
    void center_trajectory();

    void prepare_metadynamics();
    void run_metadynamics();
    void finalize_metadynamics();
    void run_metadynamics_workflow();

    void run_adaptive_production();
    SamplingRules adaptive_rules() const;
    std::filesystem::path chunk_base(int chunk) const;
    std::filesystem::path chunk_tpr(int chunk) const;
    std::filesystem::path chunk_colvar(int chunk) const;
    std::filesystem::path chunk_plumed(int chunk) const;
    bool chunk_completed(int chunk) const;
    std::vector<int> completed_chunks() const;
    std::vector<std::filesystem::path> completed_colvars() const;

    void grompp(const std::filesystem::path& mdp,
                const std::filesystem::path& coordinates,
                const std::filesystem::path& output,
                const std::filesystem::path& checkpoint = {},
                const std::filesystem::path& restraints = {},
                const std::filesystem::path& index = {});

    void mdrun(const std::filesystem::path& deffnm,
               const std::filesystem::path& tpr = {},
               const std::filesystem::path& plumed = {},
               int ntomp_override = 0);

    void trjconv(const std::filesystem::path& tpr,
                 const std::filesystem::path& in,
                 const std::filesystem::path& out,
                 const std::string& groups,
                 const std::vector<std::string>& options);

    void energy(const std::filesystem::path& edr,
                const std::filesystem::path& xvg,
                const std::string& observable);
};

} // namespace cg
