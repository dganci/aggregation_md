#include "Workflow.hpp"
#include "Dssp.hpp"
#include "FileUtils.hpp"
#include "GroUtils.hpp"
#include "StringUtils.hpp"
#include "TopologyEditor.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cg {

void Workflow::ensure_dirs() {
    for (const auto& dir : {cfg_.pdbPath().parent_path(), cfg_.cgPath().parent_path(),
                            cfg_.topologyPath().parent_path(), cfg_.packRawPath().parent_path(),
                            cfg_.packCleanPath().parent_path(), cfg_.packInputPath().parent_path(),
                            cfg_.beadsPath().parent_path(), cfg_.systemDir(), cfg_.resultDir()}) {
        std::filesystem::create_directories(dir);
    }
}

void Workflow::check_inputs() {
    if ((cfg_.stage == "all" || cfg_.stage == "prepare") && !std::filesystem::exists(cfg_.pdbPath()))
        throw std::runtime_error("Missing input PDB: " + cfg_.pdbPath().string());
}

void Workflow::coarse_grain() {
    std::string ss;
    try {
        ss = dssp_to_martinize_ss(run_mkdssp(sh_, cfg_), cfg_.seq_length);
    } catch (const std::exception& e) {
        std::cerr << "Warning: DSSP parsing failed; martinize will use -dssp. Details: " << e.what() << '\n';
    }

    std::vector<std::string> cmd = {
        cfg_.martinize, "-f", cfg_.pdbPath().string(), "-x", cfg_.cgPath().string(),
        "-o", cfg_.topologyPath().string(), "-name", cfg_.protomer_name,
        "-ff", "martini3IDP", "-p", "backbone", "-cys", "auto"
    };

    if (cfg_.use_elastic) {
        cmd.insert(cmd.end(), {"-elastic", "-ef", "700.0", "-el", "0", "-eu", "0.85", "-eunit", cfg_.elastic_units});
    }
    if (!ss.empty() && !cfg_.is_phospho && cfg_.use_dssp_ss_string) {
        cmd.insert(cmd.end(), {"-ss", ss});
    } else {
        cmd.push_back("-dssp");
    }

    sh_.run(cmd, cfg_.project_dir);
    move_generated_itps();
}

void Workflow::move_generated_itps() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping ITP move after martinize2.\n";
        return;
    }

    for (const auto& entry : std::filesystem::directory_iterator(cfg_.project_dir)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".itp") continue;
        const auto dst = cfg_.topologyPath().parent_path() / entry.path().filename();
        remove_if_exists(dst);
        std::filesystem::rename(entry.path(), dst);
    }
}

void Workflow::write_packmol_input() {
    const auto r = cfg_.packmol_cluster_radius_A;
    std::ostringstream out;
    out << "tolerance " << cfg_.packmol_tolerance_A << '\n'
        << "filetype pdb\n"
        << "structure " << cfg_.cgPath().string() << '\n'
        << "  number " << cfg_.n_prot << '\n'
        << "  inside sphere " << r << ' ' << r << ' ' << r << ' ' << r << '\n'
        << "  radius " << cfg_.packmol_radius_A << '\n'
        << "end structure\n"
        << "output " << cfg_.packRawPath().string() << '\n';
    write_text(cfg_.packInputPath(), out.str());
}

void Workflow::run_packmol() {
    sh_.runShell(cfg_.packmol + " < " + shell_quote(cfg_.packInputPath().string()));
}

void Workflow::clean_packmol_pdb() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping Packmol PDB cleaning because " << cfg_.packRawPath() << " was not generated.\n";
        return;
    }

    std::ifstream in(cfg_.packRawPath());
    if (!in) throw std::runtime_error("Cannot read " + cfg_.packRawPath().string());

    std::ostringstream out;
    for (std::string line; std::getline(in, line);) {
        if (starts_with(line, "CONECT")) continue;
        if (starts_with(line, "HETATM")) line.replace(0, 6, "ATOM  ");
        out << line << '\n';
    }
    write_text(cfg_.packCleanPath(), out.str());
}

void Workflow::run_insane() {
    sh_.run({cfg_.insane, "-o", cfg_.beadsPath().string(), "-p", cfg_.topologyPath().string(),
             "-f", cfg_.packCleanPath().string(), "-sol", "W", "-salt", std::to_string(cfg_.salt_M),
             "-pbc", "optimal", "-d", std::to_string(cfg_.box_margin_nm)});
    normalize_ions();
}

void Workflow::make_compact_box() {
    sh_.run({cfg_.gmx, "editconf", "-f", cfg_.packCleanPath().string(), "-o", cfg_.boxedPath().string(),
             "-bt", cfg_.box_type, "-d", std::to_string(cfg_.box_margin_nm), "-c"});
}

void Workflow::solvate_compact_box() {
    if (!sh_.dryRun() && !std::filesystem::exists(cfg_.martiniWaterPath())) {
        throw std::runtime_error("Missing Martini water coordinate file: " + cfg_.martiniWaterPath().string() +
                                 "\nProvide it with --martini-water-gro PATH");
    }

    sh_.run({cfg_.gmx, "solvate", "-cp", cfg_.boxedPath().string(), "-cs", cfg_.martiniWaterPath().string(),
             "-radius", std::to_string(cfg_.solvate_radius_nm), "-o", cfg_.solvatedPath().string(),
             "-p", cfg_.topologyPath().string()});
    if (!sh_.dryRun()) copy_overwrite(cfg_.solvatedPath(), cfg_.beadsPath());
}

std::string Workflow::make_ndx_and_find_water_group(const std::filesystem::path& gro_file,
                                                    const std::filesystem::path& ndx_path,
                                                    const std::filesystem::path& log_path) const {
    sh_.runShell("printf 'q\\n' | " + shell_quote(cfg_.gmx) +
                 " make_ndx -f " + shell_quote(gro_file.string()) +
                 " -o " + shell_quote(ndx_path.string()) +
                 " > " + shell_quote(log_path.string()));

    if (sh_.dryRun()) return "W";

    for (const auto& line : read_lines(log_path)) {
        const auto tokens = split_ws(line);
        if (tokens.size() >= 2 && tokens[1] == "W") return tokens[0];
    }
    throw std::runtime_error("Cannot identify W group for genion.");
}

void Workflow::run_genion(const std::filesystem::path& tpr_path,
                          const std::filesystem::path& ndx_path,
                          const std::string& water_group,
                          bool include_salt) {
    std::ostringstream cmd;
    cmd << "printf '" << water_group << "\\n' | " << shell_quote(cfg_.gmx)
        << " genion -s " << shell_quote(tpr_path.string())
        << " -n " << shell_quote(ndx_path.string())
        << " -o " << shell_quote(cfg_.beadsPath().string())
        << " -p " << shell_quote(cfg_.topologyPath().string())
        << " -pname NA+ -nname CL- -neutral";
    if (include_salt && cfg_.salt_M > 0.0) cmd << " -conc " << cfg_.salt_M;
    sh_.runShell(cmd.str());
    normalize_ions();
}

void Workflow::add_ions_to_solvated_system() {
    if (!cfg_.add_ions) {
        normalize_ions();
        return;
    }

    const auto ndx = cfg_.resultDir() / "genion_tmp.ndx";
    const auto log = cfg_.resultDir() / "make_ndx.out";
    const auto water = make_ndx_and_find_water_group(cfg_.beadsPath(), ndx, log);
    grompp(cfg_.systemDir() / "emin.mdp", cfg_.beadsPath(), cfg_.ionsTprPath(), {}, {}, ndx);
    run_genion(cfg_.ionsTprPath(), ndx, water, true);
    if (!sh_.dryRun()) { remove_if_exists(ndx); remove_if_exists(log); }
}

void Workflow::prepare_solvated_system() {
    if (cfg_.solvation_mode == "insane") {
        run_insane();
        patch_topology();
        write_all_mdps();
        return;
    }

    patch_topology();
    write_all_mdps();
    make_compact_box();
    solvate_compact_box();
    add_ions_to_solvated_system();
}

void Workflow::normalize_ions() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping GRO ion name normalization.\n";
        return;
    }

    const auto stats = normalize_gro_ion_names(cfg_.beadsPath());
    if (stats.total()) {
        std::cerr << "Normalized GRO ion names in " << cfg_.beadsPath() << ": "
                  << stats.atom_names << " atom name field(s), "
                  << stats.residue_names << " residue name field(s).\n";
    } else {
        std::cerr << "GRO ion names already match Martini ion names.\n";
    }
}

void Workflow::neutralize_if_needed() {
    if (!cfg_.is_phospho || !cfg_.neutralize_if_phospho) return;

    const auto ndx = cfg_.resultDir() / "genion_tmp.ndx";
    const auto log = cfg_.resultDir() / "make_ndx.out";
    const auto water = make_ndx_and_find_water_group(cfg_.beadsPath(), ndx, log);
    run_genion(cfg_.resultDir() / "em.tpr", ndx, water, false);
    if (!sh_.dryRun()) { remove_if_exists(ndx); remove_if_exists(log); }
}

void Workflow::patch_topology() {
    if (sh_.dryRun()) {
        std::cerr << "[dry-run] Skipping topology patch because martinize2 was not executed.\n";
        return;
    }

    patch_martini_topology(cfg_.topologyPath(), cfg_.protomer_name, cfg_.n_prot,
                           find_first_itp(cfg_.topologyPath().parent_path(), cfg_.protomer_name));
}

} // namespace cg
