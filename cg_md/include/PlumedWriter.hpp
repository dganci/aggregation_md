#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace cg {

struct Config;

struct PlumedFiles {
    std::filesystem::path plumed_dat;
    std::filesystem::path feature_schema_json;
    std::vector<std::string> feature_cols;
};

/// The set of PLUMED collective-variable definitions shared by *every* stage of
/// the pipeline, and the canonical names/order of the descriptors they produce.
struct PlumedDescriptors {
    /// The PLUMED input text defining everything below (groups, WHOLEMOLECULES,
    /// per-protomer COM, pairwise COM distances, pairwise COORDINATION, their
    /// sum, per-protomer and global radius of gyration).
    std::string text;

    std::vector<std::string> distances;
    std::vector<std::string> contacts;
    std::vector<std::string> per_chain_rg;

    /// Canonical feature order: distances, contacts, rg1..rgN, rg_com,
    /// cn_total.
    std::vector<std::string> feature_cols;
};

/// PLUMED input defining the permutation-invariant form of the descriptor set,
/// plus the names of the sorted components in canonical order.
struct SortedDescriptors {
    std::string text;
    std::vector<std::string> feature_cols;
};

/// Builds the SORT block for `d`.
SortedDescriptors build_sorted_descriptors(const PlumedDescriptors& d);

/// Builds the shared descriptor block.
PlumedDescriptors build_descriptors(const Config& cfg, const std::filesystem::path& ndx_path);

/// Writes an unbiased-production plumed.dat: the shared descriptor block plus a
/// PRINT of every column in PlumedDescriptors::feature_cols to `colvar_path`
/// every Config::plumed_stride steps.
PlumedFiles write_plumed_dat(const Config& cfg);
PlumedFiles write_plumed_dat(const Config& cfg,
                             const std::filesystem::path& plumed_path,
                             const std::filesystem::path& colvar_path);

/// Options controlling one metadynamics batch's plumed.dat.
struct MetadPlumedOptions {
    /// Emits a top-level `RESTART` directive.
    bool restart = false;

    /// Path (relative to the mdrun working directory) of the HILLS file.
    std::string hills_file = "HILLS";
};

/// Writes a metadynamics plumed.dat: the same shared descriptor block as
/// write_plumed_dat().
PlumedFiles write_metad_plumed_dat(const Config& cfg,
                                   const std::filesystem::path& plumed_path,
                                   const std::filesystem::path& ndx_path,
                                   const std::filesystem::path& model_path,
                                   const std::filesystem::path& colvar_path,
                                   const MetadPlumedOptions& options = {});

} // namespace cg
