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

PlumedFiles write_plumed_dat(const Config& cfg);
PlumedFiles write_plumed_dat(const Config& cfg,
                             const std::filesystem::path& plumed_path,
                             const std::filesystem::path& colvar_path);

PlumedFiles write_metad_plumed_dat(const Config& cfg,
                                   const std::filesystem::path& plumed_path,
                                   const std::filesystem::path& ndx_path,
                                   const std::filesystem::path& model_path,
                                   const std::filesystem::path& colvar_path);

} // namespace cg
