#pragma once

#include <cstddef>
#include <filesystem>

namespace cg {

struct GroNormalizationStats {
    std::size_t atom_names = 0;
    std::size_t residue_names = 0;
    std::size_t total() const { return atom_names + residue_names; }
};

GroNormalizationStats normalize_gro_ion_names(const std::filesystem::path& gro_path);

} // namespace cg
