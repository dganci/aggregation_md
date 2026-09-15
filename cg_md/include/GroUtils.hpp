#pragma once

#include <array>
#include <cstddef>
#include <filesystem>

namespace cg {

struct GroNormalizationStats {
    std::size_t atom_names = 0;
    std::size_t residue_names = 0;
    std::size_t total() const { return atom_names + residue_names; }
};

/// insane/genion/gmx solvate sometimes write ion residue/atom names as
/// "NA"/"CL" instead of the "NA+"/"CL-" spelling Martini's ion topology
/// (martini_v3.0.0_ions_v1.itp) expects.
GroNormalizationStats normalize_gro_ion_names(const std::filesystem::path& gro_path);

/// The three box vectors (nm) from the last line of a .gro file.
std::array<std::array<double, 3>, 3> read_gro_box_vectors(const std::filesystem::path& gro_path);

/// Shortest non-zero lattice vector of `box`, i.e. how far a particle is from
/// its own nearest periodic image.
double min_image_distance_nm(const std::array<std::array<double, 3>, 3>& box);

/// True when a chain of longest internal distance `dmax_nm` cannot reach its
/// own periodic image: `image_nm >= dmax_nm + 2 * rcut_nm`.
bool periodic_margin_ok(double dmax_nm, double image_nm, double rcut_nm);

} // namespace cg
