#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace cg {


/// Three-letter PDB residue codes for the phosphorylated amino acids this
/// pipeline knows how to coarse-grain (phosphoserine, phosphothreonine,
/// phosphotyrosine).
extern const std::vector<std::string> kKnownModifiedResidues;

/// Returns the first known modified-residue code found in `pdb_path`'s ATOM/
/// HETATM records (fixed-width PDB columns 18-20), or an empty string if none.
std::string detect_modified_residue(const std::filesystem::path& pdb_path);

/// Every DISTINCT modified residue in the structure, in order of first
/// appearance - not just the first one.
std::vector<std::string> detect_modified_residues(const std::filesystem::path& pdb_path);

/// Counts distinct residues found in `pdb_path`'s ATOM/HETATM records, i.e.
/// distinct (chain id, resSeq, insertion code) triples in fixed-width PDB
/// columns [21,27).
int count_pdb_residues(const std::filesystem::path& pdb_path);

/// Counts ATOM/HETATM records in the first model of `pdb_path`.
int count_pdb_atoms(const std::filesystem::path& pdb_path);

/// Distinct chain identifiers (PDB column 22) present in the first model.
std::vector<std::string> pdb_chain_ids(const std::filesystem::path& pdb_path);

/// Largest per-molecule maximum interatomic distance (Dmax, in the PDB's own
/// units - angstrom) in `pdb_path`, splitting its ATOM/HETATM records into
/// consecutive blocks of `atoms_per_molecule`.
double pdb_max_molecule_dmax(const std::filesystem::path& pdb_path, int atoms_per_molecule);

/// Smallest and largest residue sequence number (PDB columns 23-26) in the
/// first model, or {0, 0} when the file has no ATOM/HETATM records.
std::pair<int, int> pdb_residue_number_range(const std::filesystem::path& pdb_path);


} // namespace cg
