#include "test_framework.hpp"
#include "FileUtils.hpp"
#include "PdbUtils.hpp"

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

using namespace cg;

namespace {

struct ScratchPdb {
    std::filesystem::path path;

    explicit ScratchPdb(const std::string& name, const std::string& text)
        : path(std::filesystem::temp_directory_path() / ("cg_md_test_" + name + ".pdb")) {
        write_text(path, text);
    }
    ~ScratchPdb() { std::filesystem::remove(path); }

    ScratchPdb(const ScratchPdb&) = delete;
    ScratchPdb& operator=(const ScratchPdb&) = delete;
};

constexpr const char* kTwoResiduesOneChain =
    "ATOM      1  BB  MET A   1      12.345  1.000  2.000  1.00  0.00\n"
    "ATOM      2  SC1 MET A   1      13.345  1.000  2.000  1.00  0.00\n"
    "ATOM      3  BB  ALA A   2      14.345  1.000  2.000  1.00  0.00\n";

} // namespace

CG_TEST(count_pdb_atoms_counts_every_atom_record) {
    const ScratchPdb pdb("atoms", kTwoResiduesOneChain);
    CG_CHECK_EQ(count_pdb_atoms(pdb.path), 3);
}

CG_TEST(count_pdb_atoms_ignores_non_atom_records) {
    const ScratchPdb pdb("atoms_noise",
                         "HEADER \n"
                         "TITLE     Built with Packmol\n"
                         "REMARK   Packmol generated pdb file\n" +
                             std::string(kTwoResiduesOneChain) +
                             "TER\n"
                             "END\n");
    CG_CHECK_EQ(count_pdb_atoms(pdb.path), 3);
}

CG_TEST(count_pdb_atoms_counts_hetatm_as_well_as_atom) {
    const ScratchPdb pdb("atoms_hetatm",
                         "ATOM      1  BB  MET A   1      12.345  1.000  2.000  1.00  0.00\n"
                         "HETATM    2  P   SEP A   2      13.345  1.000  2.000  1.00  0.00\n");
    CG_CHECK_EQ(count_pdb_atoms(pdb.path), 2);
}

CG_TEST(count_pdb_atoms_stops_at_the_first_model) {
    const ScratchPdb pdb("atoms_models",
                         "MODEL        1\n" + std::string(kTwoResiduesOneChain) +
                             "ENDMDL\n"
                             "MODEL        2\n" +
                             std::string(kTwoResiduesOneChain) + "ENDMDL\n");
    CG_CHECK_EQ(count_pdb_atoms(pdb.path), 3);
}

CG_TEST(count_pdb_residues_counts_distinct_residues_not_atoms) {
    const ScratchPdb pdb("residues", kTwoResiduesOneChain);
    CG_CHECK_EQ(count_pdb_residues(pdb.path), 2);
}

CG_TEST(count_pdb_residues_spans_every_chain_in_the_model) {
    const ScratchPdb pdb("residues_chains",
                         "ATOM      1  BB  MET A   1      12.345  1.000  2.000  1.00  0.00\n"
                         "TER\n"
                         "ATOM      2  BB  MET B   1      22.345  1.000  2.000  1.00  0.00\n"
                         "ATOM      3  BB  ALA B   2      23.345  1.000  2.000  1.00  0.00\n");
    CG_CHECK_EQ(count_pdb_residues(pdb.path), 3);
}

CG_TEST(count_pdb_residues_distinguishes_insertion_codes) {
    const ScratchPdb pdb("residues_icode",
                         "ATOM      1  BB  MET A   1      12.345  1.000  2.000  1.00  0.00\n"
                         "ATOM      2  BB  ALA A   1A     13.345  1.000  2.000  1.00  0.00\n");
    CG_CHECK_EQ(count_pdb_residues(pdb.path), 2);
}

CG_TEST(pdb_chain_ids_lists_each_chain_once_in_order) {
    const ScratchPdb pdb("chains",
                         "ATOM      1  BB  MET A   1      12.345  1.000  2.000  1.00  0.00\n"
                         "ATOM      2  BB  MET B   1      22.345  1.000  2.000  1.00  0.00\n"
                         "ATOM      3  BB  ALA A   2      13.345  1.000  2.000  1.00  0.00\n");
    CG_CHECK_EQ(pdb_chain_ids(pdb.path), (std::vector<std::string>{"A", "B"}));
}

CG_TEST(pdb_chain_ids_ignores_a_blank_chain_column) {
    const ScratchPdb pdb("chains_blank",
                         "ATOM      1  BB  MET     1      12.345  1.000  2.000  1.00  0.00\n");
    CG_CHECK(pdb_chain_ids(pdb.path).empty());
}

CG_TEST(detect_modified_residue_finds_the_first_phospho_residue) {
    const ScratchPdb pdb("phospho",
                         "ATOM      1  BB  MET A   1      12.345  1.000  2.000  1.00  0.00\n"
                         "ATOM      2  BB  SEP A   2      13.345  1.000  2.000  1.00  0.00\n");
    CG_CHECK_EQ(detect_modified_residue(pdb.path), std::string("SEP"));
}

CG_TEST(detect_modified_residue_is_empty_for_an_unmodified_structure) {
    const ScratchPdb pdb("unmodified", kTwoResiduesOneChain);
    CG_CHECK(detect_modified_residue(pdb.path).empty());
}

CG_TEST(pdb_max_molecule_dmax_measures_within_molecules_not_across_them) {
    const ScratchPdb pdb("dmax",
        "ATOM      1 BB   ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
        "ATOM      2 BB   ALA A   2       3.000   0.000   0.000  1.00  0.00\n"
        "ATOM      3 BB   ALA B   1     100.000   0.000   0.000  1.00  0.00\n"
        "ATOM      4 BB   ALA B   2     104.000   0.000   0.000  1.00  0.00\n");

    CG_CHECK(std::abs(pdb_max_molecule_dmax(pdb.path, 2) - 4.0) < 1e-6);
    CG_CHECK(std::abs(pdb_max_molecule_dmax(pdb.path, 4) - 104.0) < 1e-6);
    CG_CHECK(pdb_max_molecule_dmax(pdb.path, 8) == 0.0);
    CG_CHECK(pdb_max_molecule_dmax(pdb.path, 1) == 0.0);
}

CG_TEST(pdb_residue_number_range_reads_a_fragments_own_numbering) {
    const ScratchPdb pdb("resrange",
        "ATOM      1 BB   ALA A 116       0.000   0.000   0.000  1.00  0.00\n"
        "ATOM      2 BB   ALA A 117       3.000   0.000   0.000  1.00  0.00\n"
        "ATOM      3 BB   ALA B 347       6.000   0.000   0.000  1.00  0.00\n"
        "ENDMDL\n"
        "ATOM      4 BB   ALA A 999       9.000   0.000   0.000  1.00  0.00\n");

    const auto [lo, hi] = pdb_residue_number_range(pdb.path);
    CG_CHECK_EQ(lo, 116);
    CG_CHECK_EQ(hi, 347);
}

CG_TEST(pdb_residue_number_range_is_zero_for_a_file_with_no_atoms) {
    const ScratchPdb pdb("resrange_empty", "REMARK nothing here\nEND\n");
    const auto [lo, hi] = pdb_residue_number_range(pdb.path);
    CG_CHECK_EQ(lo, 0);
    CG_CHECK_EQ(hi, 0);
}

CG_TEST(detect_modified_residues_lists_every_distinct_ptm) {
    const ScratchPdb pdb("multi_ptm",
        "ATOM      1  BB  SEP A   1      12.345  1.000  2.000  1.00  0.00\n"
        "ATOM      2  BB  ALA A   2      13.345  1.000  2.000  1.00  0.00\n"
        "ATOM      3  BB  TPO A   3      14.345  1.000  2.000  1.00  0.00\n"
        "ATOM      4  BB  SEP A   4      15.345  1.000  2.000  1.00  0.00\n"
        "ATOM      5  BB  PTR A   5      16.345  1.000  2.000  1.00  0.00\n");

    const auto found = detect_modified_residues(pdb.path);
    CG_CHECK_EQ(static_cast<int>(found.size()), 3);
    CG_CHECK_EQ(found[0], std::string("SEP"));
    CG_CHECK_EQ(found[1], std::string("TPO"));
    CG_CHECK_EQ(found[2], std::string("PTR"));
}

CG_TEST(detect_modified_residues_is_empty_for_an_unmodified_structure) {
    const ScratchPdb pdb("no_ptm", kTwoResiduesOneChain);
    CG_CHECK(detect_modified_residues(pdb.path).empty());
}
