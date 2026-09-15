#include "test_framework.hpp"
#include "GroUtils.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <cmath>
#include <filesystem>

using namespace cg;

namespace {
std::filesystem::path temp_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_" + name);
}
} // namespace

CG_TEST(normalize_gro_ion_names_rewrites_na_cl_to_martini_spelling) {
    const auto path = temp_path("ions.gro");
    write_lines(path, {
        "test system",
        "2",
        "    1NA      NA    1   1.000   1.000   1.000",
        "    2CL      CL    2   2.000   2.000   2.000",
        "   5.00000   5.00000   5.00000"
    });

    const auto stats = normalize_gro_ion_names(path);
    CG_CHECK_EQ(stats.atom_names, std::size_t{2});
    CG_CHECK_EQ(stats.residue_names, std::size_t{2});
    CG_CHECK_EQ(stats.total(), std::size_t{4});

    const auto lines = read_lines(path);
    CG_CHECK_EQ(trim(lines[2].substr(5, 5)), "NA+");
    CG_CHECK_EQ(trim(lines[2].substr(10, 5)), "NA+");
    CG_CHECK_EQ(trim(lines[3].substr(5, 5)), "CL-");
    CG_CHECK_EQ(trim(lines[3].substr(10, 5)), "CL-");

    std::filesystem::remove(path);
}

CG_TEST(normalize_gro_ion_names_is_idempotent_when_already_martini_spelled) {
    const auto path = temp_path("ions_already.gro");
    write_lines(path, {
        "test system",
        "1",
        "    1NA+    NA+    1   1.000   1.000   1.000",
        "   5.00000   5.00000   5.00000"
    });

    const auto stats = normalize_gro_ion_names(path);
    CG_CHECK_EQ(stats.total(), std::size_t{0});

    std::filesystem::remove(path);
}

CG_TEST(normalize_gro_ion_names_leaves_non_ion_residues_untouched) {
    const auto path = temp_path("protein.gro");
    write_lines(path, {
        "test system",
        "1",
        "    1SER     BB    1   1.000   1.000   1.000",
        "   5.00000   5.00000   5.00000"
    });

    const auto stats = normalize_gro_ion_names(path);
    CG_CHECK_EQ(stats.total(), std::size_t{0});

    std::filesystem::remove(path);
}

CG_TEST(read_gro_box_vectors_reads_a_rectangular_three_field_box) {
    const auto path = temp_path("box_rect.gro");
    write_lines(path, {"x", "1", "    1W        W    1   1.000   1.000   1.000",
                       "   4.00000   5.00000   6.00000"});

    const auto box = read_gro_box_vectors(path);
    CG_CHECK(box[0][0] == 4.0 && box[1][1] == 5.0 && box[2][2] == 6.0);
    CG_CHECK(box[0][1] == 0.0 && box[2][0] == 0.0 && box[2][1] == 0.0);

    std::filesystem::remove(path);
}

CG_TEST(min_image_distance_of_a_cube_is_its_edge) {
    const auto path = temp_path("box_cube.gro");
    write_lines(path, {"x", "1", "    1W        W    1   1.000   1.000   1.000",
                       "   7.00000   7.00000   7.00000"});

    CG_CHECK(std::abs(min_image_distance_nm(read_gro_box_vectors(path)) - 7.0) < 1e-9);

    std::filesystem::remove(path);
}

CG_TEST(min_image_distance_of_a_rhombic_dodecahedron_is_a_not_a_over_root_two) {
    const auto path = temp_path("box_dodec.gro");
    write_lines(path, {"x", "1", "    1W        W    1   1.000   1.000   1.000",
                       "  18.37488  18.37488  12.99300   0.00000   0.00000   0.00000"
                       "   0.00000   9.18744   9.18744"});

    const auto box = read_gro_box_vectors(path);
    CG_CHECK(std::abs(box[2][0] - 9.18744) < 1e-6);
    CG_CHECK(std::abs(box[2][1] - 9.18744) < 1e-6);

    const double image = min_image_distance_nm(box);
    CG_CHECK(std::abs(image - 18.37488) < 1e-4);
    CG_CHECK(image > 12.993 + 1.0);

    std::filesystem::remove(path);
}

CG_TEST(periodic_margin_needs_two_cutoffs_of_headroom) {
    CG_CHECK(periodic_margin_ok(8.0, 10.3, 1.1));
    CG_CHECK(periodic_margin_ok(8.0, 10.2, 1.1));
    CG_CHECK(!periodic_margin_ok(8.0, 10.1, 1.1));
    CG_CHECK(!periodic_margin_ok(8.0, 9.0, 1.1));

    CG_CHECK(periodic_margin_ok(27.1, 37.9, 1.1));
    CG_CHECK(!periodic_margin_ok(27.1, 28.6, 1.1));
}
