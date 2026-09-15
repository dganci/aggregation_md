#include "test_framework.hpp"
#include "RunDiagnostics.hpp"
#include "FileUtils.hpp"

#include <cmath>
#include <exception>
#include <filesystem>

using namespace cg;

namespace {
std::filesystem::path temp_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / ("cg_md_test_diag_" + name);
}
bool close(double a, double b, double eps = 1e-9) { return std::abs(a - b) < eps; }
} // namespace

CG_TEST(compute_structural_diagnostics_reduces_per_chain_rg_and_cluster_histogram) {
    const auto path = temp_path("colvar.dat");
    write_lines(path, {
        "#! FIELDS rg1 rg2 cn_1_2",
        "1.0 1.0 0.0",
        "2.0 2.0 1.0"
    });

    const auto diag = compute_structural_diagnostics({path}, /*n_prot=*/2, /*dt_colvar_ps=*/1000.0,
                                                      /*pair_contact_threshold=*/1.0);

    CG_CHECK_EQ(diag.n_frames, 2);
    CG_CHECK_EQ(static_cast<int>(diag.per_chain_rg_nm.size()), 2);
    CG_CHECK(close(diag.per_chain_rg_nm[0].mean, 1.5));
    CG_CHECK(close(diag.per_chain_rg_nm[1].mean, 1.5));

    CG_CHECK_EQ(static_cast<int>(diag.cluster_size_histogram.size()), 2);
    CG_CHECK_EQ(diag.cluster_size_histogram[0], 2);
    CG_CHECK_EQ(diag.cluster_size_histogram[1], 1);

    std::filesystem::remove(path);
}

CG_TEST(compute_structural_diagnostics_throws_on_missing_rg_column) {
    const auto path = temp_path("colvar_missing_rg.dat");
    write_lines(path, {"#! FIELDS cn_1_2", "0.0"});

    bool threw = false;
    try {
        compute_structural_diagnostics({path}, /*n_prot=*/2, 1000.0, 1.0);
    } catch (const std::exception&) {
        threw = true;
    }
    CG_CHECK(threw);

    std::filesystem::remove(path);
}

CG_TEST(energy_diagnostics_to_json_omits_pressure_and_density_when_absent) {
    EnergyDiagnostics diag;
    diag.has_pressure = false;
    diag.has_density = false;
    const auto json = diag.to_json();
    CG_CHECK(json.find("temperature_K") != std::string::npos);
    CG_CHECK(json.find("total_energy_kJ_per_mol") != std::string::npos);
    CG_CHECK(json.find("pressure_bar") == std::string::npos);
    CG_CHECK(json.find("density_kg_per_m3") == std::string::npos);
}

CG_TEST(energy_diagnostics_to_json_includes_pressure_and_density_when_present) {
    EnergyDiagnostics diag;
    diag.has_pressure = true;
    diag.has_density = true;
    const auto json = diag.to_json();
    CG_CHECK(json.find("pressure_bar") != std::string::npos);
    CG_CHECK(json.find("density_kg_per_m3") != std::string::npos);
}

CG_TEST(run_diagnostics_to_json_omits_structure_when_not_available) {
    RunDiagnostics diag;
    diag.has_structure = false;
    const auto json = diag.to_json("nvt");
    CG_CHECK(json.find("\"segment\":\"nvt\"") != std::string::npos);
    CG_CHECK(json.find("\"structure\"") == std::string::npos);
}
