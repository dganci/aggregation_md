#include "RunDiagnostics.hpp"
#include "ClusterAnalysis.hpp"
#include "ColvarTable.hpp"
#include "GromacsDriver.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cg {
namespace {

std::string pair_contact_name(int i, int j) { return "cn_" + std::to_string(i) + "_" + std::to_string(j); }

} // namespace

std::string EnergyDiagnostics::to_json() const {
    std::ostringstream o;
    o << "{\"temperature_K\":" << temperature_K.to_json()
      << ",\"total_energy_kJ_per_mol\":" << total_energy_kJ_per_mol.to_json();
    if (has_pressure) o << ",\"pressure_bar\":" << pressure_bar.to_json();
    if (has_density) o << ",\"density_kg_per_m3\":" << density_kg_per_m3.to_json();
    o << "}";
    return o.str();
}

std::string StructuralDiagnostics::to_json() const {
    std::ostringstream o;
    o << std::setprecision(10)
      << "{\"n_frames\":" << n_frames
      << ",\"total_time_ps\":" << total_time_ps
      << ",\"per_chain_rg_nm\":[";
    for (std::size_t i = 0; i < per_chain_rg_nm.size(); ++i) {
        if (i) o << ',';
        o << per_chain_rg_nm[i].to_json();
    }
    o << "],\"global_rg_nm\":" << global_rg_nm.to_json()
      << ",\"total_contacts\":" << total_contacts.to_json()
      << ",\"largest_cluster\":" << largest_cluster.to_json()
      << ",\"n_clusters\":" << n_clusters.to_json()
      << ",\"cluster_size_histogram\":{";
    for (std::size_t i = 0; i < cluster_size_histogram.size(); ++i) {
        if (i) o << ',';
        o << '"' << (i + 1) << "\":" << cluster_size_histogram[i];
    }
    o << "},\"largest_cluster_histogram\":{";
    for (std::size_t i = 0; i < largest_cluster_histogram.size(); ++i) {
        if (i) o << ',';
        o << '"' << (i + 1) << "\":" << largest_cluster_histogram[i];
    }
    o << "},\"pair_contact_occupancy\":{";
    for (std::size_t i = 0; i < pair_contact_occupancy.size() && i < pair_labels.size(); ++i) {
        if (i) o << ',';
        o << '"' << pair_labels[i] << "\":" << pair_contact_occupancy[i];
    }
    o << "}}";
    return o.str();
}

std::string RunDiagnostics::to_json(const std::string& segment_label) const {
    std::ostringstream o;
    o << "{\"segment\":\"" << segment_label << "\",\"energy\":" << energy.to_json();
    if (has_structure) o << ",\"structure\":" << structure.to_json();
    o << "}";
    return o.str();
}

EnergyDiagnostics compute_energy_diagnostics(const GromacsDriver& gmx,
                                             const std::filesystem::path& edr_path,
                                             const std::filesystem::path& scratch_dir,
                                             bool include_pressure,
                                             bool include_density) {
    std::filesystem::create_directories(scratch_dir);

    EnergyDiagnostics diag;
    diag.has_pressure = include_pressure;
    diag.has_density = include_density;

    std::vector<std::string> observables = {"Temperature", "Total-Energy"};
    if (include_pressure) observables.push_back("Pressure");
    if (include_density) observables.push_back("Density");

    const auto xvg = scratch_dir / "diag_energy.xvg";
    gmx.energy(edr_path, xvg, observables);

    const auto series = read_xvg_by_legend(xvg);
    const auto take = [&](const std::string& name) {
        const auto it = series.find(normalise_legend(name));
        if (it == series.end())
            throw std::runtime_error("gmx energy produced no '" + name + "' column in " + xvg.string() +
                                     " - the term is absent from the .edr for this stage.");
        return compute_time_series_stats(it->second);
    };

    diag.temperature_K = take("Temperature");
    diag.total_energy_kJ_per_mol = take("Total-Energy");
    if (include_pressure) diag.pressure_bar = take("Pressure");
    if (include_density) diag.density_kg_per_m3 = take("Density");

    return diag;
}

StructuralDiagnostics compute_structural_diagnostics(const std::vector<std::filesystem::path>& colvar_paths,
                                                     int n_prot,
                                                     double dt_colvar_ps,
                                                     double pair_contact_threshold) {
    const auto table = read_colvars(colvar_paths);

    StructuralDiagnostics diag;
    diag.n_frames = static_cast<int>(table.rows.size());
    diag.total_time_ps = table.total_time_ps(dt_colvar_ps);

    diag.per_chain_rg_nm.reserve(static_cast<std::size_t>(n_prot));
    for (int i = 1; i <= n_prot; ++i)
        diag.per_chain_rg_nm.push_back(
            compute_time_series_stats(column_series(table, "rg" + std::to_string(i), dt_colvar_ps)));

    for (const char* name : {"rg_com", "rg_global"}) {
        if (!table.col.count(name)) continue;
        diag.global_rg_nm = compute_time_series_stats(column_series(table, name, dt_colvar_ps));
        break;
    }
    if (table.col.count("cn_total"))
        diag.total_contacts = compute_time_series_stats(column_series(table, "cn_total", dt_colvar_ps));

    std::vector<std::pair<int, int>> pairs;
    std::vector<std::size_t> cols;
    for (int i = 1; i <= n_prot; ++i) {
        for (int j = i + 1; j <= n_prot; ++j) {
            const auto name = pair_contact_name(i, j);
            const auto it = table.col.find(name);
            if (it == table.col.end()) throw std::runtime_error("Missing expected pair contact column: " + name);
            pairs.emplace_back(i - 1, j - 1);
            cols.push_back(it->second);
        }
    }

    diag.cluster_size_histogram.assign(static_cast<std::size_t>(n_prot), 0);
    diag.largest_cluster_histogram.assign(static_cast<std::size_t>(n_prot), 0);
    diag.pair_contact_occupancy.assign(cols.size(), 0.0);

    std::vector<std::pair<double, double>> lcc_series;
    std::vector<std::pair<double, double>> ncluster_series;
    lcc_series.reserve(table.rows.size());
    ncluster_series.reserve(table.rows.size());

    for (std::size_t frame = 0; frame < table.rows.size(); ++frame) {
        ContactGraph adj(static_cast<std::size_t>(n_prot));
        for (std::size_t k = 0; k < cols.size(); ++k) {
            if (table.rows[frame][cols[k]] < pair_contact_threshold) continue;
            diag.pair_contact_occupancy[k] += 1.0;
            const auto [a, b] = pairs[k];
            adj[static_cast<std::size_t>(a)].push_back(b);
            adj[static_cast<std::size_t>(b)].push_back(a);
        }

        const auto sizes = connected_component_sizes(n_prot, adj);
        for (int size : sizes) ++diag.cluster_size_histogram[static_cast<std::size_t>(size - 1)];
        if (!sizes.empty()) ++diag.largest_cluster_histogram[static_cast<std::size_t>(sizes.front() - 1)];

        const double t = table.has_time ? table.rows[frame][0]
                                        : static_cast<double>(frame) * dt_colvar_ps;
        lcc_series.emplace_back(t, sizes.empty() ? 0.0 : static_cast<double>(sizes.front()));
        ncluster_series.emplace_back(t, static_cast<double>(sizes.size()));
    }

    if (!table.rows.empty())
        for (auto& occupancy : diag.pair_contact_occupancy)
            occupancy /= static_cast<double>(table.rows.size());

    diag.largest_cluster = compute_time_series_stats(lcc_series);
    diag.n_clusters = compute_time_series_stats(ncluster_series);
    diag.pair_labels.reserve(cols.size());
    for (const auto& [a, b] : pairs) diag.pair_labels.push_back(pair_contact_name(a + 1, b + 1));

    return diag;
}

void append_diagnostics_jsonl(const std::filesystem::path& path,
                              const RunDiagnostics& diag,
                              const std::string& segment_label) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::app);
    if (!out) throw std::runtime_error("Cannot write diagnostics JSONL: " + path.string());
    out << diag.to_json(segment_label) << '\n';
}

} // namespace cg
