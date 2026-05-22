#include "PlumedWriter.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <sstream>

namespace cg {
namespace {

std::string pair_label(const std::string& prefix, int i, int j) {
    return prefix + "_" + std::to_string(i) + "_" + std::to_string(j);
}

std::vector<std::pair<int, int>> pairs(int n) {
    std::vector<std::pair<int, int>> out;
    out.reserve(static_cast<std::size_t>(n * (n - 1) / 2));
    for (int i = 1; i <= n; ++i)
        for (int j = i + 1; j <= n; ++j) out.emplace_back(i, j);
    return out;
}

std::vector<std::string> labels(const std::string& prefix, int n, bool zero_based = false) {
    std::vector<std::string> out;
    out.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) out.push_back(prefix + std::to_string(i + (zero_based ? 0 : 1)));
    return out;
}

std::string feature_schema(const Config& cfg, const std::vector<std::string>& cols) {
    std::ostringstream out;
    out << "{\n"
        << "  \"system_name\": \"" << cfg.systemName() << "\",\n"
        << "  \"n_prot\": " << cfg.n_prot << ",\n"
        << "  \"atoms_per_prot\": " << cfg.atoms_per_prot << ",\n"
        << "  \"plumed_stride\": " << cfg.plumed_stride << ",\n"
        << "  \"contact_r0_nm\": " << cfg.contact_r0_nm << ",\n"
        << "  \"feature_cols\": [\n";
    for (std::size_t i = 0; i < cols.size(); ++i)
        out << "    \"" << cols[i] << "\"" << (i + 1 == cols.size() ? "" : ",") << '\n';
    out << "  ]\n}\n";
    return out.str();
}

} // namespace

PlumedFiles write_plumed_dat(const Config& cfg,
                             const std::filesystem::path& plumed_path,
                             const std::filesystem::path& colvar_path) {
    const auto pair_ids = pairs(cfg.n_prot);
    const auto group = labels("group", cfg.n_prot, true);
    const auto com = labels("c", cfg.n_prot);
    auto rg = labels("rg", cfg.n_prot);

    std::vector<std::string> dist, contact, feature_cols;
    dist.reserve(pair_ids.size());
    contact.reserve(pair_ids.size());

    std::ostringstream p;
    p << "UNITS LENGTH=nm TIME=ps\n\n";

    for (int i = 0; i < cfg.n_prot; ++i)
        p << group[static_cast<std::size_t>(i)] << ": GROUP NDX_FILE=" << (cfg.resultDir() / "index.ndx").string()
          << " NDX_GROUP=Protein" << i + 1 << '\n';

    p << "\nWHOLEMOLECULES ";
    for (int i = 0; i < cfg.n_prot; ++i)
        p << (i ? " " : "") << "ENTITY" << i << '=' << group[static_cast<std::size_t>(i)];
    p << "\n\n";

    for (int i = 0; i < cfg.n_prot; ++i)
        p << com[static_cast<std::size_t>(i)] << ": COM ATOMS=" << group[static_cast<std::size_t>(i)] << '\n';
    p << '\n';

    for (const auto& [i, j] : pair_ids) {
        const auto d = pair_label("d", i, j);
        dist.push_back(d);
        p << d << ": DISTANCE ATOMS=" << com[static_cast<std::size_t>(i - 1)] << ','
          << com[static_cast<std::size_t>(j - 1)] << '\n';
    }
    p << '\n';

    for (const auto& [i, j] : pair_ids) {
        const auto c = pair_label("cn", i, j);
        contact.push_back(c);
        p << c << ": COORDINATION GROUPA=" << group[static_cast<std::size_t>(i - 1)]
          << " GROUPB=" << group[static_cast<std::size_t>(j - 1)]
          << " SWITCH={RATIONAL R_0=" << cfg.contact_r0_nm << " NN=6 MM=12}\n";
    }
    p << "\ncn_total: COMBINE ARG=" << join(contact, ",")
      << " COEFFICIENTS=" << join(std::vector<std::string>(contact.size(), "1"), ",")
      << " PERIODIC=NO\n\n";

    for (int i = 0; i < cfg.n_prot; ++i)
        p << rg[static_cast<std::size_t>(i)] << ": GYRATION ATOMS=" << group[static_cast<std::size_t>(i)] << '\n';

    p << "\nallgrp: GROUP ATOMS=" << join(group, ",")
      << "\nrg_global: GYRATION ATOMS=allgrp\n\n";

    feature_cols.insert(feature_cols.end(), dist.begin(), dist.end());
    feature_cols.insert(feature_cols.end(), contact.begin(), contact.end());
    feature_cols.insert(feature_cols.end(), rg.begin(), rg.end());
    feature_cols.insert(feature_cols.end(), {"rg_global", "cn_total"});

    p << "PRINT STRIDE=" << cfg.plumed_stride << " FILE=" << colvar_path.string()
      << " ARG=" << join(feature_cols, ",") << '\n';

    const auto schema_path = cfg.resultDir() / "feature_schema.json";
    write_text(plumed_path, p.str());
    write_text(schema_path, feature_schema(cfg, feature_cols));
    return {plumed_path, schema_path, feature_cols};
}


PlumedFiles write_metad_plumed_dat(const Config& cfg,
                                   const std::filesystem::path& plumed_path,
                                   const std::filesystem::path& ndx_path,
                                   const std::filesystem::path& model_path,
                                   const std::filesystem::path& colvar_path) {
    const auto pair_ids = pairs(cfg.n_prot);
    const auto group = labels("group", cfg.n_prot, true);
    const auto com = labels("c", cfg.n_prot);

    std::vector<std::string> dist, contact, args, nodes;
    dist.reserve(pair_ids.size());
    contact.reserve(pair_ids.size());

    std::ostringstream p;
    p << "UNITS LENGTH=nm TIME=ps\n\n";

    for (int i = 0; i < cfg.n_prot; ++i)
        p << group[static_cast<std::size_t>(i)] << ": GROUP NDX_FILE=" << ndx_path.string()
          << " NDX_GROUP=Protein" << i + 1 << '\n';

    p << '\n';
    for (int i = 0; i < cfg.n_prot; ++i)
        p << com[static_cast<std::size_t>(i)] << ": COM ATOMS=" << group[static_cast<std::size_t>(i)] << '\n';

    p << '\n';
    for (const auto& [i, j] : pair_ids) {
        const auto d = "d" + std::to_string(i) + std::to_string(j);
        dist.push_back(d);
        p << d << ": DISTANCE ATOMS=" << com[static_cast<std::size_t>(i - 1)] << ','
          << com[static_cast<std::size_t>(j - 1)] << '\n';
    }

    p << '\n';
    for (const auto& [i, j] : pair_ids) {
        const auto c = "cn" + std::to_string(i) + std::to_string(j);
        contact.push_back(c);
        p << c << ": COORDINATION GROUPA=" << group[static_cast<std::size_t>(i - 1)]
          << " GROUPB=" << group[static_cast<std::size_t>(j - 1)]
          << " SWITCH={RATIONAL R_0=" << cfg.contact_r0_nm << "}\n";
    }

    p << "\ncn_total: COMBINE ARG=" << join(contact, ",")
      << " COEFFICIENTS=" << join(std::vector<std::string>(contact.size(), "1"), ",")
      << " PERIODIC=NO\n\n";

    p << "allgrp: GROUP ATOMS=" << join(group, ",") << "\n"
      << "rg: GYRATION ATOMS=allgrp\n\n";

    args = dist;
    args.insert(args.end(), {"rg", "cn_total"});
    p << "mycv: PYTORCH_MODEL FILE=" << model_path.string()
      << " ARG=" << join(args, ",") << "\n\n";

    for (int i = 0; i < cfg.metad_nodes; ++i) nodes.push_back("mycv.node-" + std::to_string(i));

    p << "metad: METAD ...\n"
      << "ARG=" << join(nodes, ",") << "\n"
      << "PACE=" << cfg.metad_pace << "\n"
      << "HEIGHT=" << cfg.metad_height << "\n"
      << "SIGMA=" << cfg.metad_sigma << "\n"
      << "BIASFACTOR=" << cfg.metad_biasfactor << "\n"
      << "TEMP=" << cfg.temperature_K << "\n"
      << "FILE=HILLS\n";
    if (cfg.metad_walkers > 1) p << "WALKERS_MPI\n";
    p << "... metad:\n\n";

    auto print_args = nodes;
    print_args.push_back("metad.bias");
    p << "PRINT STRIDE=" << cfg.metad_print_stride
      << " ARG=" << join(print_args, ",")
      << " FILE=" << colvar_path.string() << '\n';

    write_text(plumed_path, p.str());
    return {plumed_path, {}, args};
}

PlumedFiles write_plumed_dat(const Config& cfg) {
    return write_plumed_dat(cfg, cfg.resultDir() / "plumed.dat", cfg.resultDir() / "COLVAR");
}

} // namespace cg
