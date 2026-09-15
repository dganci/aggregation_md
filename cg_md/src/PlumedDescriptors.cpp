#include "PlumedWriter.hpp"
#include "Config.hpp"
#include "StringUtils.hpp"

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

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

} // namespace

PlumedDescriptors build_descriptors(const Config& cfg, const std::filesystem::path& ndx_path) {
    PlumedDescriptors d;

    const auto pair_ids = pairs(cfg.n_prot);
    const auto group = labels("group", cfg.n_prot, true);
    const auto com = labels("c", cfg.n_prot);
    d.per_chain_rg = labels("rg", cfg.n_prot);
    d.distances.reserve(pair_ids.size());
    d.contacts.reserve(pair_ids.size());

    std::ostringstream p;

    for (int i = 0; i < cfg.n_prot; ++i)
        p << group[static_cast<std::size_t>(i)] << ": GROUP NDX_FILE=" << ndx_path.string()
          << " NDX_GROUP=Protein" << i + 1 << '\n';

    const int merged_chains = cfg.merge_chains.empty()
        ? 1
        : 1 + static_cast<int>(std::count(cfg.merge_chains.begin(), cfg.merge_chains.end(), ','));

    p << '\n';
    if (merged_chains == 2 && cfg.atoms_per_prot > 1 && cfg.atoms_per_prot % 2 == 0) {
        const int per_chain = cfg.atoms_per_prot / 2;
        p << "# Atom order, not the plain group: chain A forward then chain B backward, so\n"
             "# that no step in the walk crosses the rod. See PlumedWriter.cpp.\n"
             "WHOLEMOLECULES";
        for (int i = 0; i < cfg.n_prot; ++i) {
            const int base = i * cfg.atoms_per_prot;
            p << " ENTITY" << i << '=' << (base + 1) << '-' << (base + per_chain);
            for (int k = cfg.atoms_per_prot; k > per_chain; --k) p << ',' << (base + k);
        }
        p << '\n';
    } else {
        if (merged_chains > 2)
            p << "# WARNING: " << merged_chains << " merged chains. Only the two-chain case has a\n"
                 "# verified atom order; this falls back to the plain group, which is correct only\n"
                 "# while every A->B junction stays inside half the minimum image distance.\n";
        p << "WHOLEMOLECULES ";
        for (int i = 0; i < cfg.n_prot; ++i)
            p << (i ? " " : "") << "ENTITY" << i << '=' << group[static_cast<std::size_t>(i)];
        p << '\n';
    }
    p << '\n';

    for (int i = 0; i < cfg.n_prot; ++i)
        p << com[static_cast<std::size_t>(i)] << ": COM ATOMS=" << group[static_cast<std::size_t>(i)] << '\n';
    p << '\n';

    for (const auto& [i, j] : pair_ids) {
        const auto label = pair_label("d", i, j);
        d.distances.push_back(label);
        p << label << ": DISTANCE ATOMS=" << com[static_cast<std::size_t>(i - 1)] << ','
          << com[static_cast<std::size_t>(j - 1)] << '\n';
    }
    p << '\n';

    for (const auto& [i, j] : pair_ids) {
        const auto label = pair_label("cn", i, j);
        d.contacts.push_back(label);
        p << label << ": COORDINATION GROUPA=" << group[static_cast<std::size_t>(i - 1)]
          << " GROUPB=" << group[static_cast<std::size_t>(j - 1)]
          << " SWITCH={RATIONAL R_0=" << cfg.contact_r0_nm
          << " NN=" << cfg.contact_nn << " MM=" << cfg.contact_mm << "}\n";
    }

    p << "\ncn_total: COMBINE ARG=" << join(d.contacts, ",")
      << " COEFFICIENTS=" << repeat_csv("1", d.contacts.size())
      << " PERIODIC=NO\n\n";

    for (int i = 0; i < cfg.n_prot; ++i)
        p << d.per_chain_rg[static_cast<std::size_t>(i)] << ": GYRATION ATOMS="
          << group[static_cast<std::size_t>(i)] << '\n';

    const double rg_coeff = 1.0 / (static_cast<double>(cfg.n_prot) * cfg.n_prot);
    p << "\nrg_com_sq: COMBINE ARG=" << join(d.distances, ",")
      << " POWERS=" << repeat_csv("2", d.distances.size())
      << " COEFFICIENTS=" << repeat_csv(to_string_fixed(rg_coeff, 8), d.distances.size())
      << " PERIODIC=NO\n"
      << "rg_com: CUSTOM ARG=rg_com_sq FUNC=sqrt(x) PERIODIC=NO\n";

    d.text = p.str();

    d.feature_cols.insert(d.feature_cols.end(), d.distances.begin(), d.distances.end());
    d.feature_cols.insert(d.feature_cols.end(), d.contacts.begin(), d.contacts.end());
    d.feature_cols.insert(d.feature_cols.end(), d.per_chain_rg.begin(), d.per_chain_rg.end());
    d.feature_cols.insert(d.feature_cols.end(), {"rg_com", "cn_total"});
    return d;
}

SortedDescriptors build_sorted_descriptors(const PlumedDescriptors& d) {
    SortedDescriptors s;
    std::ostringstream p;

    const auto emit_sort = [&](const std::string& label, const std::vector<std::string>& args) {
        if (args.empty()) return;
        p << label << ": SORT ARG=" << join(args, ",") << '\n';
        for (std::size_t i = 1; i <= args.size(); ++i)
            s.feature_cols.push_back(label + "." + std::to_string(i));
    };

    p << "# Permutation invariance: identical protomers are interchangeable, so the\n"
      << "# same physical state can appear with any labelling of the chains. Sorting\n"
      << "# each permutable block removes that degeneracy, so one configuration maps\n"
      << "# to one CV value. Must match cg_cvgen's offline sort exactly.\n";
    emit_sort("sorted_d", d.distances);
    emit_sort("sorted_cn", d.contacts);
    emit_sort("sorted_rg", d.per_chain_rg);

    s.feature_cols.push_back("rg_com");
    s.feature_cols.push_back("cn_total");

    s.text = p.str();
    return s;
}

} // namespace cg
