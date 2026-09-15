#include "PlumedWriter.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
namespace cg {
namespace {

std::string feature_schema(const Config& cfg, const std::vector<std::string>& cols) {
    std::ostringstream out;
    out << "{\n"
        << "  \"system_name\": \"" << cfg.systemName() << "\",\n"
        << "  \"n_prot\": " << cfg.n_prot << ",\n"
        << "  \"atoms_per_prot\": " << cfg.atoms_per_prot << ",\n"
        << "  \"plumed_stride\": " << cfg.plumed_stride << ",\n"
        << "  \"contact_r0_nm\": " << cfg.contact_r0_nm << ",\n"
        << "  \"contact_nn\": " << cfg.contact_nn << ",\n"
        << "  \"contact_mm\": " << cfg.contact_mm << ",\n"
        << "  \"n_features\": " << cols.size() << ",\n"
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
    const auto d = build_descriptors(cfg, cfg.resultDir() / "index.ndx");

    std::ostringstream p;
    p << "# Unbiased production descriptors for " << cfg.systemName() << ".\n"
      << "# Column order here IS the training feature order consumed by cg_cvgen\n"
      << "# and re-used verbatim as the PYTORCH_MODEL ARG list under metadynamics.\n"
      << "UNITS LENGTH=nm TIME=ps\n\n"
      << d.text << '\n'
      << "PRINT STRIDE=" << cfg.plumed_stride << " FILE=" << colvar_path.string()
      << " ARG=" << join(d.feature_cols, ",") << '\n';

    const auto schema_path = cfg.resultDir() / "feature_schema.json";
    if (!cfg.dry_run) {
        write_text(plumed_path, p.str());
        write_text(schema_path, feature_schema(cfg, d.feature_cols));
    }
    return {plumed_path, schema_path, d.feature_cols};
}

PlumedFiles write_plumed_dat(const Config& cfg) {
    return write_plumed_dat(cfg, cfg.resultDir() / "plumed.dat", cfg.resultDir() / "COLVAR");
}

} // namespace cg
