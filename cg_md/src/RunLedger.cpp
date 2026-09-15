#include "RunLedger.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <sstream>
#include <stdexcept>

namespace cg {

RunLedger RunLedger::load(const std::filesystem::path& path) {
    RunLedger ledger;
    if (!std::filesystem::exists(path)) return ledger;

    for (const auto& line : read_lines(path)) {
        const auto trimmed = trim(line);
        if (trimmed.empty() || starts_with(trimmed, "#")) continue;
        const auto eq = trimmed.find('=');
        if (eq == std::string::npos) continue;
        ledger.values_[trim(trimmed.substr(0, eq))] = trim(trimmed.substr(eq + 1));
    }
    return ledger;
}

bool RunLedger::empty() const { return values_.empty(); }

void RunLedger::require(const std::string& key, const std::string& value, const std::string& option) {
    const auto it = values_.find(key);
    if (it == values_.end()) {
        values_[key] = value;
        return;
    }
    if (it->second == value) return;

    throw std::runtime_error(
        "This run was started with " + key + " = " + it->second + " (" + option +
        "), but the current invocation gives " + value +
        ". Resuming would continue from state written under the original setting while interpreting it "
        "under the new one - which silently corrupts the simulated-time axis, the walker layout, or the "
        "descriptor set, with no error anywhere downstream. Restore " + option +
        " to its original value to continue this run, or start a fresh run in a new --project-dir "
        "(the recorded values are in this stage's *.ledger file).");
}

void RunLedger::set(const std::string& key, const std::string& value) { values_[key] = value; }

std::string RunLedger::get(const std::string& key) const {
    const auto it = values_.find(key);
    return it == values_.end() ? std::string{} : it->second;
}

void RunLedger::save(const std::filesystem::path& path) const {
    std::ostringstream out;
    out << "# Segmentation parameters this run was started with; see RunLedger.hpp.\n"
        << "# Changing any of these mid-run is refused, not silently accepted.\n";
    for (const auto& [key, value] : values_) out << key << " = " << value << '\n';
    write_text(path, out.str());
}

} // namespace cg
