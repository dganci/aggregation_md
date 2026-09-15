#pragma once

#include <filesystem>
#include <map>
#include <string>

namespace cg {

/// A tiny key=value record of the parameters a chunked/batched run was
/// *started* with, written next to its output.
class RunLedger {
public:
    /// Reads `path` if it exists; an absent file yields an empty ledger.
    static RunLedger load(const std::filesystem::path& path);

    /// True when the ledger held no entries (i.e. this is a fresh run).
    bool empty() const;

    /// Records `key` on a fresh ledger, or - if it is already present with a
    /// different value - throws explaining which option must be restored.
    void require(const std::string& key, const std::string& value, const std::string& option);

    /// Unconditionally sets `key` (for bookkeeping entries that legitimately
    /// change, as opposed to the invariants require() protects).
    void set(const std::string& key, const std::string& value);

    /// Value of `key`, or an empty string if absent.
    std::string get(const std::string& key) const;

    void save(const std::filesystem::path& path) const;

private:
    std::map<std::string, std::string> values_;
};

} // namespace cg
