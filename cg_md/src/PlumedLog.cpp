#include "PlumedLog.hpp"
#include "StringUtils.hpp"

#include <cctype>
#include <sstream>
#include <string>
#include <vector>

namespace cg {
namespace {

std::string lowercased(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (const unsigned char c : s) out += static_cast<char>(std::tolower(c));
    return out;
}

int as_count(const std::string& token) {
    std::string digits;
    for (const char c : token) {
        if (std::isdigit(static_cast<unsigned char>(c))) digits += c;
        else if (!digits.empty()) break;
        else if (c != '+') return -1;
    }
    if (digits.empty()) return -1;
    try {
        return std::stoi(digits);
    } catch (const std::exception&) {
        return -1;
    }
}

} // namespace

int plumed_walker_count(const std::string& log_text) {
    std::istringstream in(log_text);
    for (std::string line; std::getline(in, line);) {
        const auto tokens = split_ws(lowercased(line));

        for (std::size_t i = 0; i + 1 < tokens.size(); ++i) {
            if (tokens[i] == "walkers:" || (tokens[i] == "walkers" && tokens[i + 1] == ":")) {
                const int n = as_count(tokens[i] == "walkers:" ? tokens[i + 1]
                                                              : (i + 2 < tokens.size() ? tokens[i + 2] : ""));
                if (n > 0) return n;
            }
            if (tokens[i] == "multiple" && tokens[i + 1] == "walkers" && i > 0) {
                const int n = as_count(tokens[i - 1]);
                if (n > 0) return n;
            }
        }
    }
    return -1;
}

bool is_plumed_backup(const std::string& filename, const std::string& target) {
    static const std::string kPrefix = "bck.";
    if (!starts_with(filename, kPrefix)) return false;
    const auto suffix = "." + target;
    if (!ends_with(filename, suffix)) return false;
    const auto middle = filename.substr(kPrefix.size(),
                                        filename.size() - kPrefix.size() - suffix.size());
    if (middle.empty()) return false;
    if (middle == "last") return true;
    return middle.find_first_not_of("0123456789") == std::string::npos;
}

} // namespace cg
