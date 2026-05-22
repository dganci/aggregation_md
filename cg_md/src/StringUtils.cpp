#include "StringUtils.hpp"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>

namespace cg {

std::string trim(const std::string& s) {
    const auto b = std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); });
    const auto e = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c); }).base();
    return b < e ? std::string(b, e) : std::string{};
}

bool starts_with(const std::string& s, const std::string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

bool ends_with(const std::string& s, const std::string& suffix) {
    return suffix.size() <= s.size() && std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

std::vector<std::string> split_ws(const std::string& s) {
    std::istringstream in(s);
    std::vector<std::string> out;
    for (std::string x; in >> x;) out.push_back(std::move(x));
    return out;
}

std::string join(const std::vector<std::string>& xs, const std::string& sep) {
    std::ostringstream out;
    for (std::size_t i = 0; i < xs.size(); ++i) out << (i ? sep : "") << xs[i];
    return out.str();
}

std::string shell_quote(const std::string& s) {
    std::string out{"'"};
    for (char c : s) out += (c == '\'') ? "'\\''" : std::string(1, c);
    return out + "'";
}

std::string basename_without_ext(const std::string& path) {
    const auto sep = path.find_last_of("/\\");
    const auto name = sep == std::string::npos ? path : path.substr(sep + 1);
    const auto dot = name.find_last_of('.');
    return dot == std::string::npos ? name : name.substr(0, dot);
}

std::string to_string_fixed(double x, int precision) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << x;
    return out.str();
}

std::string zero_padded(int value, int width) {
    std::ostringstream out;
    out << std::setw(width) << std::setfill('0') << value;
    return out.str();
}

} // namespace cg
