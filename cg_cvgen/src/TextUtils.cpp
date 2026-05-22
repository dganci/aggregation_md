#include "TextUtils.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace cgcv {

std::string shell_quote(const std::string& value) {
    std::string out = "'";
    for (const char c : value) out += c == '\'' ? "'\\''" : std::string(1, c);
    out += "'";
    return out;
}

std::string json_escape(const std::string& value) {
    std::string out;
    out.reserve(value.size() + 8);
    for (const char c : value) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default: out += c; break;
        }
    }
    return out;
}

std::string join(const std::vector<std::string>& values, const std::string& sep) {
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) out << (i ? sep : "") << values[i];
    return out.str();
}

namespace {
std::vector<std::string> split_csv(std::string s) {
    std::replace(s.begin(), s.end(), ';', ',');
    std::stringstream ss(s);
    std::vector<std::string> out;
    for (std::string item; std::getline(ss, item, ',');) {
        item.erase(item.begin(), std::find_if(item.begin(), item.end(), [](unsigned char ch){ return !std::isspace(ch); }));
        item.erase(std::find_if(item.rbegin(), item.rend(), [](unsigned char ch){ return !std::isspace(ch); }).base(), item.end());
        if (!item.empty()) out.push_back(item);
    }
    return out;
}
}

std::vector<int> parse_int_list(const std::string& value) {
    std::vector<int> out;
    for (const auto& item : split_csv(value)) {
        int x = 0;
        try { x = std::stoi(item); }
        catch (...) { throw std::runtime_error("Invalid integer list: " + value); }
        if (x <= 0) throw std::runtime_error("Integer list values must be > 0: " + value);
        out.push_back(x);
    }
    if (out.empty()) throw std::runtime_error("Empty integer list: " + value);
    return out;
}

std::vector<int> parse_layers(const std::string& value) { return parse_int_list(value); }

std::string to_json_array(const std::vector<int>& values) {
    std::ostringstream out;
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) out << (i ? ", " : "") << values[i];
    out << ']';
    return out.str();
}

std::filesystem::path default_backend_path(const char* argv0) {
    const auto exe = std::filesystem::absolute(argv0 ? argv0 : "cg_cvgen");
    const auto bin_dir = exe.parent_path();
    const std::vector<std::filesystem::path> candidates = {
        std::filesystem::current_path() / "scripts" / "cvgen_backend.py",
        bin_dir / "scripts" / "cvgen_backend.py",
        bin_dir.parent_path() / "scripts" / "cvgen_backend.py",
        bin_dir.parent_path().parent_path() / "scripts" / "cvgen_backend.py"
    };
    for (const auto& p : candidates) if (std::filesystem::exists(p)) return p;
    return candidates.front();
}

} // namespace cgcv
