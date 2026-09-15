#include "FileUtils.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace cg {

std::string read_text(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read " + path.string());
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

std::string read_tail(const std::filesystem::path& path, std::size_t max_bytes) {
    std::error_code ec;
    const auto size = std::filesystem::file_size(path, ec);
    if (ec) return {};

    std::ifstream in(path, std::ios::binary);
    if (!in) return {};

    const auto want = static_cast<std::streamoff>(std::min<std::uintmax_t>(size, max_bytes));
    in.seekg(-want, std::ios::end);
    if (!in) in.seekg(0, std::ios::beg);

    std::string out(static_cast<std::size_t>(want), '\0');
    in.read(out.data(), want);
    out.resize(static_cast<std::size_t>(in.gcount()));
    return out;
}

std::vector<std::string> read_lines(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read " + path.string());
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) lines.push_back(std::move(line));
    return lines;
}

void write_text(const std::filesystem::path& path, const std::string& text) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write " + path.string());
    out << text;
}

void write_lines(const std::filesystem::path& path, const std::vector<std::string>& lines) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write " + path.string());
    for (const auto& line : lines) out << line << '\n';
}

void copy_overwrite(const std::filesystem::path& from, const std::filesystem::path& to) {
    if (!to.parent_path().empty()) std::filesystem::create_directories(to.parent_path());
    std::filesystem::copy_file(from, to, std::filesystem::copy_options::overwrite_existing);
}

void remove_if_exists(const std::filesystem::path& path) {
    if (std::filesystem::exists(path)) std::filesystem::remove(path);
}

} // namespace cg
