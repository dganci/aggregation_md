#include "FileUtils.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cgcv {

void require_file(const std::filesystem::path& path, const std::string& label) {
    if (!std::filesystem::is_regular_file(path)) throw std::runtime_error("Missing " + label + ": " + path.string());
}

void ensure_dir(const std::filesystem::path& path) {
    std::error_code ec;
    std::filesystem::create_directories(path, ec);
    if (ec) throw std::runtime_error("Cannot create directory: " + path.string() + " (" + ec.message() + ")");
}

void write_text(const std::filesystem::path& path, const std::string& text) {
    ensure_dir(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write file: " + path.string());
    out << text;
}

std::string read_text(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read file: " + path.string());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

} // namespace cgcv
