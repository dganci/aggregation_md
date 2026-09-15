#include "CvReadiness.hpp"
#include "SamplingMonitor.hpp"

#include <exception>
#include <fstream>
#include <string>

namespace cg {

CvReadiness read_cv_readiness(const std::filesystem::path& path) {
    CvReadiness cv;
    std::ifstream in(path);
    if (!in) {
        cv.note = "no output file";
        return cv;
    }
    std::string key, value;
    bool saw_its = false;
    while (in >> key >> value) {
        try {
            if (key == "plateau") cv.plateau = std::stoi(value) != 0;
            else if (key == "its1_ps") { cv.its1_ps = std::stod(value); saw_its = true; }
            else if (key == "total_time_ps") cv.total_time_ps = std::stod(value);
            else if (key == "note") { std::getline(in, cv.note); cv.note = value + cv.note; }
        } catch (const std::exception&) {
            cv.note = "malformed value for " + key;
            return cv;
        }
    }
    cv.ok = saw_its;
    if (!cv.ok && cv.note.empty()) cv.note = "no its1_ps in output";
    return cv;
}

void apply_cv_readiness(SamplingMetrics& m, const SamplingRules& rules, const CvReadiness& cv) {
    m.cv_plateau_found = cv.ok && cv.plateau;
    m.cv_its1_ps = cv.ok ? cv.its1_ps : 0.0;
    const double total_ps = m.total_time_us * 1'000'000.0;
    m.cv_time_over_its = (cv.ok && cv.its1_ps > 0.0) ? total_ps / cv.its1_ps : 0.0;
    decide_stop(m, rules);
}

} // namespace cg
