// LLmap — Pangenome mapping-determinism prior implementation.

#include "core/pangenome_mappability.h"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace llmap::core {

void PangenomeMappability::AddSample(std::string_view ref_id, std::uint64_t pos,
                                     float confidence) {
    confidence = std::clamp(confidence, 0.0f, 1.0f);
    const std::uint64_t win = pos / window_bp_;
    Cell& c = track_[std::string(ref_id)][win];
    c.sum += confidence;
    ++c.n;
}

float PangenomeMappability::PopulationSupport(std::string_view ref_id,
                                              std::uint64_t pos,
                                              float neutral) const {
    const auto rit = track_.find(std::string(ref_id));
    if (rit == track_.end()) return neutral;
    const auto wit = rit->second.find(pos / window_bp_);
    if (wit == rit->second.end() || wit->second.n == 0) return neutral;
    return static_cast<float>(wit->second.sum /
                              static_cast<double>(wit->second.n));
}

std::uint32_t PangenomeMappability::SampleCount(std::string_view ref_id,
                                                std::uint64_t pos) const {
    const auto rit = track_.find(std::string(ref_id));
    if (rit == track_.end()) return 0;
    const auto wit = rit->second.find(pos / window_bp_);
    return wit == rit->second.end() ? 0u : wit->second.n;
}

std::size_t PangenomeMappability::n_windows() const noexcept {
    std::size_t n = 0;
    for (const auto& [ref, wins] : track_) n += wins.size();
    return n;
}

bool PangenomeMappability::Save(const std::string& path) const {
    std::ofstream out(path);
    if (!out) return false;
    out << "# ref\tstart\tend\tM\tn_samples\twindow_bp=" << window_bp_ << "\n";
    for (const auto& [ref, wins] : track_) {
        for (const auto& [win, cell] : wins) {
            if (cell.n == 0) continue;
            const std::uint64_t start = win * window_bp_;
            const double m = cell.sum / static_cast<double>(cell.n);
            out << ref << '\t' << start << '\t' << (start + window_bp_) << '\t'
                << m << '\t' << cell.n << '\n';
        }
    }
    return static_cast<bool>(out);
}

bool PangenomeMappability::Load(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::map<std::string, std::map<std::uint64_t, Cell>> loaded;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string ref;
        std::uint64_t start = 0, end = 0;
        double m = 0.0;
        std::uint32_t n = 0;
        if (!(ss >> ref >> start >> end >> m >> n)) continue;
        Cell c;
        c.sum = m * static_cast<double>(n);  // reconstruct sum from mean × n
        c.n = n;
        loaded[ref][start / window_bp_] = c;
    }
    if (loaded.empty()) return false;
    track_ = std::move(loaded);
    return true;
}

}  // namespace llmap::core
