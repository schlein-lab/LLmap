// LLmap — Per-position splice determinism implementation.

#include "mapping/splice_determinism.h"

#include <algorithm>
#include <cstdint>
#include <utility>

namespace llmap::mapping {

namespace {

// Tokenise a CIGAR into (length, op) pairs.
std::vector<std::pair<std::uint64_t, char>> Parse(std::string_view cigar) {
    std::vector<std::pair<std::uint64_t, char>> ops;
    std::uint64_t num = 0;
    for (const char c : cigar) {
        if (c >= '0' && c <= '9') {
            num = num * 10 + static_cast<std::uint64_t>(c - '0');
        } else {
            ops.emplace_back(num, c);
            num = 0;
        }
    }
    return ops;
}

bool ConsumesRefMatch(char op) { return op == 'M' || op == '=' || op == 'X'; }

}  // namespace

void DeterminismAccumulator::MarkIncluded(const std::string& ref_id,
                                          std::uint64_t pos0,
                                          std::string_view cigar) {
    auto& posmap = positions[ref_id];
    auto& juncs = junctions[ref_id];
    std::uint64_t cur = pos0;
    for (const auto& [len, op] : Parse(cigar)) {
        if (ConsumesRefMatch(op)) {
            for (std::uint64_t p = cur; p < cur + len; ++p) ++posmap[p].included;
            cur += len;
        } else if (op == 'N') {
            ++juncs[{cur, cur + len}];  // donor = cur, acceptor = cur + len
            cur += len;
        } else if (op == 'D') {
            cur += len;  // ref consumed, neither included nor excluded
        }
        // I / S / H / P: no reference advance.
    }
}

void DeterminismAccumulator::MarkExcluded(const std::string& ref_id,
                                          std::uint64_t pos0,
                                          std::string_view cigar) {
    const auto ref_it = positions.find(ref_id);
    if (ref_it == positions.end()) return;
    auto& posmap = ref_it->second;
    std::uint64_t cur = pos0;
    for (const auto& [len, op] : Parse(cigar)) {
        if (ConsumesRefMatch(op)) {
            cur += len;
        } else if (op == 'N') {
            // Only existing (exonic-somewhere) positions inside the intron.
            auto lo = posmap.lower_bound(cur);
            const auto hi = posmap.lower_bound(cur + len);
            for (; lo != hi; ++lo) ++lo->second.excluded;
            cur += len;
        } else if (op == 'D') {
            cur += len;
        }
    }
}

double Determinism(const PositionCounts& c) noexcept {
    const std::uint32_t total = c.included + c.excluded;
    if (total == 0) return 0.0;
    const std::uint32_t modal = std::max(c.included, c.excluded);
    return 100.0 * static_cast<double>(modal) / static_cast<double>(total);
}

std::vector<JunctionUsage> SortedJunctions(const DeterminismAccumulator& acc,
                                           const std::string& ref_id) {
    std::vector<JunctionUsage> out;
    const auto it = acc.junctions.find(ref_id);
    if (it == acc.junctions.end()) return out;
    out.reserve(it->second.size());
    for (const auto& [key, n] : it->second) {  // map is ordered by (donor,acceptor)
        out.push_back(JunctionUsage{key.first, key.second, n});
    }
    return out;
}

}  // namespace llmap::mapping
