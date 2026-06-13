// LLmap — WaveCollapse merge-ambiguity guard implementation.

#include "mapper/merge_ambiguity.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

namespace llmap::mapper {

namespace {

// Floor division (so negative and positive offsets bucket symmetrically).
std::int64_t FloorDiv(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b, r = a % b;
    if (r != 0 && ((r < 0) != (b < 0))) --q;
    return q;
}

}  // namespace

DiagonalClarity AssessDiagonal(std::span<const std::int64_t> offsets,
                               const AmbiguityConfig& cfg) {
    DiagonalClarity out;
    out.shared = static_cast<std::uint32_t>(offsets.size());
    if (offsets.empty()) return out;

    // Dominant diagonal: largest window of width 2*tol over the sorted offsets.
    std::vector<std::int64_t> v(offsets.begin(), offsets.end());
    std::sort(v.begin(), v.end());
    const std::int64_t width = 2 * cfg.offset_tolerance;
    std::size_t best = 0, l = 0;
    for (std::size_t r = 0; r < v.size(); ++r) {
        while (v[r] - v[l] > width) ++l;
        best = std::max(best, r - l + 1);
    }
    out.on_diagonal = static_cast<std::uint32_t>(best);
    out.dominant_fraction =
        static_cast<double>(best) / static_cast<double>(out.shared);

    // Offset entropy (diagnostic, D(pos)-style): low ⇒ collapsed, high ⇒ blurry.
    const std::int64_t bw = std::max<std::int64_t>(1, width + 1);
    std::map<std::int64_t, std::uint32_t> hist;
    for (const std::int64_t o : v) ++hist[FloorDiv(o, bw)];
    if (out.shared > 1) {
        double h = 0.0;
        for (const auto& [bucket, cnt] : hist) {
            const double p = static_cast<double>(cnt) / out.shared;
            h -= p * std::log(p);
        }
        out.offset_entropy =
            std::clamp(h / std::log(static_cast<double>(out.shared)), 0.0, 1.0);
    }

    out.clear = out.shared >= cfg.min_anchors &&
                out.dominant_fraction >= cfg.min_dominant_fraction;
    return out;
}

}  // namespace llmap::mapper
