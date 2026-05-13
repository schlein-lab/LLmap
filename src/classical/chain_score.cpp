#include "classical/chain.h"

#include <algorithm>
#include <cmath>

namespace llmap::classical {

int32_t AnchorPairScore(
    const Anchor& a,
    const Anchor& b,
    const ChainConfig& config)
{
    // Anchors must be on same ref and strand
    if (a.ref_id != b.ref_id || a.same_strand != b.same_strand) {
        return std::numeric_limits<int32_t>::min();
    }

    // b must come after a in both coordinates for forward chaining
    int64_t ref_gap = static_cast<int64_t>(b.ref_pos) - static_cast<int64_t>(a.ref_pos);
    int64_t query_gap = static_cast<int64_t>(b.query_pos) - static_cast<int64_t>(a.query_pos);

    if (a.same_strand) {
        // Forward strand: both gaps must be positive
        if (ref_gap <= 0 || query_gap <= 0) {
            return std::numeric_limits<int32_t>::min();
        }
    } else {
        // Reverse strand: ref_gap positive, query_gap negative
        if (ref_gap <= 0 || query_gap >= 0) {
            return std::numeric_limits<int32_t>::min();
        }
        query_gap = -query_gap;
    }

    // Check gap limits
    if (static_cast<uint64_t>(ref_gap) > config.max_gap_ref ||
        static_cast<uint64_t>(query_gap) > config.max_gap_query) {
        return std::numeric_limits<int32_t>::min();
    }

    // Score = match - gap_penalty * |gap_difference|
    // Gap difference = abs(ref_gap - query_gap)
    int64_t gap_diff = std::abs(ref_gap - query_gap);

    // Minimap2-style scoring: min(gap) bases are matches
    int64_t min_gap = std::min(ref_gap, query_gap);
    int64_t score = min_gap * config.match_score - gap_diff * config.gap_penalty;

    return static_cast<int32_t>(std::clamp(
        score,
        static_cast<int64_t>(std::numeric_limits<int32_t>::min()),
        static_cast<int64_t>(std::numeric_limits<int32_t>::max())));
}

bool IsColinear(
    const Anchor& a,
    const Anchor& b,
    const ChainConfig& config)
{
    return AnchorPairScore(a, b, config) > std::numeric_limits<int32_t>::min();
}

}  // namespace llmap::classical
