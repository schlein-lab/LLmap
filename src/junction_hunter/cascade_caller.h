// LLmap — junction_hunter: variable-tier cascade caller.

#pragma once

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/cascade_pair_index.h"

#include <cstdint>
#include <vector>

namespace llmap::junction_hunter {

struct CascadeConfig {
    /// k for each cascade tier, smallest first. tier_class[i] uses k = k_values[i].
    std::vector<std::uint8_t> k_values;
    /// Minimum tier-0 anchor hits (across LcrUp ∪ LcrDown ∪ Interior)
    /// to promote a (read, pair) tuple to higher tiers.
    std::uint16_t tier1_min_anchor_hits{3};
    /// Consensus min agreement count for class call per read position.
    std::uint8_t consensus_min{3};
    /// Monotonicity threshold (Spearman ρ; ≥ this = monotonic).
    float monotonicity_min{0.95f};

    /// Long-read preset (HiFi/ONT 10 kb+ reads).
    /// k = {11, 17, 25, 35, 51, 71, 101, 125} — broad k spectrum, room
    /// for both cheap membership tests and PSV-grade discrimination.
    static CascadeConfig LongReadPreset();

    /// Short-read preset (Illumina paired-end 150 bp).
    /// k = {11, 13, 15, 17, 19, 21, 25} — only small k that fit in
    /// a 150 bp read with enough remaining positions to vote.
    static CascadeConfig ShortReadPreset();
};

JunctionRecord CascadeCall(std::string_view read_id,
                           std::string_view read_seq,
                           const CascadePairIndex& pair_index,
                           const CascadeConfig& cfg);

}  // namespace llmap::junction_hunter
