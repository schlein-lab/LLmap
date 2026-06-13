// LLmap — Multi-level consensus resolution for the adaptive cascade.
//
// Operator refinement: a read is "sufficiently mapped" when N consecutive
// cascade LEVELS converge on the same placement. Cross-validation across
// independent levels is more robust than a single-stage criterion — it directly
// guards against the false-early-exit risk (a read called easy that was actually
// spliced/ambiguous). The levels already exist as the bucket pyramid (coarse L2
// → fine L0) and the feature cascade (chain → EM → …); agreement of the last K
// is the stop criterion.
//
// Pure function over the ordered per-level placements; the dispatcher records
// each level's result and asks this when to stop. Decoupled (stdlib only).

#pragma once

#include <cstdint>
#include <span>
#include <string>

namespace llmap::mapper {

// One level's placement decision for a read.
struct LevelPlacement {
    bool          mapped{false};   // did this level place the read?
    std::string   ref_id;          // chosen reference (empty if unmapped)
    std::uint64_t pos{0};          // chosen position (0-based)
};

struct ConsensusConfig {
    std::uint32_t agree_k{2};        // how many consecutive levels must agree
    std::uint64_t pos_tolerance{10}; // bp slack when comparing positions
};

struct ConsensusResult {
    bool        resolved{false};      // last agree_k levels converged
    std::uint32_t agreeing_levels{0}; // length of the trailing agreeing run
};

// Resolved iff the LAST `agree_k` level placements are the same (same ref, same
// pos within tolerance, all mapped). Walking from the most recent level back,
// count the trailing run of agreeing placements; resolved when it reaches
// agree_k. Fewer than agree_k levels, or a disagreement/unmapped in the trailing
// window, → not resolved (escalate / keep applying features).
[[nodiscard]] ConsensusResult IsConsensusResolved(
    std::span<const LevelPlacement> levels, const ConsensusConfig& cfg = {});

}  // namespace llmap::mapper
