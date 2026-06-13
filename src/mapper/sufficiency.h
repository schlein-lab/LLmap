// LLmap — Early-exit sufficiency criterion for the adaptive feature cascade.
//
// The matched-data profile showed the bottleneck is EXTENSION, not the EM:
// ~17 WFA extensions per read, ~94% discarded on identity — we pay the full
// extend-everything price even when a read has one obvious answer. The operator's
// adaptive cascade fixes this: before extending all chains, ask whether the read
// is PROVABLY resolved by a single dominant chain; if so, extend only that one
// (1 extension instead of ~17) and exit. Otherwise escalate (extend top-k +
// WaveCollapse + bucket-routed features).
//
// This is NOT triage / a bypass that denies WaveCollapse to easy reads: the
// early-exit fires only on PROVABLE sufficiency — one colinear chain covering
// the whole read, dominating the runner-up, with no intron-gap signature and no
// second competitive locus. For such reads the full EM is provably redundant
// (it would return posterior≈1 on the same placement). The honesty is preserved
// for the reads that need it — they escalate. Lossless holds: no read is
// devalued, the easy read just doesn't pay the hard read's price.
//
// Pure function over compact chain summaries (no dep on the chainer's types) so
// it composes anywhere and is unit-testable; the pipeline fills ChainSummary
// from its chains and calls this between chaining and extension.

#pragma once

#include <cstdint>
#include <span>

namespace llmap::mapper {

// Compact, pre-extension view of one chain (filled from the chainer's output).
struct ChainSummary {
    std::int32_t  score{0};          // chain score (pre-extension)
    std::uint32_t query_start{0};    // covered read span [start, end)
    std::uint32_t query_end{0};
    bool          has_intron_gap{false};  // intron-sized ref gap inside the chain
                                          // (transcript signature → must escalate)
};

struct SufficiencyConfig {
    // The best chain must beat the runner-up by this factor (score-wise).
    float dominance_factor{1.5f};
    // The best chain must cover at least this fraction of the read.
    float min_query_coverage{0.90f};
};

struct SufficiencyResult {
    bool        resolved{false};   // provably sufficient → extend only `best_idx`
    std::size_t best_idx{0};
    const char* reason{"escalate"};  // why escalated (or "single-dominant" when resolved)
};

// Decide whether the read is provably resolved by one dominant chain. Resolved
// iff: a unique best chain covers >= min_query_coverage of the read, dominates
// the runner-up (best.score >= dominance_factor * second.score), and carries no
// intron-gap signature. Returns the best chain to extend; on escalation, the
// caller extends top-k and runs the full cascade.
[[nodiscard]] SufficiencyResult IsProvablyResolved(
    std::span<const ChainSummary> chains, std::uint32_t query_len,
    const SufficiencyConfig& cfg = {});

}  // namespace llmap::mapper
