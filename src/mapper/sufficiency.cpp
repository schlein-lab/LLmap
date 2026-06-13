// LLmap — Early-exit sufficiency criterion implementation.

#include "mapper/sufficiency.h"

namespace llmap::mapper {

SufficiencyResult IsProvablyResolved(std::span<const ChainSummary> chains,
                                     std::uint32_t query_len,
                                     const SufficiencyConfig& cfg) {
    SufficiencyResult r;
    if (chains.empty()) {
        r.reason = "no-chains";
        return r;
    }
    if (query_len == 0) {
        r.reason = "no-query-length";
        return r;
    }

    // Best + runner-up by score.
    std::size_t best_idx = 0;
    std::int32_t best_score = chains[0].score;
    std::int32_t second_score = 0;
    bool have_second = false;
    for (std::size_t i = 1; i < chains.size(); ++i) {
        const std::int32_t s = chains[i].score;
        if (s > best_score) {
            second_score = best_score;
            have_second = true;
            best_score = s;
            best_idx = i;
        } else if (!have_second || s > second_score) {
            second_score = s;
            have_second = true;
        }
    }
    const ChainSummary& best = chains[best_idx];

    // Intron-gap signature → this is a spliced read; it must escalate to the
    // splice path (with its own budget), never early-exit as "easy".
    if (best.has_intron_gap) {
        r.reason = "intron-signature";
        return r;
    }

    // The best chain must cover (almost) the whole read.
    const std::uint32_t covered =
        best.query_end > best.query_start ? best.query_end - best.query_start : 0;
    if (static_cast<double>(covered) <
        static_cast<double>(query_len) * cfg.min_query_coverage) {
        r.reason = "partial-coverage";
        return r;
    }

    // A competing locus of comparable score ⇒ genuinely ambiguous ⇒ escalate.
    if (have_second && second_score > 0 &&
        static_cast<double>(best_score) <
            cfg.dominance_factor * static_cast<double>(second_score)) {
        r.reason = "competing-locus";
        return r;
    }

    r.resolved = true;
    r.best_idx = best_idx;
    r.reason = "single-dominant";
    return r;
}

}  // namespace llmap::mapper
