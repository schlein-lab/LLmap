// LLmap — Multi-level consensus resolution implementation.

#include "mapper/consensus_resolution.h"

namespace llmap::mapper {

namespace {

bool SamePlacement(const LevelPlacement& a, const LevelPlacement& b,
                   std::uint64_t tol) {
    if (!a.mapped || !b.mapped) return false;
    if (a.ref_id != b.ref_id) return false;
    const std::uint64_t d = a.pos > b.pos ? a.pos - b.pos : b.pos - a.pos;
    return d <= tol;
}

}  // namespace

ConsensusResult IsConsensusResolved(std::span<const LevelPlacement> levels,
                                    const ConsensusConfig& cfg) {
    ConsensusResult r;
    const std::uint32_t k = cfg.agree_k == 0 ? 1 : cfg.agree_k;
    if (levels.size() < k) return r;

    const LevelPlacement& last = levels.back();
    if (!last.mapped) return r;  // an unmapped last level can't resolve

    // Count the trailing run of placements agreeing with the most recent one.
    std::uint32_t run = 1;
    for (std::size_t i = levels.size() - 1; i > 0; --i) {
        if (SamePlacement(levels[i - 1], last, cfg.pos_tolerance)) {
            ++run;
        } else {
            break;
        }
    }
    r.agreeing_levels = run;
    r.resolved = run >= k;
    return r;
}

}  // namespace llmap::mapper
