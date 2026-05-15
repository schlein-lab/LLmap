// LLmap -- sorted-by-start interval index for per-base annotation lookup.
//
// Annotation intervals come from one or more layers (taxonomy, specific
// loci, agent). At lookup time we want the highest-priority (highest layer)
// interval whose [start, end) covers a given position.
//
// Implementation: per-contig vector of intervals sorted by start. Lookup is
// a binary search to find candidates, then a linear scan over overlapping
// intervals (typically <=5 on overlap, since intervals from different
// layers compose by priority -- we don't expect deep stacking).
//
// Overlay semantics: an interval's ParamOverride is composed over a base.
// AnnotationStore::ParamsAt(pos) walks layers high -> low and composes.

#pragma once

#include "annot/annot_types.h"

#include <span>
#include <unordered_map>
#include <vector>

namespace llmap::annot {

class AnnotationIndex {
public:
    AnnotationIndex() = default;

    // Replace the entire interval set with the given list.
    // Intervals from different ref_ids are partitioned internally.
    void Build(std::vector<AnnotationInterval> intervals);

    // All intervals that cover (ref_id, pos), sorted by layer (highest first).
    // Empty if no annotation present.
    std::vector<const AnnotationInterval*> IntervalsAt(uint32_t ref_id,
                                                       uint32_t pos) const;

    // Compose the ParamOverride for (ref_id, pos): layer-by-layer overlay
    // starting from highest layer down to lowest. Caller supplies a default
    // ParamOverride that fills any field unset by any covering interval.
    ParamOverride ParamsAt(uint32_t ref_id, uint32_t pos,
                           const ParamOverride& defaults = {}) const;

    // Iterate over every interval; useful for serialisation.
    std::vector<AnnotationInterval> Intervals() const;

    size_t Size() const { return total_count_; }
    bool Empty() const { return total_count_ == 0; }

private:
    // Per-contig sorted-by-start lists.
    std::unordered_map<uint32_t, std::vector<AnnotationInterval>> by_ref_;
    size_t total_count_ = 0;
};

}  // namespace llmap::annot
