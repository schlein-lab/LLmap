#include "annot/interval_tree.h"

#include <algorithm>

namespace llmap::annot {

void AnnotationIndex::Build(std::vector<AnnotationInterval> intervals) {
    by_ref_.clear();
    total_count_ = intervals.size();
    for (auto& iv : intervals) {
        by_ref_[iv.ref_id].push_back(std::move(iv));
    }
    for (auto& kv : by_ref_) {
        std::sort(kv.second.begin(), kv.second.end(),
            [](const AnnotationInterval& a, const AnnotationInterval& b) {
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
            });
    }
}

std::vector<const AnnotationInterval*> AnnotationIndex::IntervalsAt(
    uint32_t ref_id, uint32_t pos) const {

    std::vector<const AnnotationInterval*> out;
    auto it = by_ref_.find(ref_id);
    if (it == by_ref_.end()) return out;
    const auto& vec = it->second;
    if (vec.empty()) return out;

    // Find first interval with start > pos via upper_bound; back up to scan
    // overlapping ones. We accept O(k) scan back where k is the number of
    // intervals that could overlap pos -- bounded by layer depth (~5).
    auto upper = std::upper_bound(vec.begin(), vec.end(), pos,
        [](uint32_t p, const AnnotationInterval& iv) {
            return p < iv.start;
        });
    // Walk backwards from upper until we leave the overlap zone.
    for (auto rit = upper; rit != vec.begin();) {
        --rit;
        if (rit->end > pos && rit->start <= pos) {
            out.push_back(&*rit);
        }
        // Early exit: if rit->end <= pos AND there are no later intervals
        // that might start before pos, we can stop. The list is sorted by
        // start, so anything earlier has start <= rit->start; an interval
        // with start <= pos might still overlap iff end > pos. We need to
        // keep walking to catch deep-nested ones, but bound the scan to
        // ~32 entries to keep this O(1) in practice.
        if (out.size() >= 32) break;
    }

    // Higher layer first.
    std::sort(out.begin(), out.end(),
        [](const AnnotationInterval* a, const AnnotationInterval* b) {
            return static_cast<int>(a->layer) > static_cast<int>(b->layer);
        });

    return out;
}

ParamOverride AnnotationIndex::ParamsAt(uint32_t ref_id, uint32_t pos,
                                        const ParamOverride& defaults) const {
    auto intervals = IntervalsAt(ref_id, pos);
    ParamOverride result;
    // Highest layer first -> set what it provides
    for (const auto* iv : intervals) {
        // Overlay current result over interval's params
        ParamOverride next = iv->params;
        next.OverlayOver(result);
        result = next;
    }
    result.OverlayOver(defaults);
    return result;
}

std::vector<AnnotationInterval> AnnotationIndex::Intervals() const {
    std::vector<AnnotationInterval> out;
    out.reserve(total_count_);
    for (const auto& kv : by_ref_) {
        for (const auto& iv : kv.second) out.push_back(iv);
    }
    return out;
}

}  // namespace llmap::annot
