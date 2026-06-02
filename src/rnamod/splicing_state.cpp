// LLmap — SplicingStateClassifier implementation.
//
// Pure-functional logic: takes a vector of ObservedJunction values and
// returns a SplicingStateInference. No side effects, no I/O, fully
// testable. Read-state is inferred by inspecting the geometry + flags
// of the observed junctions:
//
//   * any back-splice (spliceosome_class == 3)       → BackSplicedCircular
//   * any retained intron                             → IntronRetained
//                                                       (or PartiallySpliced
//                                                        if multiple)
//   * no junctions at all                             → Unspliced
//   * all junctions annotated + canonical-class       → Canonical
//   * mixed annotated/un-annotated + long ref-gap    → RecursiveSpliced
//   * non-canonical-class junctions                   → NovelSplicingState
//
// Trans-splicing detection requires multi-chrom evidence; that comes
// from the ChimericDetail layer (AlignmentRecord), not from this
// classifier. We surface a TransSpliced state only when caller passes
// donor/acceptor pos from different chroms — encoded via flag at the
// junction level isn't possible without breaking ObservedJunction
// trivial-aggregate-ness. So TransSpliced is set externally by the
// chimera detector.

#include "rnamod/splicing_state.h"

#include <algorithm>

namespace llmap::splicing {

const char* SplicingStateName(SplicingState s) noexcept {
    switch (s) {
        case SplicingState::Unknown:                return "unknown";
        case SplicingState::Canonical:              return "canonical";
        case SplicingState::Unspliced:              return "unspliced";
        case SplicingState::PartiallySpliced:       return "partially_spliced";
        case SplicingState::IntronRetained:         return "intron_retained";
        case SplicingState::Lariat:                 return "lariat";
        case SplicingState::RecursiveSpliced:       return "recursive_spliced";
        case SplicingState::TransSpliced:           return "trans_spliced";
        case SplicingState::AlternativeCassetteIn:  return "alt_cassette_in";
        case SplicingState::AlternativeCassetteOut: return "alt_cassette_out";
        case SplicingState::Alt5ss:                 return "alt_5ss";
        case SplicingState::Alt3ss:                 return "alt_3ss";
        case SplicingState::MutuallyExclusive:      return "mutually_exclusive";
        case SplicingState::BackSplicedCircular:    return "back_spliced_circular";
        case SplicingState::HalfSplicedCotrans:     return "half_spliced_cotrans";
        case SplicingState::NovelSplicingState:     return "novel_splicing_state";
    }
    return "unknown";
}

SplicingStateInference SplicingStateClassifier::Classify(
    std::span<const ObservedJunction> junctions) const {

    SplicingStateInference inf;

    if (junctions.empty()) {
        inf.dominant = SplicingState::Unspliced;
        inf.confidence = 1.0f;
        return inf;
    }

    // ----- Tally features ------------------------------------------------
    std::size_t n_back   = 0;
    std::size_t n_retain = 0;
    std::size_t n_anno   = 0;
    std::size_t n_noncan = 0;
    std::uint32_t recursive_pos = 0;
    bool has_recursive = false;
    std::uint64_t max_ref_gap = 0;

    for (std::size_t i = 0; i < junctions.size(); ++i) {
        const auto& j = junctions[i];
        if (j.spliceosome_class == 3) {
            ++n_back;
        }
        if (j.is_retained) {
            ++n_retain;
            inf.retained_intron_indices.push_back(static_cast<std::uint32_t>(i));
        }
        if (j.is_annotated) {
            ++n_anno;
        }
        if (j.spliceosome_class == 2) {
            ++n_noncan;
        }

        // ref-gap (acceptor - donor); enough to flag recursive splicing if
        // a single junction spans >100 kb AND is un-annotated.
        const std::uint64_t lo = std::min(j.donor_genomic_pos,
                                            j.acceptor_genomic_pos);
        const std::uint64_t hi = std::max(j.donor_genomic_pos,
                                            j.acceptor_genomic_pos);
        const std::uint64_t gap = hi - lo;
        if (gap > max_ref_gap) max_ref_gap = gap;
        if (gap > 100'000 && !j.is_annotated) {
            has_recursive = true;
            recursive_pos = static_cast<std::uint32_t>(
                (j.donor_genomic_pos + j.acceptor_genomic_pos) / 2);
        }
    }

    // ----- Decide dominant ----------------------------------------------
    if (n_back > 0) {
        inf.dominant = SplicingState::BackSplicedCircular;
        inf.confidence = 0.95f;
        return inf;
    }
    if (n_retain >= 2) {
        inf.dominant = SplicingState::PartiallySpliced;
        inf.additional.push_back(SplicingState::IntronRetained);
        inf.confidence = 0.85f;
        return inf;
    }
    if (n_retain == 1) {
        inf.dominant = SplicingState::IntronRetained;
        inf.confidence = 0.90f;
        return inf;
    }
    if (has_recursive) {
        inf.dominant = SplicingState::RecursiveSpliced;
        inf.recursive_splice_site_pos = recursive_pos;
        inf.confidence = 0.70f;
        return inf;
    }
    if (n_noncan > 0 && n_noncan == junctions.size()) {
        inf.dominant = SplicingState::NovelSplicingState;
        inf.confidence = 0.50f;
        return inf;
    }
    if (n_anno == junctions.size()) {
        inf.dominant = SplicingState::Canonical;
        inf.confidence = 0.95f;
        return inf;
    }

    // Default — at least some junctions, mostly annotated, no special
    // signature → canonical-ish with reduced confidence.
    inf.dominant = SplicingState::Canonical;
    inf.confidence = 0.60f;
    return inf;
}

}  // namespace llmap::splicing
