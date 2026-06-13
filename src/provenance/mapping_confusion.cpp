// LLmap — Mapping-confusion provenance detector implementation.

#include "provenance/mapping_confusion.h"

#include <algorithm>

namespace llmap::provenance {

const char* MappingConfusionTag(MappingConfusion m) noexcept {
    switch (m) {
        case MappingConfusion::None:       return "none";
        case MappingConfusion::Paralog:    return "para";
        case MappingConfusion::Numt:       return "numt";
        case MappingConfusion::Pseudogene: return "pseudo";
        case MappingConfusion::Rdna:       return "rdna";
    }
    return "none";
}

MappingConfusionCall ClassifyMappingConfusion(const MappingConfusionEvidence& e) {
    // 1. Processed pseudogene: at a known parent/pseudogene pair, an intronless
    //    read where the parent gene is intron-bearing is pseudogene-derived.
    if (e.at_pseudogene_parent_locus && !e.read_is_spliced) {
        return MappingConfusionCall{MappingConfusion::Pseudogene, 0.70f};
    }

    // 2. NUMT artifact: at an mt-homologous locus, if the read matches a nuclear
    //    NUMT copy better than the mitochondrial reference, it is nuclear DNA
    //    masquerading as mt → the classic fake heteroplasmy. (The opposite case,
    //    real mt, is host biology → IsRealMtHeteroplasmy / Layer-3, not here.)
    if (e.at_mt_homologous_locus &&
        e.identity_to_nuclear_numt > e.identity_to_mt + e.mt_identity_margin) {
        const float conf =
            std::clamp((e.identity_to_nuclear_numt - e.identity_to_mt) * 5.0f,
                       0.0f, 1.0f);
        return MappingConfusionCall{MappingConfusion::Numt, conf};
    }

    // 3. Paralog / segdup: an ambiguous PSV posterior means the assigned copy is
    //    not confidently the true one.
    if (e.has_psv && e.psv_posterior < e.paralog_posterior_threshold) {
        return MappingConfusionCall{MappingConfusion::Paralog,
                                    1.0f - e.psv_posterior};
    }

    // 4. rDNA / satellite / HOR array: non-unique placement (low MAPQ in a
    //    known repeat array).
    if (e.in_repeat_array && e.mapq < 10) {
        return MappingConfusionCall{MappingConfusion::Rdna, 0.50f};
    }

    return MappingConfusionCall{};  // None → host / correctly placed
}

bool IsRealMtHeteroplasmy(const MappingConfusionEvidence& e) {
    // Genuine mt read: sits at an mt-homologous locus and matches the mt
    // reference at least as well as (margin) any nuclear NUMT copy. Mutually
    // exclusive with a Numt call (which requires nuclear to win by the margin).
    return e.at_mt_homologous_locus &&
           e.identity_to_mt >= e.identity_to_nuclear_numt + e.mt_identity_margin;
}

}  // namespace llmap::provenance
