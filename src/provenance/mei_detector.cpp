// LLmap — Mobile-element-insertion (MEI) signature detector implementation.

#include "provenance/mei_detector.h"

#include <algorithm>

namespace llmap::provenance {

const char* TeFamilyTag(TeFamily f) noexcept {
    switch (f) {
        case TeFamily::Alu:  return "alu";
        case TeFamily::L1:   return "l1";
        case TeFamily::Sva:  return "sva";
        case TeFamily::Herv: return "herv";
        case TeFamily::Unknown: return "te";
    }
    return "te";
}

MeiCall ClassifyMei(const MeiEvidence& e) {
    MeiCall c;
    c.te_family = e.te_family;
    c.has_polyA = e.has_polyA;
    c.has_tsd = e.tsd_length >= 4 && e.tsd_length <= 20;  // canonical TSD window

    // Core signature: a split read whose clipped portion is a TE consensus.
    if (!e.is_split || !e.clip_matches_te) return c;  // not an MEI

    float conf = 0.40f;                       // split + TE-consensus clip
    if (c.has_polyA) conf += 0.30f;           // poly-A/T tail — strong L1/Alu hallmark
    if (c.has_tsd) conf += 0.20f;             // target-site duplication
    if (e.five_prime_truncated) conf += 0.10f;  // L1 5'-truncation
    if (e.at_reference_te) conf -= 0.15f;     // bare clip at a ref TE ⇒ likely confusion
    conf = std::clamp(conf, 0.0f, 1.0f);
    c.confidence = conf;

    // A NOVEL insertion needs at least one insertion hallmark beyond the bare
    // TE-clip (poly-A or a canonical-length TSD); a bare clip alone defers to the
    // Layer-1 reference-TE-confusion path (mapping_confusion `mei`).
    c.is_mei = (c.has_polyA || c.has_tsd) && conf >= 0.60f;
    return c;
}

}  // namespace llmap::provenance
