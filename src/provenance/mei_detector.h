// LLmap — Mobile-element-insertion (MEI) signature detector — Block 2 (i).
//
// Detects a NON-reference TE insertion (polymorphic or somatic Alu/L1/SVA/HERV)
// from the read signature, and routes it to the Layer-3 bioconfounder `bio:mei`
// — real host biology that looks like an SV, NOT a contamination artifact (same
// lossless doctrine as VDJ / numt-vs-mthet / bio:sv). A confident MEI is also a
// FINDING (somatic L1 retrotransposition is clinically relevant), not a
// throwaway flag.
//
// CRITICAL split from reference-TE confusion (mapping_confusion `mei` Layer-1):
//   * A read landing in / mis-mapped across the ~1M REFERENCE TE copies →
//     Layer-1 `mei` ARTIFACT (mapping_confusion).
//   * A read spanning a NOVEL insertion absent from the reference (flank +
//     TE-consensus clip + TSD + poly-A) → Layer-3 `bio:mei` (this module).
// The diagnostic signature for a novel insertion: a split/soft-clipped read
// whose clipped portion matches a TE consensus, a target-site duplication (TSD,
// ~4-20 bp flanking repeat), and a poly-A/poly-T tail (L1/Alu hallmark),
// optionally 5'-truncation (L1).
//
// Pure classifier over an explicit evidence struct (like the other detector
// classes) — decoupled from the TE-family catalog (Block 2 (ii)), which fills
// clip_matches_te / te_family.

#pragma once

#include <cstdint>

namespace llmap::provenance {

enum class TeFamily : std::uint8_t { Unknown = 0, Alu, L1, Sva, Herv };

[[nodiscard]] const char* TeFamilyTag(TeFamily f) noexcept;  // "alu"/"l1"/"sva"/"herv"/"te"

struct MeiEvidence {
    bool     is_split{false};            // soft-clip / SA: genomic flank + clipped portion
    bool     clip_matches_te{false};     // clipped portion matches a TE consensus (catalog ii)
    TeFamily te_family{TeFamily::Unknown};
    bool     has_polyA{false};           // poly-A/T tail at the junction (L1/Alu hallmark)
    std::uint32_t tsd_length{0};         // target-site duplication length (bp); 4-20 typical
    bool     five_prime_truncated{false};// L1 5'-truncation signature
    // True iff the insertion site is itself an ANNOTATED reference TE: then a
    // bare TE-clip is more likely reference-confusion (Layer-1) than a novel
    // insertion, so we demand stronger insertion hallmarks before calling MEI.
    bool     at_reference_te{false};
};

struct MeiCall {
    bool     is_mei{false};        // non-reference MEI → Layer-3 `bio:mei` + finding
    TeFamily te_family{TeFamily::Unknown};
    float    confidence{0.0f};     // [0,1]
    bool     has_tsd{false};       // TSD length in the canonical 4-20 bp window
    bool     has_polyA{false};
};

// A novel MEI requires the split + TE-consensus clip plus at least one insertion
// hallmark (poly-A or a canonical-length TSD). Confidence accumulates over the
// hallmarks present; at_reference_te discounts it (bare clip there = likely
// reference-confusion, deferred to mapping_confusion's Layer-1 `mei`).
[[nodiscard]] MeiCall ClassifyMei(const MeiEvidence& e);

}  // namespace llmap::provenance
