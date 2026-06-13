// LLmap — TE-scale WaveCollapse spread-mass (Block 2 (iii)).
//
// THE intellectually-satisfying core, where LLmap structurally beats minimap2/
// BWA on transposable elements. A read landing inside a TE family (~1.1M Alu,
// ~500k L1 near-identical copies) has NO honest unique placement. Standard
// mappers either pick one copy (a FAKED unique hit, often wrong, with a
// misleading MAPQ) or drop it at MAPQ0. LLmap's lossless wave-particle doctrine:
//
//   * keep probability MASS over the candidate copies (don't collapse to argmax),
//   * if part of the read has a UNIQUE anchor (a flank reaching out of the TE),
//     collapse via that anchor to the true locus,
//   * report MAPQ = the HONEST ambiguity (entropy of the spread), not a fake 60.
//
// This turns the ~2–6.5% repeat-ambiguous reads the matched HG002 test surfaced
// from silent mis-placements into explicit, quantified positional uncertainty —
// and feeds the provenance `mei`/`paralog` Layer-1 confidence with a real number.
//
// Pure (stdlib only). The EM/integration layer supplies the candidate placements
// (from the minimizer index + TE family catalog); this resolves the mass.

#pragma once

#include "provenance/te_family_catalog.h"

#include <cstdint>
#include <string>
#include <vector>

namespace llmap::provenance {

// One candidate placement of a read (a TE-family copy it aligns to, or a
// unique-anchored locus).
struct TePlacement {
    std::string   chrom;
    std::uint64_t pos{0};
    std::int32_t  align_score{0};   // alignment score at this placement
    bool          has_unique_anchor{false};  // a read flank maps uniquely here
    TeClass       te_class{TeClass::OtherRepeat};
    std::string   subfamily;

    // Population prior (Operator's directive): how consistently this placement is
    // sensibly mappable across hundreds of pangenome genomes — from the per-locus
    // mappability baseline. 1.0 = sharply mappable everywhere (trustworthy → mass
    // boost); ~0 = this placement is rarely the right one across the population;
    // <0 sentinel (default) = no baseline → ignored (reference-only behaviour).
    // Blurriness that is CONSISTENT across the pangenome is a real locus property
    // (honest spread); blurriness only HERE is an anomaly (the caller flags it).
    float         population_support{-1.0f};
};

struct SpreadMass {
    std::vector<float> posterior;   // per-placement, sums to 1 (the kept mass)
    int   collapsed{-1};            // index collapsed-to via a unique anchor, else -1
    float confidence{0.0f};         // 1 - normalised entropy ∈ [0,1]
    std::uint8_t mapq{0};           // HONEST MAPQ from the best posterior (not faked)
    bool  is_ambiguous{false};      // true when no unique anchor → mass stays spread
};

// Resolve the kept mass over candidate placements. With a unique anchor, collapse
// to the best-scoring anchored placement (confident). Without, softmax the scores
// (temperature τ) and keep the spread — confidence/MAPQ then reflect the TRUE
// multi-copy ambiguity. Empty input → zero-confidence, ambiguous.
[[nodiscard]] SpreadMass ResolveTeSpreadMass(const std::vector<TePlacement>& cands,
                                             float temperature = 1.0f);

}  // namespace llmap::provenance
