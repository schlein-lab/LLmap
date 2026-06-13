// LLmap — Chimera provenance detector (Layer-1 `chim`) — class D.
//
// Thin wrapper over the Block-7 chimera detector (src/chimera). Feeds the
// resolver's ReadEvidence.is_chimera: a read whose aligned parts form a
// ligation/PCR chimera, ONT fused-read, or cross-chromosomal artifact is an
// ARTIFACT origin (Layer-1 `chim`).
//
// Critical separation (same lossless doctrine as Numt vs real mt-heteroplasmy):
// a split read is not automatically an ARTIFACT. The matched HG002 long-vs-short
// test made this concrete — HiFi shows ~1.7 % SA-tag splits vs ~0 % on Illumina,
// the OPPOSITE of an artifact signal: long reads simply SPAN real structural
// breakpoints. So we route three ways, not one:
//   * VDJ class-switch / recombination (kind 'V')          → Layer-3 `bio:vdj`
//     (real immune biology — IGH/TCR somatic rearrangement)
//   * intra-chromosomal split (kind 'I', SV-spanning)       → Layer-3 `bio:sv`
//     (real structural variant — deletion/inversion/dup breakpoint a long read
//      crosses; overwhelmingly biology, not a ligation artifact)
//   * cross-chromosomal split (kind 'X')                    → Layer-1 `chim`
//     (artifact-leaning — ligation/PCR chimera; a true translocation is rarer)
// Only the last is a Layer-1 artifact origin; the first two are host biology
// (flag, never bucket — same rule as VDJ / numt-vs-mthet). Routing the resolver
// reads off the three flags below.

#pragma once

#include "chimera/chimera_detector.h"

#include <cstdint>
#include <span>

namespace llmap::provenance {

struct ChimeraCall {
    bool  is_chimera{false};            // cross-chrom artifact → Layer-1 `chim`
    bool  is_vdj_recombination{false};  // real immune biology → Layer-3 `bio:vdj`
    bool  is_sv_spanning{false};        // real structural variant → Layer-3 `bio:sv`
    float confidence{0.0f};             // [0,1] for the call
    char  kind{'.'};                    // 'I' intra / 'X' cross-chrom / 'V' vdj
};

// Classify a read's aligned parts into exactly one of: artifact chimera
// (cross-chrom 'X' → is_chimera, Layer-1), VDJ recombination ('V' →
// is_vdj_recombination, Layer-3), SV-spanning intra-chromosomal split ('I' →
// is_sv_spanning, Layer-3). No chimeric hypothesis → all false (host / partial).
// cfg defaults match Block-7. The intra-chrom 'I' default to real SV is the
// lossless choice (prefer flagging real biology over bucketing it as artifact);
// a rare intra-region PCR chimera is the documented false-positive.
[[nodiscard]] ChimeraCall ClassifyChimera(
    std::span<const chimera::AlignedPart> parts,
    const chimera::VdjLocusMask& vdj_mask,
    const chimera::ChimeraConfig& cfg = {});

}  // namespace llmap::provenance
