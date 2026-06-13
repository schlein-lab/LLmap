// LLmap — Splice-site (GT/AG) boundary snapping for Transcript-Mode.
//
// The seed-chain path resolves an exon boundary only to ~seed-window precision:
// the upstream exon's alignment ends a few bp short of the true donor and the
// downstream exon starts a few bp after the true acceptor, leaving an unaligned
// read gap (q_gap) in between. The real-data minimap2 benchmark showed the two
// costs of that: (1) the intron N length is inflated by the unaligned flank, so
// the junction is not at the true splice site; (2) the residual q_gap makes the
// joiner's geometry gate fragile on error-bearing reads (junctions don't merge,
// → only 1/24 reads spliced vs minimap2's 18/24).
//
// Snapping fixes both. Given two adjacent same-(ref,strand) sub-chains with an
// intron-sized reference gap and a small read gap [a.query_end, b.query_start],
// it picks the single read split point `query_split` and the matching reference
// donor/acceptor positions such that the intron is canonical — sense-strand GT
// at the donor (intron 5') and AG at the acceptor (intron 3'). Extending the
// upstream exon to query_split and the downstream exon back to query_split
// closes the q_gap to zero and sets the N to the true intron length. On the '-'
// strand the forward-reference motif is the reverse-complement (CT..AC) so the
// transcript still reads GT..AG in sense.
//
// Pure (stdlib only); the joiner injects the reference sequence and applies the
// returned positions (growing/shrinking the boundary M ops). A no-op (snapped =
// false, boundary kept as-is) when no canonical site sits in the window, so the
// junction confidence / jM tag reflects the non-canonical call rather than a
// silent forced merge.

#pragma once

#include "mapping/chain_spliced.h"

#include <cstdint>
#include <string_view>

namespace llmap::mapping {

struct SpliceSnap {
    bool          snapped{false};       // canonical GT..AG found in the window
    std::uint64_t donor_ref_pos{0};     // first intron base (a-side), genomic
    std::uint64_t acceptor_ref_pos{0};  // one-past last intron base (b-side)
    std::uint32_t query_split{0};       // read point: exon-a ends / exon-b begins
    float         motif_score{0.0f};    // 1.0 canonical, 0.0 non-canonical (for jM)
};

// Snap the boundary between adjacent sub-chains `a` (upstream) and `b`
// (downstream) — same ref+strand, `a.ref_end <= b.ref_start` — to the canonical
// splice site nearest the seed boundary. Searches read split points in
// [a.query_end - window, b.query_start + window]; for each, the donor follows
// the upstream exon's end and the acceptor the downstream exon's start, so the
// read stays contiguous (q_gap = 0). Returns the smallest-displacement canonical
// placement, or {snapped=false} with the boundary unchanged if none is found.
[[nodiscard]] SpliceSnap SnapJunction(const LinearSubChain& a,
                                      const LinearSubChain& b,
                                      std::string_view ref_seq,
                                      char strand,
                                      std::uint32_t window = 30);

}  // namespace llmap::mapping
