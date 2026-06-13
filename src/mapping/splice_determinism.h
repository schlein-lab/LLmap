// LLmap — Per-position splice determinism D(pos) + junction usage.
//
// Operator metric (2026-06-13): overlay all spliced read mappings and, per
// genomic position, quantify how *determined* the splice state is — i.e. how
// consistently reads treat that position. A constitutive exon base that every
// read includes scores 100; a cassette-exon base excluded in 9/10 reads scores
// 90 (modal state = excluded, 90 % agree); a base near a blurry 5'/3' site that
// reads splice at slightly different bp scores lower still. Formally
//
//     D(pos) = 100 * max(included, excluded) / (included + excluded)
//
// = 100 * (1 - normalised splice entropy) at that position. This is the
// splicing analogue of LLmap's lossless doctrine: ambiguity is *quantified*,
// not hidden, and it is the natural per-position aggregate of the WaveCollapse
// collapse confidence. Variant-useful: a variant at D≈100 sits in deterministic
// context, one at D≈76 in a natively variable / blurry locus — different
// functional prior (the basis for a future GTEx sQTL/eQTL comparison).
//
// Tool-agnostic: consumes any spliced SAM (llmap OR minimap2 -ax splice), so the
// same metric compares mappers and, later, against GENCODE boundaries.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace llmap::mapping {

// Per-position tally across reads.
struct PositionCounts {
    std::uint32_t included{0};  // reads aligning exon body (M/=/X) over this pos
    std::uint32_t excluded{0};  // reads skipping this pos inside an intron (N)
};

// One observed intron and how many reads splice exactly there.
struct JunctionUsage {
    std::uint64_t donor{0};      // first intron base, 0-based genomic
    std::uint64_t acceptor{0};   // one-past last intron base
    std::uint32_t n_reads{0};
};

// Accumulator over a set of spliced alignments (one reference at a time is fine;
// keyed by ref_id so multi-contig input works).
struct DeterminismAccumulator {
    // ref_id → (genomic pos → counts), sparse (only covered positions stored).
    std::map<std::string, std::map<std::uint64_t, PositionCounts>> positions;
    // ref_id → (packed (donor<<32|acceptor-ish key) → usage). Stored as a map of
    // (donor,acceptor) → count via an ordered map keyed by the pair.
    std::map<std::string, std::map<std::pair<std::uint64_t, std::uint64_t>,
                                   std::uint32_t>> junctions;

    // PASS 1: walk one alignment's CIGAR from 0-based genomic `pos0`, marking
    // every M/=/X position as INCLUDED and recording each N span as a junction.
    // This builds the set of positions that are exonic in at least one read —
    // the only positions where the included/excluded contrast is meaningful.
    void MarkIncluded(const std::string& ref_id, std::uint64_t pos0,
                      std::string_view cigar);

    // PASS 2 (run after all MarkIncluded): for each N span in the CIGAR,
    // increment EXCLUDED for the already-present (exonic-somewhere) positions
    // inside the intron. Range-restricted to existing map entries, so a 100 kb
    // intron costs O(exonic positions within it), not O(intron length).
    void MarkExcluded(const std::string& ref_id, std::uint64_t pos0,
                      std::string_view cigar);
};

// D(pos) in [0, 100] from a position's counts. Returns 0 for an uncovered
// position (included+excluded == 0).
[[nodiscard]] double Determinism(const PositionCounts& c) noexcept;

// Flatten the accumulated junctions for one ref into a sorted vector (by donor
// then acceptor), for TSV / PSI reporting.
[[nodiscard]] std::vector<JunctionUsage> SortedJunctions(
    const DeterminismAccumulator& acc, const std::string& ref_id);

}  // namespace llmap::mapping
