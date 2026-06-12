// LLmap — Transcript-Mode spliced stage.
//
// Post-alignment stage that turns the classical chainer's per-read output
// (several LinearSubChains, one per exon block, because the linear chainer's
// hard gap limit breaks a spliced read at every intron) into spliced
// alignments: adjacent sub-chains separated by an intron-sized reference gap
// are joined (chain_spliced::JoinSplicedChains) into one alignment with an
// N-op CIGAR and per-junction confidence.
//
// Design notes:
//   * Dependency-light: this module knows only `mapping` types + stdlib. The
//     junction scorer (real one is annot::SpliceSiteDb) is INJECTED as a
//     functor, so the module needs no annot dependency and tests can pin a
//     fake scorer. cmd_align wires the real SpliceSiteDb.
//   * Lossless: every input sub-chain ends up in exactly one emitted
//     alignment (merged group or singleton — never silently dropped).
//   * Strand-aware: for '-' alignments the donor/acceptor intron motifs are
//     reverse-complemented and their roles swapped, so the scorer always sees
//     sense-strand GT/AG-class motifs.
//   * The classical chainer is NOT modified; DNA-mode is unaffected.

#pragma once

#include "mapping/chain_spliced.h"

#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace llmap::mapping {

// Junction scorer: given the 2-bp intron-side donor & acceptor motifs (already
// strand-normalised to the transcript sense) plus short intron flanks, return
// a junction probability in [0, 1]. Real implementation wraps
// annot::SpliceSiteDb::ScoreJunction.
using JunctionScorer =
    std::function<float(std::string_view donor2, std::string_view acceptor2,
                        std::string_view intron3p, std::string_view intron5p)>;

// Reference-sequence lookup by ref_id. Returns an empty view when unknown
// (junctions on that ref then fall back to the default 0.5 probability).
using RefSeqLookup = std::function<std::string_view(std::string_view ref_id)>;

struct TranscriptStageConfig {
    JoinerConfig joiner{};
};

// One alignment emitted after the spliced stage (a merged group or a
// surviving singleton).
struct SplicedAlignment {
    std::string ref_id;
    std::uint64_t ref_start{0};   // half-open span over all sub-chains
    std::uint64_t ref_end{0};
    std::uint32_t query_start{0};
    std::uint32_t query_end{0};
    char strand{'+'};
    std::int32_t score{0};
    std::string cigar;            // spliced (N ops) for merges, plain otherwise
    bool is_spliced{false};       // true when >1 sub-chain merged
    // genomic (donor_ref_pos, acceptor_ref_pos) per junction
    std::vector<std::pair<std::uint64_t, std::uint64_t>> junctions;
    std::vector<float> junction_conf;
};

// Run the spliced stage over ONE read's linear sub-chains. Groups by
// (ref_id, strand), sorts each group by ref_start, scores the inter-block
// gaps from reference motifs, joins, and emits SplicedAlignments.
[[nodiscard]] std::vector<SplicedAlignment> ApplyTranscriptStage(
    std::span<const LinearSubChain> sub_chains,
    const RefSeqLookup& ref_lookup,
    const JunctionScorer& score_junction,
    const TranscriptStageConfig& cfg = {});

// Strand-aware extraction of the sense-normalised donor/acceptor 2-bp motifs
// (and short intron flanks) for the gap between sub-chains `a` and `b` on the
// same reference `ref_seq`. `a.ref_end <= b.ref_start` is assumed (caller
// sorts). Returns false if positions fall outside `ref_seq`.
[[nodiscard]] bool ExtractJunctionMotifs(const LinearSubChain& a,
                                         const LinearSubChain& b,
                                         std::string_view ref_seq,
                                         std::string& donor2,
                                         std::string& acceptor2,
                                         std::string& intron3p,
                                         std::string& intron5p);

}  // namespace llmap::mapping
