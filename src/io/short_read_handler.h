// LLmap — short-read pipeline branch.
//
// Most public RNA-seq data is short-read (GTEx bulk, TCGA, Recount3,
// ENCODE, ENA / SRA). Long-read is more interesting biology but
// short-read is the volume. LLmap must handle both seriously; this
// module is the short-read-specific entry point.
//
// Key differences from the long-read path:
//
//   * k=21 default (vs k=51 long). Short reads are typically 75-300 bp;
//     a 51-mer barely fits and would lose most of the signal. 21-mer
//     hits get aggregated into equivalence classes downstream.
//
//   * Paired-end aware. Junction inference now uses split-read evidence
//     across mate pairs in addition to per-read CIGAR. A mate whose
//     mapped position implies an exon-spanning insert size flags a
//     candidate junction even when neither mate's CIGAR shows N.
//
//   * Equivalence classes. When per-read unique mapping is impossible
//     (sequence-identical paralogs, low MAPQ), we emit a Salmon/Kallisto
//     -style equivalence class: a set of transcripts the read could
//     belong to + the relative weight. Lossless: the read is preserved
//     as part of the class, never dropped.
//
// This module does NOT re-implement an aligner. It wraps the existing
// MinimizerIndex / classical pipeline and adapts their output into the
// short-read-friendly shape. The CLI flag `--mode short` activates this
// branch.

#pragma once

#include "anchor/anchor_record.h"
#include "core/transcript_kind.h"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::io {

// ===========================================================================
// Config.
// ===========================================================================

struct ShortReadConfig {
    bool paired_end{true};

    /// Expected fragment size (bp) — used to flag splice candidates
    /// where the implied insert grossly exceeds this value.
    std::uint16_t expected_insert_size{300};

    /// Short k for short reads.
    std::uint8_t kmer_short{21};
    std::uint8_t kmer_short_alt{31};

    /// If true, emit equivalence classes when unique mapping fails
    /// instead of dropping into Tentative.
    bool emit_equivalence_classes{true};

    /// Minimum confidence to keep an equivalence class member.
    float min_member_weight{0.05f};
};

// ===========================================================================
// Per-read input + output records.
// ===========================================================================

struct FastqRecordPair {
    std::string read_id;
    std::string seq_r1;
    std::string seq_r2;  // empty when single-end
    std::string qual_r1;
    std::string qual_r2;
};

/// One member of an equivalence class: a transcript the read could map
/// to, with a soft weight ∈ [0,1]. Sum of weights across a class
/// equals 1.0 (renormalised on emit).
struct EquivalenceClassMember {
    std::string transcript_id;
    float weight{0.0f};
};

/// Output for a single short read (or pair).
struct ShortReadOutcome {
    std::string read_id;
    std::vector<EquivalenceClassMember> members;  // empty ⇒ Unmapped

    /// Implied insert size when both mates mapped to the same transcript.
    /// 0 ⇒ unknown / single-end.
    std::uint32_t implied_insert_size{0};

    /// True if the implied insert exceeds 1.5× expected and the mates
    /// straddle an annotated exon-exon junction — flags a splicing event.
    bool junction_inferred{false};
};

// ===========================================================================
// Pipeline-level helpers.
// ===========================================================================

/// Build equivalence class from a vector of (transcript_id, score) hits.
/// Renormalises weights so they sum to 1.0; drops members below
/// cfg.min_member_weight.
std::vector<EquivalenceClassMember>
BuildEquivalenceClass(
    const std::vector<std::pair<std::string, float>>& hits,
    const ShortReadConfig& cfg);

/// Compute implied insert size from two mate positions on the same
/// transcript; returns 0 when mates don't share a transcript.
std::uint32_t ImpliedInsertSize(std::uint32_t pos1, std::uint32_t pos2);

/// Process one paired/single-end record. The caller has already done
/// the alignment work (passed in via `aln_hits` — transcript-id ×
/// raw score). This function compiles a ShortReadOutcome by applying
/// equivalence-class semantics + insert-size logic.
ShortReadOutcome BuildOutcome(
    const FastqRecordPair& rec,
    const std::vector<std::pair<std::string, float>>& aln_hits_r1,
    const std::vector<std::pair<std::string, float>>& aln_hits_r2,
    std::optional<std::uint32_t> implied_insert,
    bool junction_seen,
    const ShortReadConfig& cfg);

}  // namespace llmap::io
