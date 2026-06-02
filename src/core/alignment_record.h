// LLmap — lossless alignment record.
//
// The lossless guarantee is enforced at construction: every record must
// have a defined AlignmentStatus and the optional fields required by
// that status must be populated. `make_*` factory functions are the
// only sanctioned construction path; they return records whose
// is_lossless_consistent() returns true at the point of construction.
//
// This file used to carry a 3-status model
//   {Mapped, Tentative, Unmapped}
// — which forced every spliced / chimeric / partial / dark-transcriptome
// observation into one of those three buckets, losing biological
// information. The Transcript-Mode plan (~/.claude/plans/
// memoized-crafting-lecun.md, Block 1) expands the status set to 10 to
// make these distinct cases first-class:
//
//   Mapped                — primary converged ≥ 0.99 (standard mRNA)
//   MappedSterile         — sterile germline transcript (Iγ→Sγ→Cγ, no V-D-J)
//   MappedPreMrna         — pre-mRNA / retained-intron transcript
//   MappedCircular        — circRNA back-splice junction detected
//   PartialMappedExact    — only part of read 100%-matches a known anchor;
//                            the rest is preserved verbatim in the record
//   ChimericIntraRegion   — two or more anchors within the same locus
//                            (e.g. VDJ class-switch IGHM-CH1 → IGHG4-CH1)
//   ChimericInterChrom    — anchors on different chromosomes
//                            (rare; translocations; preserved, not dropped)
//   DarkNovel             — read forms or joins a cluster but no DB anchor
//                            matches; cluster-internal anchor is recorded
//   Tentative             — EM did not converge but Top-K candidates exist
//   Unmapped              — nothing — explicit RejectionReason required
//
// All ten branches share is_lossless_consistent() invariants documented
// at the bottom of this header.

#pragma once

#include "core/transcript_kind.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>
#include <stdexcept>

namespace llmap {

// ===========================================================================
// Enums
// ===========================================================================

enum class AlignmentStatus : std::uint8_t {
    Mapped = 0,
    MappedSterile,           // sterile germline transcript
    MappedPreMrna,           // pre-mRNA with retained introns
    MappedCircular,          // circRNA back-splice
    PartialMappedExact,      // partial-but-exact (100%) match
    ChimericIntraRegion,     // multi-anchor within a region
    ChimericInterChrom,      // multi-anchor across chromosomes
    DarkNovel,               // cluster-anchor only, no DB hit
    Tentative,               // EM didn't converge, distribution preserved
    Unmapped,                // explicit reject; carries RejectionReason
};

/// True for every status that represents a successful alignment of some
/// kind, regardless of which Mapped* variant it is. Useful in output
/// writers that need a quick yes/no on "is this record alignable".
constexpr bool IsAnyMapped(AlignmentStatus s) noexcept {
    switch (s) {
        case AlignmentStatus::Mapped:
        case AlignmentStatus::MappedSterile:
        case AlignmentStatus::MappedPreMrna:
        case AlignmentStatus::MappedCircular:
        case AlignmentStatus::PartialMappedExact:
        case AlignmentStatus::ChimericIntraRegion:
        case AlignmentStatus::ChimericInterChrom:
        case AlignmentStatus::DarkNovel:
            return true;
        case AlignmentStatus::Tentative:
        case AlignmentStatus::Unmapped:
            return false;
    }
    return false;
}

enum class RejectionReason : std::uint8_t {
    NoSeeds,
    LowSeedDensity,
    ChainScoreBelowThreshold,
    AmbiguousNoAnchor,
    FailedWfa2Extension,
    DidNotConverge,
    HostContamination,
    LowComplexity,
};

/// Stable string for telemetry / output (lossless.summary.json).
const char* AlignmentStatusName(AlignmentStatus s) noexcept;
const char* RejectionReasonName(RejectionReason r) noexcept;

// ===========================================================================
// Per-hit records — shared by several status branches.
// ===========================================================================

struct CigarString {
    std::string ops;
};

struct PsvObservation {
    std::uint64_t psv_id;
    char allele;
    float quality;
};

struct AlignmentHit {
    std::string target_id;
    std::uint64_t start{0};
    std::uint64_t end{0};
    CigarString cigar;
    std::int32_t score{0};
    std::uint32_t nm{0};
    bool is_reverse{false};
    std::vector<PsvObservation> psv_calls;
};

struct TentativeTarget {
    std::string target_id;
    std::uint64_t approx_start{0};
    std::uint64_t approx_end{0};
    std::uint32_t n_seeds{0};
    std::int32_t partial_chain_score{0};
    float sequence_identity_estimate{0.0f};
    float final_probability{0.0f};
};

struct AnchorEvidence {
    bool five_prime_anchored{false};
    bool three_prime_anchored{false};
    float anchor_confidence{0.0f};
};

struct ParalogCall {
    // inter_paralog: P(paralog_k) keyed by paralog id; vector ordered for
    // cache friendliness on small (typically ≤ 6) paralog sets.
    std::vector<std::pair<std::string, float>> inter_paralog;
    std::optional<float> p_canonical;
    std::optional<float> p_dup;
    std::uint32_t n_discriminating_psvs{0};
};

// ===========================================================================
// Status-specific payloads.
//
// Each new status carries its own minimal struct; pulling them out keeps
// the AlignmentRecord clean and lets the invariant checker reason about
// "this optional must be present for status X" purely from the type.
// ===========================================================================

/// Payload for ChimericIntraRegion / ChimericInterChrom.
///
/// Each Part is an AlignmentHit on a separate anchor; combined they
/// account for the full read (or a labelled fraction of it). `kind`
/// classifies the chimeric flavour:
///   'I' — intra-region (anchors within ~1 Mb on same chrom)
///   'X' — cross-chromosomal
///   'V' — VDJ class-switch / recombination (IG/TR loci); always intra
struct ChimericDetail {
    std::vector<AlignmentHit> parts;
    std::vector<float> part_probabilities;   // Wave-Collapse mass per part
    char kind{'.'};
    bool vdj_class_switch_detected{false};   // true ⇒ kind=='V'
    std::optional<std::uint64_t> genomic_distance_bp;  // only for intra-chrom
};

/// Payload for PartialMappedExact: the (start,end) of the matched
/// substring inside the read, the anchor id, and the SHA-1 of the
/// unmatched tail so downstream tools can detect duplicates without
/// re-storing the sequence.
struct PartialMatchDetail {
    std::uint32_t read_offset{0};            // 0-based start in read
    std::uint32_t match_length{0};
    std::string anchor_id;
    std::string unmatched_tail;              // raw seq of the rest
};

/// Payload for DarkNovel: the cluster-internal anchor id (computed in
/// Stage-1 Self-WaveCollapse) plus the cluster id and the size of the
/// cluster the read joined. No DB anchor matched — this is by design.
struct DarkNovelDetail {
    std::uint32_t cluster_id{0};
    std::string cluster_anchor_id;           // computed anchor id
    std::uint32_t cluster_size{0};           // number of reads in cluster
    std::vector<std::string> cluster_anchor_ids;  // peer anchors
};

/// Payload for MappedSterile: which I-promoter / S-region / C-gene the
/// read crossed. All three identifiers needed for a confident sterile-
/// germline call; missing any of them ⇒ MappedSterile invariant breaks.
struct SterileDetail {
    std::string i_promoter_id;               // e.g. "I_gamma4_promoter"
    std::string s_region_id;                 // e.g. "S_gamma4"
    std::string c_gene_id;                   // e.g. "IGHG4"
};

/// Payload for MappedPreMrna: per-read intron-retention map.
/// retained_introns is a vector of (read_start, read_end) ranges in the
/// READ coordinate (not genomic) for each retained intron.
struct PreMrnaDetail {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> retained_introns;
    float pct_introns_retained{0.0f};        // 0..1 retained / total
};

/// Payload for MappedCircular: back-splice junction position in genomic
/// coordinates, host gene, and (optionally) the host's transcript id.
struct CircularDetail {
    std::uint64_t backsplice_donor_pos{0};
    std::uint64_t backsplice_acceptor_pos{0};   // < donor in linear genome
    std::string host_gene_id;
    std::optional<std::string> host_transcript_id;
};

/// AID-mediated C→U editing footprint (IG locus specific).
///
/// Returned separately from the Mapped* status so it can ride along on
/// any successful alignment that happens to fall in a switch-region.
/// See [[ighg4_sgamma4_identical_in_tandem_dup]].
struct AidFootprint {
    std::uint32_t n_c_to_u_events{0};
    std::vector<std::uint32_t> c_to_u_positions_in_read;
    std::string switch_region_id;            // e.g. "S_gamma4"
};

// ===========================================================================
// The lossless record.
// ===========================================================================

struct AlignmentRecord {
    // --- universal --------------------------------------------------------
    std::string read_id;
    std::uint32_t read_len{0};
    AlignmentStatus status{AlignmentStatus::Unmapped};

    /// TranscriptKind annotation. Filled for any Mapped* status from the
    /// matched AnchorRecord; left as Unknown for Tentative/Unmapped.
    core::TranscriptKind transcript_kind{core::TranscriptKind::Unknown};

    /// Free-form extension label when transcript_kind == NovelUnclassified.
    std::optional<core::CustomKindTag> custom_kind;

    // --- Mapped (canonical) ----------------------------------------------
    std::optional<AlignmentHit> primary;
    std::vector<AlignmentHit> alternatives;

    // --- Tentative -------------------------------------------------------
    std::vector<TentativeTarget> tentative_targets;
    std::vector<float> confidence_scores;

    // --- shared anchor + paralog evidence --------------------------------
    AnchorEvidence anchor_resolution;
    std::optional<ParalogCall> paralog_assignment;

    // --- status-specific payloads (only one populated per record) --------
    std::optional<ChimericDetail>    chimeric;        // ChimericIntra/Inter
    std::optional<PartialMatchDetail> partial_match;  // PartialMappedExact
    std::optional<DarkNovelDetail>   dark_novel;      // DarkNovel
    std::optional<SterileDetail>     sterile;         // MappedSterile
    std::optional<PreMrnaDetail>     pre_mrna;        // MappedPreMrna
    std::optional<CircularDetail>    circular;        // MappedCircular

    // --- orthogonal annotations (can ride on any Mapped* status) ---------
    std::optional<AidFootprint> aid_footprint;

    // --- single-cell / phasing -------------------------------------------
    std::optional<std::string> cell_barcode;
    std::optional<std::string> umi;
    std::optional<std::uint8_t> haplotype;

    // --- rejection -------------------------------------------------------
    std::optional<RejectionReason> rejection_reason;

    // --- WaveCollapse provenance -----------------------------------------
    std::uint32_t collapsed_at_iteration{0};
    std::uint8_t collapsed_at_level{0};
    std::uint32_t cluster_id{0};
    bool is_cluster_representative{false};

    // ---------------------------------------------------------------------
    // Lossless invariant
    //
    //   Mapped              ⇒ primary.has_value()
    //   MappedSterile       ⇒ primary.has_value() && sterile.has_value()
    //   MappedPreMrna       ⇒ primary.has_value() && pre_mrna.has_value()
    //   MappedCircular      ⇒ primary.has_value() && circular.has_value()
    //   PartialMappedExact  ⇒ partial_match.has_value()
    //                          && partial_match->match_length > 0
    //   ChimericIntraRegion ⇒ chimeric.has_value()
    //                          && chimeric->parts.size() >= 2
    //                          && chimeric->kind ∈ {'I','V'}
    //   ChimericInterChrom  ⇒ chimeric.has_value()
    //                          && chimeric->parts.size() >= 2
    //                          && chimeric->kind == 'X'
    //   DarkNovel           ⇒ dark_novel.has_value()
    //                          && !dark_novel->cluster_anchor_id.empty()
    //   Tentative           ⇒ !tentative_targets.empty()
    //                          && rejection_reason == DidNotConverge
    //   Unmapped            ⇒ rejection_reason.has_value()
    //
    // The check is deliberately one function rather than per-status type
    // safety because AlignmentRecord must remain trivially-copyable
    // (it's the single output unit; we don't want to enforce class
    // hierarchies on it).
    // ---------------------------------------------------------------------
    [[nodiscard]] bool is_lossless_consistent() const noexcept;
};

// ===========================================================================
// Factory functions — the only sanctioned construction path.
//
// Every factory ends with assert(rec.is_lossless_consistent()) in debug
// builds; in release builds the invariant relies on the caller, but the
// factories are designed so that "give me what I need" mechanically
// satisfies the invariant.
// ===========================================================================

[[nodiscard]] AlignmentRecord make_mapped(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    std::vector<AlignmentHit> alternatives = {},
    core::TranscriptKind kind = core::TranscriptKind::MatureMrna
);

[[nodiscard]] AlignmentRecord make_mapped_sterile(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    SterileDetail sterile,
    std::vector<AlignmentHit> alternatives = {}
);

[[nodiscard]] AlignmentRecord make_mapped_pre_mrna(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    PreMrnaDetail pre_mrna,
    std::vector<AlignmentHit> alternatives = {}
);

[[nodiscard]] AlignmentRecord make_mapped_circular(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    CircularDetail circular
);

[[nodiscard]] AlignmentRecord make_partial_exact(
    std::string read_id,
    std::uint32_t read_len,
    PartialMatchDetail detail
);

[[nodiscard]] AlignmentRecord make_chimeric(
    std::string read_id,
    std::uint32_t read_len,
    ChimericDetail detail
);

[[nodiscard]] AlignmentRecord make_dark_novel(
    std::string read_id,
    std::uint32_t read_len,
    DarkNovelDetail detail
);

[[nodiscard]] AlignmentRecord make_tentative(
    std::string read_id,
    std::uint32_t read_len,
    std::vector<TentativeTarget> targets,
    RejectionReason reason
);

[[nodiscard]] AlignmentRecord make_unmapped(
    std::string read_id,
    std::uint32_t read_len,
    RejectionReason reason
);

}  // namespace llmap
