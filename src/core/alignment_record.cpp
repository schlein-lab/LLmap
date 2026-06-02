// LLmap — AlignmentRecord factories + lossless invariant.
//
// One source-of-truth per status branch. Each factory mechanically
// populates the fields the lossless invariant requires; in debug builds
// the assertion below catches accidental drift.
//
// File deliberately kept under the 400-LOC modular cap by delegating
// status/rejection-name tables to small lookup loops.

#include "core/alignment_record.h"

#include <cassert>

namespace llmap {

// ===========================================================================
// Stable strings — never rename existing entries (output-format contract).
// ===========================================================================

const char* AlignmentStatusName(AlignmentStatus s) noexcept {
    switch (s) {
        case AlignmentStatus::Mapped:              return "MAPPED";
        case AlignmentStatus::MappedSterile:       return "MAPPED_STERILE";
        case AlignmentStatus::MappedPreMrna:       return "MAPPED_PREMRNA";
        case AlignmentStatus::MappedCircular:      return "MAPPED_CIRCULAR";
        case AlignmentStatus::PartialMappedExact:  return "PARTIAL";
        case AlignmentStatus::ChimericIntraRegion: return "CHIM_INTRA";
        case AlignmentStatus::ChimericInterChrom:  return "CHIM_INTER";
        case AlignmentStatus::DarkNovel:           return "DARK";
        case AlignmentStatus::Tentative:           return "TENT";
        case AlignmentStatus::Unmapped:            return "UNMAPPED";
    }
    return "UNKNOWN";
}

const char* RejectionReasonName(RejectionReason r) noexcept {
    switch (r) {
        case RejectionReason::NoSeeds:                  return "NoSeeds";
        case RejectionReason::LowSeedDensity:           return "LowSeedDensity";
        case RejectionReason::ChainScoreBelowThreshold: return "ChainScoreBelowThreshold";
        case RejectionReason::AmbiguousNoAnchor:        return "AmbiguousNoAnchor";
        case RejectionReason::FailedWfa2Extension:      return "FailedWfa2Extension";
        case RejectionReason::DidNotConverge:           return "DidNotConverge";
        case RejectionReason::HostContamination:        return "HostContamination";
        case RejectionReason::LowComplexity:            return "LowComplexity";
    }
    return "Unknown";
}

// ===========================================================================
// is_lossless_consistent
//
// One switch per status, each branch encoding the invariant documented
// at the bottom of alignment_record.h. Mismatch ⇒ false.
// ===========================================================================

bool AlignmentRecord::is_lossless_consistent() const noexcept {
    if (read_id.empty()) return false;

    switch (status) {
        case AlignmentStatus::Mapped:
            return primary.has_value()
                && !rejection_reason.has_value();

        case AlignmentStatus::MappedSterile:
            return primary.has_value()
                && sterile.has_value()
                && !sterile->i_promoter_id.empty()
                && !sterile->s_region_id.empty()
                && !sterile->c_gene_id.empty()
                && !rejection_reason.has_value();

        case AlignmentStatus::MappedPreMrna:
            return primary.has_value()
                && pre_mrna.has_value()
                // retained_introns may be empty if pct=0 (edge case),
                // but pct must be in [0,1]
                && pre_mrna->pct_introns_retained >= 0.0f
                && pre_mrna->pct_introns_retained <= 1.0f
                && !rejection_reason.has_value();

        case AlignmentStatus::MappedCircular:
            return primary.has_value()
                && circular.has_value()
                // back-splice acceptor must come before donor in linear
                // genomic coords (this is what makes it a circle)
                && circular->backsplice_acceptor_pos
                    < circular->backsplice_donor_pos
                && !circular->host_gene_id.empty()
                && !rejection_reason.has_value();

        case AlignmentStatus::PartialMappedExact:
            return partial_match.has_value()
                && partial_match->match_length > 0
                && !partial_match->anchor_id.empty()
                && !rejection_reason.has_value();

        case AlignmentStatus::ChimericIntraRegion:
            return chimeric.has_value()
                && chimeric->parts.size() >= 2
                && chimeric->part_probabilities.size()
                    == chimeric->parts.size()
                && (chimeric->kind == 'I' || chimeric->kind == 'V')
                && !rejection_reason.has_value();

        case AlignmentStatus::ChimericInterChrom:
            return chimeric.has_value()
                && chimeric->parts.size() >= 2
                && chimeric->part_probabilities.size()
                    == chimeric->parts.size()
                && chimeric->kind == 'X'
                && !rejection_reason.has_value();

        case AlignmentStatus::DarkNovel:
            return dark_novel.has_value()
                && !dark_novel->cluster_anchor_id.empty()
                && dark_novel->cluster_size > 0
                && !rejection_reason.has_value();

        case AlignmentStatus::Tentative:
            return !primary.has_value()
                && !tentative_targets.empty()
                && rejection_reason.has_value();

        case AlignmentStatus::Unmapped:
            return !primary.has_value()
                && rejection_reason.has_value();
    }
    return false;
}

// ===========================================================================
// Factories
// ===========================================================================

AlignmentRecord make_mapped(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    std::vector<AlignmentHit> alternatives,
    core::TranscriptKind kind
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::Mapped;
    r.transcript_kind = kind;
    r.primary = std::move(primary);
    r.alternatives = std::move(alternatives);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_mapped_sterile(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    SterileDetail sterile,
    std::vector<AlignmentHit> alternatives
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::MappedSterile;
    r.transcript_kind = core::TranscriptKind::SterileGermline;
    r.primary = std::move(primary);
    r.alternatives = std::move(alternatives);
    r.sterile = std::move(sterile);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_mapped_pre_mrna(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    PreMrnaDetail pre_mrna,
    std::vector<AlignmentHit> alternatives
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::MappedPreMrna;
    r.transcript_kind = core::TranscriptKind::PreMrna;
    r.primary = std::move(primary);
    r.alternatives = std::move(alternatives);
    r.pre_mrna = std::move(pre_mrna);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_mapped_circular(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    CircularDetail circular
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::MappedCircular;
    r.transcript_kind = core::TranscriptKind::CircularRna;
    r.primary = std::move(primary);
    r.circular = std::move(circular);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_partial_exact(
    std::string read_id,
    std::uint32_t read_len,
    PartialMatchDetail detail
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::PartialMappedExact;
    r.partial_match = std::move(detail);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_chimeric(
    std::string read_id,
    std::uint32_t read_len,
    ChimericDetail detail
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = (detail.kind == 'X')
        ? AlignmentStatus::ChimericInterChrom
        : AlignmentStatus::ChimericIntraRegion;
    r.chimeric = std::move(detail);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_dark_novel(
    std::string read_id,
    std::uint32_t read_len,
    DarkNovelDetail detail
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::DarkNovel;
    r.dark_novel = std::move(detail);
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_tentative(
    std::string read_id,
    std::uint32_t read_len,
    std::vector<TentativeTarget> targets,
    RejectionReason reason
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::Tentative;
    r.tentative_targets = std::move(targets);
    r.rejection_reason = reason;
    assert(r.is_lossless_consistent());
    return r;
}

AlignmentRecord make_unmapped(
    std::string read_id,
    std::uint32_t read_len,
    RejectionReason reason
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::Unmapped;
    r.rejection_reason = reason;
    assert(r.is_lossless_consistent());
    return r;
}

}  // namespace llmap
