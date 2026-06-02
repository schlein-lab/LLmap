// LLmap — ChimeraDetector + VdjLocusMask implementation.

#include "chimera/chimera_detector.h"

#include <algorithm>
#include <cmath>

namespace llmap::chimera {

// ===========================================================================
// VdjLocusMask
// ===========================================================================

void VdjLocusMask::LoadGrch38Defaults() {
    // Ensembl 110 / GRCh38 coordinates.
    Add({"chr14", 105'586'437, 106'879'844, "IGH"});
    Add({"chr2",   88'857'160,  90'238'368, "IGK"});
    Add({"chr22", 22'026'076,  22'922'913, "IGL"});
    Add({"chr14", 21'621'904,  22'552'132, "TRA"});
    Add({"chr14", 22'422'546,  22'550'064, "TRD"});  // nested in TRA
    Add({"chr7",  142'299'011, 142'813'287, "TRB"});
    Add({"chr7",  38'240'024,  38'368'055,  "TRG"});
}

void VdjLocusMask::Add(VdjRegion r) {
    regions_.push_back(std::move(r));
}

std::string VdjLocusMask::LocusAt(std::string_view chrom,
                                    std::uint64_t pos) const {
    for (const auto& r : regions_) {
        if (r.chrom == chrom && pos >= r.start && pos < r.end) {
            return r.locus_id;
        }
    }
    return {};
}

bool VdjLocusMask::BothInSameVdjLocus(
    std::string_view chrom_a, std::uint64_t pos_a,
    std::string_view chrom_b, std::uint64_t pos_b) const {
    const std::string la = LocusAt(chrom_a, pos_a);
    const std::string lb = LocusAt(chrom_b, pos_b);
    return !la.empty() && la == lb;
}

// ===========================================================================
// DistancePrior — Plan-Block 7 table.
//
// log-prior approximations (relative scale, used for tie-breaking only):
//   intra-IGH switch-context     log p ≈ -0.7   (most common chimerism)
//   intra-chrom < 1 Mb            log p ≈ -2.0
//   intra-chrom 1 Mb..50 Mb       log p ≈ -4.0  (distance penalty)
//   intra-chrom > 50 Mb           log p ≈ -6.0
//   cross-chrom                   log p ≈ -8.0  (translocations rare)
// ===========================================================================

float DistancePrior(std::uint64_t genomic_distance_bp,
                     bool cross_chrom) {
    if (cross_chrom) return -8.0f;
    if (genomic_distance_bp <= 1'000'000)         return -2.0f;
    if (genomic_distance_bp <= 50'000'000ULL)     return -4.0f;
    return -6.0f;
}

// ===========================================================================
// Helper utilities for Analyze().
// ===========================================================================

namespace {

float PartCoverage(const AlignedPart& p, std::uint32_t read_len) {
    if (read_len == 0) return 0.0f;
    return static_cast<float>(p.read_length) / static_cast<float>(read_len);
}

bool AnyPartIsSwitchRegion(std::span<const AlignedPart> parts) {
    for (const auto& p : parts) {
        if (p.is_switch_region) return true;
    }
    return false;
}

}  // namespace

// ===========================================================================
// Analyze
// ===========================================================================

std::vector<ChimericHypothesis>
Analyze(std::span<const AlignedPart> parts,
        const VdjLocusMask& vdj_mask,
        const ChimeraConfig& cfg) {

    std::vector<ChimericHypothesis> out;
    if (parts.size() < 2) return out;  // not chimeric

    // Estimate the read length as the max read_offset + read_length.
    std::uint32_t total_read_len = 0;
    for (const auto& p : parts) {
        const std::uint32_t end = p.read_offset + p.read_length;
        if (end > total_read_len) total_read_len = end;
    }
    if (total_read_len == 0) return out;

    // Filter parts below min_part_coverage.
    std::vector<AlignedPart> kept;
    for (const auto& p : parts) {
        if (PartCoverage(p, total_read_len) >= cfg.min_part_coverage) {
            kept.push_back(p);
        }
    }
    if (kept.size() < 2) return out;

    // Combined coverage must reach min_combined_coverage to even
    // consider the read chimeric.
    float combined = 0.0f;
    for (const auto& p : kept) combined += PartCoverage(p, total_read_len);
    if (combined < cfg.min_combined_coverage) return out;

    // Sort by read_offset so the parts come in 5'→3' order.
    std::sort(kept.begin(), kept.end(),
              [](const AlignedPart& a, const AlignedPart& b) {
                  return a.read_offset < b.read_offset;
              });

    // Establish chrom-uniformity.
    bool same_chrom = true;
    for (std::size_t i = 1; i < kept.size(); ++i) {
        if (kept[i].ref_chrom != kept[0].ref_chrom) {
            same_chrom = false;
            break;
        }
    }

    // Distance metric (only meaningful when same_chrom).
    std::uint64_t span_bp = 0;
    if (same_chrom) {
        std::uint64_t lo = kept.front().ref_start;
        std::uint64_t hi = kept.front().ref_end;
        for (const auto& p : kept) {
            if (p.ref_start < lo) lo = p.ref_start;
            if (p.ref_end   > hi) hi = p.ref_end;
        }
        if (hi > lo) span_bp = hi - lo;
    }

    ChimericHypothesis h;
    h.parts = kept;
    h.part_probabilities.resize(kept.size(), 1.0f / kept.size());
    h.genomic_distance_bp = same_chrom ? span_bp : 0;
    h.vdj_class_switch_detected = AnyPartIsSwitchRegion(kept);

    // ----- Kind classification -----
    if (!same_chrom) {
        h.kind = 'X';
        h.log_prior = DistancePrior(0, /*cross_chrom=*/true);
    } else {
        const bool in_vdj = vdj_mask.BothInSameVdjLocus(
            kept.front().ref_chrom, kept.front().ref_start,
            kept.back().ref_chrom,  kept.back().ref_end);
        if (in_vdj && h.vdj_class_switch_detected) {
            h.kind = 'V';
            h.log_prior = -0.7f;  // VDJ class-switch is biologically common
        } else if (span_bp <= cfg.intra_region_max_bp) {
            h.kind = 'I';
            h.log_prior = DistancePrior(span_bp, false);
        } else {
            h.kind = 'I';
            h.log_prior = DistancePrior(span_bp, false);
        }
    }

    // Emit unless cross-chrom is suppressed.
    if (h.kind != 'X' || cfg.emit_low_prior_chimeric) {
        out.push_back(std::move(h));
    }

    return out;
}

}  // namespace llmap::chimera
