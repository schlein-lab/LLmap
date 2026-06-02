// LLmap — Chimera detector + VDJ-recombination mask.
//
// Reads that align to two (or more) distinct loci in the same read get
// classified into ChimericIntraRegion / ChimericInterChrom /
// (VDJ class-switch sub-type of IntraRegion). Distance-aware priors
// drive the discrimination:
//
//   intra-IGH (anchor parts both in IGH locus, possibly with a
//   detected switch-region between them) → strong prior, ChimericIntraRegion
//   with kind='V' (VDJ recombination / class-switch)
//
//   intra-chrom non-IG (same chromosome, < 1 Mb apart, no switch
//   evidence) → ChimericIntraRegion with kind='I'
//
//   cross-chromosomal → ChimericInterChrom with kind='X'
//
// Per Plan-Block 7 the chimera detector is the chimera-side input
// generator to the Lossless-Schema. The actual AlignmentRecord
// construction lives elsewhere (Block 8 BAM/Parquet writers); we just
// produce ChimericHypothesis structs that the caller turns into
// AlignmentRecord::chimeric.
//
// VDJ mask: simple chromosomal range table for IGH (chr14), IGK
// (chr2p), IGL (chr22), TRA (chr14), TRB (chr7), TRG (chr7), TRD
// (chr14). When both anchor parts of a chimeric pair fall in the same
// VDJ region AND the configured switch-region tag is present in the
// catalog/anchor metadata, we set kind='V'.

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace llmap::chimera {

// ===========================================================================
// VdjLocusMask
// ===========================================================================

struct VdjRegion {
    std::string chrom;
    std::uint64_t start{0};   // 0-based inclusive
    std::uint64_t end{0};     // exclusive
    std::string locus_id;     // "IGH" / "IGK" / "IGL" / "TRA" / ...
};

class VdjLocusMask {
public:
    /// Built-in mask for human GRCh38 IGH/IGK/IGL/TR* loci. Always
    /// succeeds; numbers are from Ensembl 110.
    void LoadGrch38Defaults();

    /// Add a custom VDJ-locus region. Useful for non-human assemblies
    /// or tests.
    void Add(VdjRegion r);

    /// Returns the VDJ locus id containing (chrom, pos), or empty
    /// string when none.
    [[nodiscard]] std::string LocusAt(std::string_view chrom,
                                       std::uint64_t pos) const;

    /// True iff the two positions sit inside the SAME VDJ locus.
    [[nodiscard]] bool BothInSameVdjLocus(
        std::string_view chrom_a, std::uint64_t pos_a,
        std::string_view chrom_b, std::uint64_t pos_b) const;

    [[nodiscard]] std::size_t size() const noexcept { return regions_.size(); }

private:
    std::vector<VdjRegion> regions_;
};

// ===========================================================================
// Chimeric inputs.
//
// AlignedPart mirrors the relevant bits of AlignmentHit so this
// module doesn't drag in the AlignmentRecord header.
// ===========================================================================

struct AlignedPart {
    std::string ref_chrom;
    std::uint64_t ref_start{0};   // 0-based inclusive
    std::uint64_t ref_end{0};
    std::uint32_t read_offset{0};
    std::uint32_t read_length{0};
    float score{0.0f};
    /// Whether the part falls within a known switch region (Sγ/Sμ/Sα/…).
    bool is_switch_region{false};
};

// ===========================================================================
// Output — one ChimericHypothesis.
// ===========================================================================

struct ChimericHypothesis {
    std::vector<AlignedPart> parts;          // ordered along read
    std::vector<float> part_probabilities;
    /// 'I' intra-region / 'X' cross-chromosomal / 'V' VDJ class-switch
    char kind{'.'};
    /// True iff a switch region sits between two of the parts (=
    /// hallmark of class-switch recombination).
    bool vdj_class_switch_detected{false};
    /// Implied genomic distance between the most distant parts (bp).
    /// Set only when all parts share the same chrom.
    std::uint64_t genomic_distance_bp{0};
    /// Combined log-prior — kept on a log scale so callers can compare
    /// chimeric vs non-chimeric without log-sum-exp gymnastics.
    float log_prior{0.0f};
};

// ===========================================================================
// Detector configuration.
// ===========================================================================

struct ChimeraConfig {
    /// Minimum coverage (fraction of read length) each chimeric part
    /// must contribute. Parts below this are ignored.
    float min_part_coverage{0.10f};

    /// Combined coverage threshold below which the read should NOT be
    /// flagged chimeric (it's a partial match instead).
    float min_combined_coverage{0.90f};

    /// Emit ChimericInterChrom hypotheses even at low prior. Plan-
    /// Block 7 requires this be true by default — lossless preserves
    /// the read, downstream consumers can filter.
    bool emit_low_prior_chimeric{true};

    /// Genomic distance below which kind='I' (intra-region) is emitted.
    /// Above this (but same chrom) the prior decays toward inter-chrom.
    std::uint64_t intra_region_max_bp{1'000'000};
};

// ===========================================================================
// Public API.
// ===========================================================================

/// Analyse a set of aligned parts (from one read) and emit chimeric
/// hypotheses. Returns at most ONE hypothesis per input combination —
/// the detector is purely a classifier on the already-aligned parts,
/// not an alignment pass.
///
/// Caller composes the final AlignmentRecord::chimeric from the
/// hypothesis (Block 8).
[[nodiscard]] std::vector<ChimericHypothesis>
Analyze(std::span<const AlignedPart> parts,
        const VdjLocusMask& vdj_mask,
        const ChimeraConfig& cfg = {});

/// Convenience: distance-aware prior. Returns a log-prior used inside
/// Analyze() and exposed here for tests + Wave-Collapse weighting.
[[nodiscard]] float DistancePrior(std::uint64_t genomic_distance_bp,
                                   bool cross_chrom);

}  // namespace llmap::chimera
