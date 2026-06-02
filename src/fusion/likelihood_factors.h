// LLmap — Multi-Signal Fusion Engine: the differentiating layer.
//
// This module is where LLmap takes the categorical jump past minimap2 /
// IsoQuant / FLAIR / Bambu. Every existing long-read aligner uses each
// signal in isolation (sequence-only seed-extend, MAPQ threshold filter,
// post-hoc expression normalisation). LLmap fuses 10 signals into a
// single multiplicative likelihood that Wave-Collapse propagates:
//
//   L(r | b) = ∏ L_i(r | b),  i ∈ {1..10}
//
// Each L_i ∈ [floor, 1] (floor = 1e-6 to preserve lossless semantics —
// no signal may zero out a candidate; rare biological realities exist).
//
// The 10 factors, per the Block-4.5 design in
// ~/.claude/plans/memoized-crafting-lecun.md:
//
//   1. L_sequence              Forward pairHMM with modification tolerance
//   2. L_modification          REPIC m6A / REDIportal A→I match
//   3. L_depth_coverage        per-exon Negative-Binomial
//   4. L_expression_prior      tissue × cell-type TPM sigmoid
//   5. L_phasing               HP-tag consistency
//   6. L_pseudogene_compatibility  biotype-aware
//   7. L_junction              Splice-site PWM × DB evidence
//   8. L_barcode_context       single-cell cell-type prior
//   9. L_mapq_signal           sigmoid(MAPQ - 10) / 5  — NEVER a filter
//  10. L_length_plausibility   per TranscriptKind length window
//
// All numerical defaults are literature-calibrated (Pratanwanich 2021,
// Robinson 2010, Burge & Sharp 1997, LRGASP 2024 etc.). They are
// overridable per run via a Config file or the CLI.

#pragma once

#include "anchor/anchor_record.h"
#include "core/transcript_kind.h"
#include "fusion/expression_db.h"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace llmap::fusion {

// ===========================================================================
// Per-factor lossless floor — no factor may zero out a candidate.
// 1e-6 is small enough that a single failing factor still gets dominated
// by anything else passing, but large enough that 10 simultaneous floors
// stay above 1e-60 (representable in float32).
// ===========================================================================
inline constexpr float kFloor = 1.0e-6f;

// ===========================================================================
// LikelihoodFactors — one struct per (read, anchor) pair.
// ===========================================================================

struct LikelihoodFactors {
    float L_sequence{1.0f};
    float L_modification{1.0f};
    float L_depth_coverage{1.0f};
    float L_expression_prior{1.0f};
    float L_phasing{1.0f};
    float L_pseudogene_compatibility{1.0f};
    float L_junction{1.0f};
    float L_barcode_context{1.0f};
    float L_mapq_signal{1.0f};
    float L_length_plausibility{1.0f};

    /// Multiplicative product with lossless floor. Returns a value in
    /// [kFloor, 1].
    [[nodiscard]] float Product() const noexcept;
};

// ===========================================================================
// Disable mask — CLI exposure to turn factors off (set L_x = 1.0).
// ===========================================================================

struct FactorDisableMask {
    bool sequence{false};
    bool modification{false};
    bool depth_coverage{false};
    bool expression_prior{false};
    bool phasing{false};
    bool pseudogene_compatibility{false};
    bool junction{false};
    bool barcode_context{false};
    bool mapq_signal{false};
    bool length_plausibility{false};
};

// ===========================================================================
// Inputs.
//
// We pass them as discrete optional fields rather than a giant Context
// struct so callers can leave any input nullopt / empty when they don't
// have it — and the corresponding factor degrades gracefully to 1.0
// (neutral). This is *the* mechanism for "no signal = no penalty";
// minimap2's equivalent (MAPQ < threshold → drop) is what we explicitly
// reject.
// ===========================================================================

struct ReadContext {
    std::string_view read_id;
    std::uint32_t read_length{0};
    std::string_view read_sequence;
    /// Sequencing platform: affects sequence PWM and base-quality model.
    /// One of: "hifi" / "isoseq" / "ont_rna004" / "ont_legacy" /
    /// "illumina" / "unknown". Unknown ⇒ defaults to "ont_rna004".
    std::string_view platform{"ont_rna004"};
    /// MAPQ from upstream seed-extend stage; 0..60. -1 ⇒ unknown.
    std::int32_t mapq{-1};
    /// Haplotype (0/1) if upstream phasing known, else nullopt.
    std::optional<std::uint8_t> haplotype;
    /// Cell-type label if known from single-cell barcode demux.
    std::string_view cell_type;
};

struct TissueContext {
    std::string_view label;   ///< "lymph" / "pbmc" / "b_cell" / "" (unknown)
};

/// Observed RNA-modification calls on this read at this position window.
/// Index 0 = position pos_in_read=0; index k = a possible mod at
/// pos_in_read=k. nullopt ⇒ no call from upstream tool at that pos.
/// Vector length matches read length when populated.
struct ObservedModificationCalls {
    /// (pos_in_read, confidence ∈ [0,1]) pairs for positions where
    /// a tool reported a modification. Empty ⇒ no mod data available.
    std::vector<std::pair<std::uint32_t, float>> m6a_calls;
    std::vector<std::pair<std::uint32_t, float>> a_to_i_calls;
    std::vector<std::pair<std::uint32_t, float>> c_to_u_calls;
};

// ===========================================================================
// ComputeFactors — main entry point.
//
// Self-contained; reads only what the caller passes. The factor
// implementations (one per L_x) live in the .cpp; each is independently
// testable.
// ===========================================================================

LikelihoodFactors ComputeFactors(
    const ReadContext& read,
    const anchor::AnchorRecord& anchor,
    const ObservedModificationCalls& mods,
    const TissueContext& tissue,
    const FactorDisableMask& disabled = {});

/// Overload that consults the ExpressionDb for L_expression_prior /
/// L_depth_coverage / L_barcode_context. When `expr_db` is null, falls
/// through to the no-expression version above.
LikelihoodFactors ComputeFactorsWithExpression(
    const ReadContext& read,
    const anchor::AnchorRecord& anchor,
    const ObservedModificationCalls& mods,
    const TissueContext& tissue,
    const ExpressionDb* expr_db,
    const FactorDisableMask& disabled = {});

}  // namespace llmap::fusion
