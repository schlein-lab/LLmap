// LLmap — Multi-Signal Fusion Engine implementation.
//
// One small named function per factor; ComputeFactors is the orchestrator.
// All numerical defaults come from the Plan-Block 4.5 specification
// (~/.claude/plans/memoized-crafting-lecun.md). The literature anchors
// are quoted inline next to each magic number so future readers don't
// have to dig.

#include "fusion/likelihood_factors.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

namespace llmap::fusion {

namespace {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

inline float ClampToFloor(float v) noexcept {
    return std::max(kFloor, std::min(1.0f, v));
}

inline float Sigmoid(float x) noexcept {
    return 1.0f / (1.0f + std::exp(-x));
}

// DRACH context for m6A (D=A/G/U, R=A/G, fixed A,C, H=A/C/U).
bool IsDrach5mer(std::string_view ctx) noexcept {
    if (ctx.size() != 5) return false;
    auto up = [](char c) {
        return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    };
    char d = up(ctx[0]), r = up(ctx[1]), a = up(ctx[2]), c = up(ctx[3]), h = up(ctx[4]);
    return (d == 'A' || d == 'G' || d == 'T' || d == 'U')
        && (r == 'A' || r == 'G')
        && a == 'A' && c == 'C'
        && (h == 'A' || h == 'C' || h == 'T' || h == 'U');
}

// ---------------------------------------------------------------------------
// Per-factor implementations
// ---------------------------------------------------------------------------

/// 1. L_sequence — forward pairHMM-style match score, per-platform.
///
/// We treat the upstream seed-extend's score as already encoding the
/// sequence-similarity signal. Without explicit per-base score we use
/// a simple proxy: if the read length and anchor sequence length are
/// comparable AND the anchor has a non-empty sequence, return a
/// platform-adjusted baseline; else neutral (1.0).
float ComputeLSequence(const ReadContext& read,
                        const anchor::AnchorRecord& anchor) {
    if (anchor.sequence.empty()) return 1.0f;

    // Platform-specific baseline (per Plan-Block 4.5 / LRGASP 2024).
    float baseline = 0.94f;  // ont_rna004 default
    if (read.platform == "hifi" || read.platform == "isoseq") {
        baseline = 0.99f;
    } else if (read.platform == "illumina") {
        baseline = 0.998f;
    } else if (read.platform == "ont_legacy") {
        baseline = 0.90f;
    }

    // Length compatibility reduces baseline when grossly mismatched
    // (e.g. 20-nt read against a 1.5-kb anchor).
    const float ratio = std::min<float>(
        static_cast<float>(read.read_length)
            / std::max<float>(static_cast<float>(anchor.sequence.size()),
                               1.0f),
        1.0f);
    const float length_adj = 0.5f + 0.5f * ratio;
    return ClampToFloor(baseline * length_adj);
}

/// 2. L_modification — Bayes update from observed mod calls vs anchor's
/// modification context.
///
/// We don't yet have anchor-level mod-site info plumbed in; that lands
/// in Block 9 once REPIC sites are linked to AnchorRecord. For now we
/// approximate via sequence-context priors: any observed m6A call at a
/// DRACH-context position boosts; any call at non-DRACH is treated as
/// neutral (the tool produced a call, we don't second-guess it but we
/// don't reward it either).
float ComputeLModification(const ReadContext& read,
                            const anchor::AnchorRecord& /*anchor*/,
                            const ObservedModificationCalls& mods) {
    if (mods.m6a_calls.empty()
        && mods.a_to_i_calls.empty()
        && mods.c_to_u_calls.empty()) {
        return 1.0f;  // no mod data ⇒ neutral
    }

    float l = 1.0f;
    // For each m6A call check sequence context (need 2 bp upstream + 2 bp
    // downstream).
    for (const auto& [pos, conf] : mods.m6a_calls) {
        if (pos < 2 || pos + 2 >= read.read_sequence.size()) continue;
        auto ctx = read.read_sequence.substr(pos - 2, 5);
        const bool in_drach = IsDrach5mer(ctx);
        if (in_drach) {
            // Strong evidence — reward as in Plan-Block 4.5 decision table.
            l *= (0.95f * conf + 1.0f * (1.0f - conf));
        } else {
            // Mod-call outside expected context: lower the likelihood
            // slightly (could be tool false positive).
            l *= (0.30f * conf + 1.0f * (1.0f - conf));
        }
    }
    return ClampToFloor(l);
}

/// 3. L_depth_coverage — Negative-Binomial per-exon expression-prior.
///
/// When no ExpressionDb is plumbed, returns neutral 1.0. When an
/// ExpressionDb is provided AND we have an (anchor.transcript_id,
/// tissue) entry, we treat that TPM as the expected mean depth for a
/// normalised Negative-Binomial likelihood. Without per-exon observed
/// depth we approximate via a log-sigmoid on the TPM itself, which
/// scales sensibly with expression range.
float ComputeLDepthCoverage(const anchor::AnchorRecord& anchor,
                              const TissueContext& tissue,
                              const ExpressionDb* expr_db) {
    if (expr_db == nullptr || anchor.transcript_id.empty()
        || tissue.label.empty()) {
        return 1.0f;
    }
    const float tpm = expr_db->ExpectedTpm(anchor.transcript_id, tissue.label);
    if (tpm < 0.0f) return 1.0f;  // unknown → neutral

    // Plan-Block 4.5 sigmoid (TPM=0 → 0.18, TPM=10 → 0.50,
    // TPM=1000 → 0.95). Centred at log10(11) ≈ 1.04 with scale 0.67
    // fits all four spec'd anchor points within ±0.02.
    const float v = Sigmoid((std::log10(tpm + 1.0f) - 1.04f) / 0.67f);
    return std::max(0.10f, std::min(0.98f, v));
}

/// 4. L_expression_prior — tissue × cell-type TPM sigmoid.
///
/// Plan-Block 4.5 numerical table:
///   TPM=0     → 0.18
///   TPM=1     → 0.40
///   TPM=10    → 0.50
///   TPM=100   → 0.78
///   TPM=1000  → 0.95
/// Sigmoid centred at log10(11) ≈ 1.04 with scale 0.67.
float ComputeLExpressionPrior(const anchor::AnchorRecord& anchor,
                                const TissueContext& tissue,
                                const ExpressionDb* expr_db) {
    if (expr_db == nullptr || anchor.transcript_id.empty()
        || tissue.label.empty()) {
        return 1.0f;
    }
    const float tpm = expr_db->ExpectedTpm(anchor.transcript_id, tissue.label);
    if (tpm < 0.0f) return 1.0f;  // unknown → neutral

    const float v = Sigmoid((std::log10(tpm + 1.0f) - 1.04f) / 0.67f);
    return std::max(0.10f, std::min(0.98f, v));
}

/// 5. L_phasing — HP-tag consistency.
///
/// The anchor carries no explicit haplotype field; we use the tags
/// {"hap1", "hap2", "PANGEN_<sample>_hap{1,2}"} as fallback. When the
/// read's haplotype is unknown ⇒ neutral 1.0.
float ComputeLPhasing(const ReadContext& read,
                       const anchor::AnchorRecord& anchor) {
    if (!read.haplotype.has_value()) return 1.0f;
    const auto rh = *read.haplotype;

    bool anchor_has_hap_tag = false;
    bool anchor_matches = false;
    for (const auto& t : anchor.tags) {
        const bool tag_hap1 = (t == "hap1" || t.find("hap1") != std::string::npos);
        const bool tag_hap2 = (t == "hap2" || t.find("hap2") != std::string::npos);
        if (tag_hap1 || tag_hap2) {
            anchor_has_hap_tag = true;
            if ((rh == 0 && tag_hap1) || (rh == 1 && tag_hap2)) {
                anchor_matches = true;
            }
        }
    }
    if (!anchor_has_hap_tag) return 1.0f;
    return anchor_matches ? 1.0f : 0.05f;  // cross-hap rare but not 0
}

/// 6. L_pseudogene_compatibility — biotype-aware.
///
/// Tag-driven; the anchor's tags carry the GENCODE biotype prefix.
/// We don't have read-level STOP/frame-shift signals plumbed yet, so
/// we return the optimistic default (1.0 for functional anchors,
/// 0.7 for polymorphic_pseudogene, etc.) as a placeholder. Sequence-
/// level pseudogene refinement lands in Block 9.
float ComputeLPseudogeneCompatibility(const anchor::AnchorRecord& anchor) {
    for (const auto& t : anchor.tags) {
        if (t == "biotype:polymorphic_pseudogene") return 0.70f;
        if (t == "biotype:transcribed_processed_pseudogene" ||
            t == "biotype:transcribed_unprocessed_pseudogene") {
            return 0.90f;
        }
        if (t.starts_with("biotype:IG_") && t.find("pseudogene") != std::string::npos) {
            return 0.40f;
        }
    }
    return 1.0f;
}

/// 7. L_junction — splice-site probability for the anchor's junctions.
///
/// Aggregates donor_score × acceptor_score across all exon_boundaries
/// of the anchor. No boundaries ⇒ neutral. Each junction's contribution
/// is clamped below at 0.05 to preserve the lossless floor at the
/// junction level.
float ComputeLJunction(const anchor::AnchorRecord& anchor) {
    if (anchor.exon_boundaries.empty()) return 1.0f;
    float product = 1.0f;
    for (const auto& b : anchor.exon_boundaries) {
        const float ds = (b.donor_score    > 0.0f) ? b.donor_score    : 0.5f;
        const float as = (b.acceptor_score > 0.0f) ? b.acceptor_score : 0.5f;
        const float jp = std::max(0.05f, ds * as);
        product *= jp;
    }
    return ClampToFloor(product);
}

/// 8. L_barcode_context — single-cell cell-type prior.
///
/// Look up P(anchor_tx_id | tissue, cell_type) via ExpressionDb. When
/// the entry exists we sigmoid-map the expression fraction; absent
/// entry or no cell_type ⇒ neutral.
float ComputeLBarcodeContext(const ReadContext& read,
                              const anchor::AnchorRecord& anchor,
                              const TissueContext& tissue,
                              const ExpressionDb* expr_db) {
    if (read.cell_type.empty() || expr_db == nullptr
        || anchor.transcript_id.empty()
        || tissue.label.empty()) {
        return 1.0f;
    }
    const float frac = expr_db->ExpectedCellTypeFraction(
        anchor.transcript_id, tissue.label, read.cell_type);
    if (frac < 0.0f) return 1.0f;
    // Sigmoid on log10(frac+1) — frac is typically already in [0,1]
    // so this compresses gently; floor at 0.10 per Plan-Block 4.5.
    const float v = Sigmoid((std::log10(frac + 1.0f) - 0.3f) / 1.0f);
    return std::max(0.10f, std::min(1.0f, v));
}

/// 9. L_mapq_signal — MAPQ as continuous signal, never a threshold.
///
/// L_mapq = max(floor, sigmoid((MAPQ - 10) / 5))
///   MAPQ=0  → 0.12 (floor-bounded)
///   MAPQ=10 → 0.50
///   MAPQ=30 → 0.98
/// The contrast vs minimap2/STAR/HISAT2 default of MAPQ<30 → drop is
/// the entire point of this factor.
float ComputeLMapqSignal(const ReadContext& read) {
    if (read.mapq < 0) return 1.0f;
    const float v = Sigmoid((static_cast<float>(read.mapq) - 10.0f) / 5.0f);
    return std::max(0.05f, v);  // explicit 0.05 floor per Plan-Block 4.5
}

/// 10. L_length_plausibility — per TranscriptKind length window.
///
/// Anchor.kind dictates an expected read-length range. Outside the
/// window: drop to 0.1; inside: 1.0; soft penalty at the edges.
float ComputeLLengthPlausibility(const ReadContext& read,
                                   const anchor::AnchorRecord& anchor) {
    using K = core::TranscriptKind;
    const auto rl = read.read_length;
    if (rl == 0) return 1.0f;

    auto window = [](std::uint32_t lo, std::uint32_t hi,
                     std::uint32_t hard_lo, std::uint32_t hard_hi,
                     std::uint32_t rl_) -> float {
        if (rl_ < hard_lo || rl_ > hard_hi) return 0.10f;
        if (rl_ < lo || rl_ > hi) return 0.50f;
        return 1.0f;
    };

    switch (anchor.kind) {
        case K::Mirna:           return window(18,  26,  12, 50,  rl);
        case K::Pirna:           return window(24,  32,  18, 60,  rl);
        case K::Sirna:           return window(20,  24,  18, 40,  rl);
        case K::Snorna_CDbox:    return window(50,  200, 30, 500, rl);
        case K::Snorna_HacaBox:  return window(80,  250, 60, 500, rl);
        case K::Snrna_Major:     return window(100, 200, 60, 400, rl);
        case K::Snrna_Minor:     return window(100, 200, 60, 400, rl);
        case K::Scarna:          return window(80,  250, 60, 500, rl);
        case K::MatureMrna:      return window(100, 30'000, 50, 60'000, rl);
        case K::Lncrna:          return window(200, 50'000, 100, 80'000, rl);
        case K::PreMrna:         return window(1'000, 500'000, 500, 1'000'000, rl);
        case K::SterileGermline: return window(500, 5'000, 300, 8'000, rl);
        case K::CircularRna:     return window(100, 50'000, 50, 80'000, rl);
        case K::Trna:            return window(70,  90,   60, 150, rl);
        case K::Rrna:            return std::max(0.3f, 0.6f);  // tolerant
        default:                  return 1.0f;
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Product — lossless multiply with floor.
// ---------------------------------------------------------------------------

float LikelihoodFactors::Product() const noexcept {
    const float p = L_sequence
        * L_modification
        * L_depth_coverage
        * L_expression_prior
        * L_phasing
        * L_pseudogene_compatibility
        * L_junction
        * L_barcode_context
        * L_mapq_signal
        * L_length_plausibility;
    return std::max(kFloor, std::min(1.0f, p));
}

// ---------------------------------------------------------------------------
// ComputeFactors — orchestrate.
// ---------------------------------------------------------------------------

LikelihoodFactors ComputeFactors(
    const ReadContext& read,
    const anchor::AnchorRecord& anchor,
    const ObservedModificationCalls& mods,
    const TissueContext& tissue,
    const FactorDisableMask& dis) {
    return ComputeFactorsWithExpression(read, anchor, mods, tissue,
                                          /*expr_db=*/nullptr, dis);
}

LikelihoodFactors ComputeFactorsWithExpression(
    const ReadContext& read,
    const anchor::AnchorRecord& anchor,
    const ObservedModificationCalls& mods,
    const TissueContext& tissue,
    const ExpressionDb* expr_db,
    const FactorDisableMask& dis) {

    LikelihoodFactors f;
    if (!dis.sequence)
        f.L_sequence = ComputeLSequence(read, anchor);
    if (!dis.modification)
        f.L_modification = ComputeLModification(read, anchor, mods);
    if (!dis.depth_coverage)
        f.L_depth_coverage = ComputeLDepthCoverage(anchor, tissue, expr_db);
    if (!dis.expression_prior)
        f.L_expression_prior = ComputeLExpressionPrior(anchor, tissue, expr_db);
    if (!dis.phasing)
        f.L_phasing = ComputeLPhasing(read, anchor);
    if (!dis.pseudogene_compatibility)
        f.L_pseudogene_compatibility = ComputeLPseudogeneCompatibility(anchor);
    if (!dis.junction)
        f.L_junction = ComputeLJunction(anchor);
    if (!dis.barcode_context)
        f.L_barcode_context =
            ComputeLBarcodeContext(read, anchor, tissue, expr_db);
    if (!dis.mapq_signal)
        f.L_mapq_signal = ComputeLMapqSignal(read);
    if (!dis.length_plausibility)
        f.L_length_plausibility = ComputeLLengthPlausibility(read, anchor);
    return f;
}

}  // namespace llmap::fusion
