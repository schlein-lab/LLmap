// LLmap — Spliced-aware likelihood + kernel implementation.
//
// All scoring functions follow the Plan-Block 5 contract: per-segment
// products, junction-probability multiplication, lossless floor-clamp
// from fusion::kFloor.

#include "wavecollapse/spliced_likelihood.h"

#include <algorithm>
#include <cmath>

namespace llmap::wavecollapse {

namespace {

// Floor matches the Multi-Signal Fusion engine's per-factor floor.
constexpr float kFloor = fusion::kFloor;

inline float Sigmoid(float x) noexcept {
    return 1.0f / (1.0f + std::exp(-x));
}

}  // namespace

// ---------------------------------------------------------------------------
// SplicedLikelihoodFactored
// ---------------------------------------------------------------------------

float SplicedLikelihoodFactored(
    const std::vector<float>& per_segment_likelihoods,
    const std::vector<float>& junction_probs) {

    // No segments at all → unmappable; fall through with floor so EM
    // can still move the bucket to a better state without an
    // information-theoretic singularity.
    if (per_segment_likelihoods.empty()) return kFloor;

    // Segments without junctions = single segment, no junctions allowed.
    if (junction_probs.size() + 1 != per_segment_likelihoods.size()
        && !junction_probs.empty()) {
        // Shape mismatch — treat as low-confidence.
        return kFloor;
    }

    float product = 1.0f;
    for (float ls : per_segment_likelihoods) {
        product *= std::clamp(ls, kFloor, 1.0f);
    }
    for (float lj : junction_probs) {
        product *= std::clamp(lj, kFloor, 1.0f);
    }
    return std::max(kFloor, std::min(1.0f, product));
}

// ---------------------------------------------------------------------------
// JunctionProbability
// ---------------------------------------------------------------------------

float JunctionProbability(const anchor::ExonBoundary& b,
                            const annot::JunctionEvidence& ev,
                            float tissue_compat) {

    // Splice-site PWM scores already in [0,1]; missing values (0) get
    // a neutral mid-point so brand-new anchors aren't penalised.
    const float ds = (b.donor_score    > 0.0f) ? b.donor_score    : 0.5f;
    const float as = (b.acceptor_score > 0.0f) ? b.acceptor_score : 0.5f;
    float p = ds * as;

    // DB-evidence multiplier — Plan-Block 4.5 sigmoidal table.
    // Annotated junctions get a modest boost; un-supported keep PWM
    // baseline. Floor kept above zero per lossless contract.
    float db_mult = 0.60f;  // default for un-supported junctions
    if (ev.in_gencode) db_mult = 1.05f;
    else if (ev.in_chessdb) db_mult = 0.95f;
    else if (ev.in_gtex)    db_mult = 0.90f;
    else if (ev.in_circ_db) db_mult = 0.95f;
    p *= db_mult;

    // Tissue-compatibility scaling — never below 0.20.
    p *= std::max(0.20f, std::min(1.0f, tissue_compat));

    // Back-splice junctions (spliceosome_class == 3) are valid only for
    // circRNA loci; we don't penalise them here (the AlignmentRecord
    // status differentiates), but we clamp aggressively to avoid runaway
    // products.
    if (b.spliceosome_class == 3) {
        p = std::min(p, 0.95f);
    }

    return std::max(kFloor, std::min(1.0f, p));
}

// ---------------------------------------------------------------------------
// KSpliced — Gaussian × junction × tissue
// ---------------------------------------------------------------------------

float KSpliced(std::uint64_t delta_pos_bp,
                float gauss_scale_bp,
                float junction_prob,
                float tissue_compat) {

    if (gauss_scale_bp <= 0.0f) return kFloor;

    // Gaussian on genomic distance.
    const float x = static_cast<float>(delta_pos_bp) / gauss_scale_bp;
    const float gauss = std::exp(-0.5f * x * x);

    // Composite kernel: classical Gaussian × junction × tissue, with
    // all factors clamped to [floor, 1] to prevent any zero from
    // collapsing the product.
    const float jp = std::clamp(junction_prob, kFloor, 1.0f);
    const float tc = std::clamp(tissue_compat, kFloor, 1.0f);
    return std::max(kFloor, std::min(1.0f, gauss * jp * tc));
}

// ---------------------------------------------------------------------------
// TissueCompatibility
// ---------------------------------------------------------------------------

float TissueCompatibility(float tpm) {
    if (tpm < 0.0f) return 1.0f;  // unknown TPM → neutral

    // Compressed sigmoid in [0.20, 1.0] — matches Plan-Block 4.5
    // requirement that the kernel never penalises as hard as the
    // sequence likelihood.
    const float v = Sigmoid((std::log10(tpm + 1.0f) - 1.0f) / 0.7f);
    return std::max(0.20f, std::min(1.0f, v));
}

}  // namespace llmap::wavecollapse
