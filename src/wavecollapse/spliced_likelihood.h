// LLmap — Spliced-aware likelihood + neighborhood kernel for
// Transcript-Mode Wave-Collapse.
//
// The classical EM iterator (src/reference_collapse/em_iterator.cpp)
// implements
//
//     P_{t+1}(b | r) = (1-γ) P_t + γ · Z⁻¹ [L(r|b) · λ_t · π_AI · π_bio · K]
//
// where L is sequence likelihood and K is a Gaussian on genome distance.
// For Transcript-Mode we extend BOTH terms:
//
//   * L_spliced(r | b) — sequence likelihood factorised over exon
//     segments, with per-junction probability multiplied in.
//   * K_spliced(b, b') — Gaussian on positional distance × splice-site
//     canonicality × DB-evidence × tissue-compatibility.
//
// This module is COMPLETELY STANDALONE — it does not modify
// em_iterator.cpp, so DNA-mode WaveCollapse continues to behave
// bit-for-bit identically. The Transcript-Mode dispatcher (lands in
// Block 6) calls these new functions when `mode == Transcript`.
//
// Everything here is pure-functional + side-effect-free. Each helper
// takes plain inputs (no mutation, no globals) so tests can pin the
// exact numerical behaviour.

#pragma once

#include "anchor/anchor_record.h"
#include "annot/junction_db.h"
#include "fusion/likelihood_factors.h"

#include <cstdint>
#include <string_view>

namespace llmap::wavecollapse {

// ===========================================================================
// Per-segment likelihood for a spliced read.
//
// Inputs:
//   per_segment_likelihoods — L(r_seg | anchor_seg) for each exon segment.
//   junction_probs          — L_junc(j) ∈ [0,1] for each junction between
//                              segments. junctions.size() == segments - 1.
//
// Output: ∏ per_segment_likelihoods × ∏ junction_probs, floor-clamped.
//
// This is the multiplicative form prescribed by Plan-Block 5.
// ===========================================================================

[[nodiscard]] float SplicedLikelihoodFactored(
    const std::vector<float>& per_segment_likelihoods,
    const std::vector<float>& junction_probs);

// ===========================================================================
// Junction probability — combines splice-site PWM evidence (already
// embedded in the anchor's ExonBoundary::donor_score / acceptor_score)
// with multi-source DB evidence and optional tissue prior.
//
// Tissue prior is a simple per-tissue scaling factor; pass 1.0 when no
// tissue context is available.
// ===========================================================================

[[nodiscard]] float JunctionProbability(
    const anchor::ExonBoundary& boundary,
    const annot::JunctionEvidence& db_evidence,
    float tissue_compatibility = 1.0f);

// ===========================================================================
// K_spliced(b, b') — extended neighborhood kernel.
//
// Inputs:
//   delta_pos_bp     — genomic distance between bucket centres.
//   gauss_scale_bp   — Gaussian σ for positional decay
//                        (matches existing em_iterator config; pass
//                         the same value the DNA-mode kernel uses).
//   junction_prob    — probability that the junction between the two
//                        buckets is real (canonical + annotated etc.).
//                        Use 1.0 if the two buckets are inside the
//                        same exon (no junction between them).
//   tissue_compat    — extra multiplicative factor for tissue context.
//
// Output:
//   K_spliced ∈ [0, 1].
// ===========================================================================

[[nodiscard]] float KSpliced(
    std::uint64_t delta_pos_bp,
    float gauss_scale_bp,
    float junction_prob,
    float tissue_compat = 1.0f);

// ===========================================================================
// Tissue-compatibility helper — turns a TPM into a [0.2, 1.0] factor
// suitable for use as `tissue_compatibility` in the kernel. Same shape
// as L_expression_prior but more compressed (we don't want the kernel
// to penalise as hard as the sequence likelihood does).
// ===========================================================================

[[nodiscard]] float TissueCompatibility(float tpm);

}  // namespace llmap::wavecollapse
