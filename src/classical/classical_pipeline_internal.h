// LLmap — ClassicalPipeline internal header for split implementation.
#pragma once

#include "classical/classical_pipeline.h"

#include <chrono>
#include <cmath>

namespace llmap::classical {

// Helper to measure time (shared across implementation files)
class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    float ElapsedUs() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<float, std::micro>(now - start_).count();
    }

    float ElapsedMs() const {
        return ElapsedUs() / 1000.0f;
    }

private:
    std::chrono::steady_clock::time_point start_;
};

// Compute MAPQ from alignment scores and confidence metrics.
// Uses a formula similar to minimap2: MAPQ reflects the probability
// that the alignment is wrong, based on score gap to secondary.
//
// Parameters:
//   primary_score: Alignment score of the primary (best) alignment
//   secondary_score: Alignment score of the best secondary (0 if none)
//   identity: Alignment identity [0.0, 1.0]
//   num_secondaries: Number of secondary alignments found
//   query_len: Length of the query sequence
//
// Returns: MAPQ in range [0, 60]
inline uint32_t ComputeMapq(
    int32_t primary_score,
    int32_t secondary_score,
    float identity,
    uint32_t num_secondaries,
    uint32_t query_len) {

    constexpr uint32_t kMaxMapq = 60;

    // No alignment or invalid score
    if (primary_score <= 0) {
        return 0;
    }

    // Low identity = low confidence
    if (identity < 0.5f) {
        return 0;
    }

    // Unique alignment (no secondaries): base MAPQ on score and identity
    if (num_secondaries == 0 || secondary_score <= 0) {
        // Score-based component: higher score = higher confidence
        // Normalize by query length for fair comparison
        float score_density = static_cast<float>(primary_score) /
                             static_cast<float>(std::max(query_len, 1u));

        // Identity-based component
        float identity_factor = identity;

        // Combined confidence - scale to MAPQ
        // High-scoring unique alignments get MAPQ=60
        float confidence = score_density * identity_factor;

        if (confidence >= 1.0f && identity >= 0.9f) {
            return kMaxMapq;  // Maximum confidence
        } else if (confidence >= 0.7f && identity >= 0.8f) {
            return 50;
        } else if (confidence >= 0.5f && identity >= 0.7f) {
            return 40;
        } else {
            return std::min(kMaxMapq,
                std::max(10u, static_cast<uint32_t>(confidence * 40.0f)));
        }
    }

    // Multi-mapping: compute MAPQ from score gap
    // Formula: MAPQ = -10 * log10(P_error)
    // P_error estimated from score difference using Phred-like scaling
    int32_t score_diff = primary_score - secondary_score;

    if (score_diff <= 0) {
        // Primary and secondary have same score = ambiguous mapping
        return 0;
    }

    // Estimate probability that secondary is correct:
    // P(secondary correct) ~ exp(-score_diff / scale)
    // For gap-affine scoring, a typical scale is ~20-50 points
    constexpr float kScoreScale = 30.0f;
    float p_secondary = std::exp(-static_cast<float>(score_diff) / kScoreScale);

    // Cap P_secondary considering multiple hits dilute confidence
    if (num_secondaries > 1) {
        // More secondaries = more uncertainty, but diminishing returns
        float multi_factor = 1.0f + 0.3f * std::log2(static_cast<float>(num_secondaries));
        p_secondary = std::min(0.5f, p_secondary * multi_factor);
    }

    // Convert to MAPQ: MAPQ = -10 * log10(P_error)
    // Very small p_secondary means high confidence
    if (p_secondary < 1e-6f) {
        return kMaxMapq;  // Effectively unique
    }

    // Clamp p_secondary to avoid log10(0)
    p_secondary = std::max(p_secondary, 1e-6f);

    float mapq_raw = -10.0f * std::log10(p_secondary);

    // Apply identity penalty for lower-confidence alignments
    if (identity < 0.9f) {
        mapq_raw *= identity;
    }

    // Clamp to valid range [0, 60]
    if (mapq_raw < 0.0f) return 0;
    if (mapq_raw > static_cast<float>(kMaxMapq)) return kMaxMapq;
    return static_cast<uint32_t>(mapq_raw);
}

}  // namespace llmap::classical
