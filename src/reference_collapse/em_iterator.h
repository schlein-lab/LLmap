// LLmap — EM Iterator for WaveCollapse updates.
//
// Implements one iteration step of the WaveCollapse EM algorithm:
//   P_{t+1}(b|r) = (1-γ) P_t(b|r) + γ · Z⁻¹ [
//       L(r|b) · λ_t(b) · π_AI(b|r) · π_bio(b) · Σ_{b'∈N(b)} K(b,b') · P_t(b'|r)
//   ]
//
// Components:
//   - P-update: main probability update with all likelihood terms
//   - λ-update: coverage prior λ_t(b) = Σ_r P_t(b|r)
//   - K-smoothing: spatial coupling via Gaussian kernel
//
// CPU fallback implementation; CUDA version in em_iterator.cu (future).

#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector>

#include "core/wave_state.h"

namespace llmap {

// Smoothing kernel configuration for K(b, b')
struct SmoothingKernelConfig {
    float sigma_genome_bp = 100'000.0f;  // Gaussian σ in base pairs
    float min_weight = 0.001f;           // Truncate kernel below this
    std::uint32_t max_neighbors = 20;    // Max neighbors per bucket
};

// Configuration for one EM iteration
struct EmIteratorConfig {
    float gamma = 0.3f;                   // Damping (γ), higher = faster but less stable
    float tau_collapse = 0.99f;           // Collapse threshold
    bool apply_smoothing = true;          // Enable K(b,b') spatial coupling
    SmoothingKernelConfig kernel_config;

    // Weight multipliers for likelihood terms (for ablation/tuning)
    float weight_seq_likelihood = 1.0f;   // L(r|b)
    float weight_coverage = 1.0f;         // λ_t(b)
    float weight_ai_prior = 1.0f;         // π_AI(b|r)
    float weight_bio_prior = 1.0f;        // π_bio(b)
};

// Statistics from one EM iteration
struct EmIterationStats {
    std::uint32_t reads_processed = 0;
    std::uint32_t reads_collapsed = 0;
    std::uint32_t reads_active = 0;
    float max_prob_delta = 0.0f;          // Largest probability change
    float mean_entropy = 0.0f;            // Mean entropy across reads
    float iteration_time_ms = 0.0f;
};

// Pre-computed neighbor structure for K-smoothing
// Stores which buckets are neighbors of each bucket
struct BucketNeighborhood {
    std::vector<std::uint32_t> offsets;   // [num_buckets + 1], CSR style
    std::vector<std::uint32_t> neighbors; // neighbor bucket IDs
    std::vector<float> weights;           // K(b, b') weights

    [[nodiscard]] std::span<const std::uint32_t>
    GetNeighbors(std::uint32_t bucket_id) const;

    [[nodiscard]] std::span<const float>
    GetWeights(std::uint32_t bucket_id) const;

    [[nodiscard]] std::size_t NumBuckets() const noexcept {
        return offsets.empty() ? 0 : offsets.size() - 1;
    }
};

// Build neighborhood structure from bucket positions
// positions: [num_buckets] genomic position for each bucket
BucketNeighborhood BuildNeighborhood(
    std::span<const std::uint64_t> positions,
    const SmoothingKernelConfig& config);

// Likelihood function type: L(r|b) sequence likelihood
// Takes (read_idx, bucket_id) and returns likelihood score [0, 1]
using LikelihoodFn = std::function<float(std::uint32_t, std::uint32_t)>;

// AI prior function type: π_AI(b|r) from foundation model
// Takes (read_idx, bucket_id) and returns prior [0, 1]
using AiPriorFn = std::function<float(std::uint32_t, std::uint32_t)>;

// Biology prior array: π_bio(b) indexed by bucket_id
// If empty, all buckets have uniform prior 1.0

// Main EM iterator class
class EmIterator {
public:
    explicit EmIterator(const EmIteratorConfig& config = {});
    ~EmIterator();

    EmIterator(const EmIterator&) = delete;
    EmIterator& operator=(const EmIterator&) = delete;
    EmIterator(EmIterator&&) noexcept;
    EmIterator& operator=(EmIterator&&) noexcept;

    // --- Setup before iterations ---

    // Set likelihood function L(r|b)
    void SetLikelihoodFn(LikelihoodFn fn);

    // Set AI prior function π_AI(b|r)
    void SetAiPriorFn(AiPriorFn fn);

    // Set biology prior array π_bio(b), indexed by bucket_id
    void SetBiologyPrior(std::span<const float> bio_prior);

    // Set neighborhood structure for K-smoothing
    void SetNeighborhood(const BucketNeighborhood& neighborhood);

    // --- Per-iteration operations ---

    // Compute coverage prior λ_t(b) = Σ_r P_t(b|r)
    // Returns vector indexed by bucket_id
    std::vector<float> ComputeCoveragePrior(
        const WaveState& state,
        std::uint32_t num_buckets) const;

    // Apply K-smoothing to probabilities
    // Updates probabilities in-place, weighted by neighbor coupling
    void ApplyKSmoothing(
        WaveState& state,
        const std::vector<float>& coverage_prior) const;

    // Perform one full EM iteration step
    // Updates state in-place, returns stats
    EmIterationStats Step(
        WaveState& state,
        std::uint32_t num_buckets);

    // Perform collapse check and mark converged reads
    // Returns number of newly collapsed reads
    std::uint32_t CheckAndCollapse(WaveState& state);

    // --- Accessors ---

    const EmIteratorConfig& Config() const { return config_; }
    void SetConfig(const EmIteratorConfig& config) { config_ = config; }

private:
    EmIteratorConfig config_;

    LikelihoodFn likelihood_fn_;
    AiPriorFn ai_prior_fn_;
    std::vector<float> bio_prior_;
    const BucketNeighborhood* neighborhood_{nullptr};

    // Compute unnormalized update for one read
    // Returns vector of (bucket_idx_in_read, new_prob) pairs
    void ComputeReadUpdate(
        const WaveState& state,
        std::uint32_t read_idx,
        const std::vector<float>& coverage_prior,
        std::vector<float>& out_probs) const;

    // Compute smoothed probability contribution from neighbors
    float ComputeNeighborContribution(
        const WaveState& state,
        std::uint32_t read_idx,
        std::uint32_t bucket_id) const;
};

// Compute Shannon entropy of a probability distribution
float ComputeEntropy(std::span<const float> probs);

// Utility: normalize probability array in-place
void NormalizeProbabilities(std::span<float> probs);

}  // namespace llmap
