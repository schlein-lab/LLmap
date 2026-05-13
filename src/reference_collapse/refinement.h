// LLmap — Refinement: coarse→fine bucket expansion during WaveCollapse.
//
// After EM iterations converge at a coarse level (L0 or L1), high-probability
// buckets are expanded to their children at the finer level (L1 or L2).
// Probabilities are re-initialized proportional to the parent probability,
// distributed uniformly or weighted by AI prior among children.
//
// Workflow:
//   1. Check refinement readiness (enough reads above threshold)
//   2. Build child index from BucketPyramid (once per level transition)
//   3. For each read above threshold: expand candidates to child buckets
//   4. Update WaveState level tracking
//   5. Continue EM at the new resolution

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "core/bucket_pyramid.h"
#include "core/wave_state.h"

namespace llmap {

// Configuration for refinement step
struct RefinementConfig {
    float expansion_threshold = 0.01f;   // Expand buckets with P >= this
    std::uint32_t max_candidates = 50;   // Max candidates after expansion
    bool weight_by_ai_prior = false;     // Use AI prior to weight children
    bool preserve_relative_probs = true; // Children get P proportional to parent
};

// Statistics from one refinement step
struct RefinementStats {
    std::uint32_t reads_refined = 0;
    std::uint32_t reads_skipped = 0;      // Already collapsed or wrong level
    std::uint32_t buckets_expanded = 0;
    std::uint32_t total_children_added = 0;
    float avg_expansion_factor = 0.0f;    // Avg children per expanded bucket
    float refinement_time_ms = 0.0f;
};

// Pre-computed parent→children mapping for fast expansion
// Built from BucketPyramid for a specific level transition
class ChildIndex {
public:
    ChildIndex() = default;

    // Build L0→L1 child index
    static ChildIndex BuildL0ToL1(const BucketPyramid& pyramid);

    // Build L1→L2 child index
    static ChildIndex BuildL1ToL2(const BucketPyramid& pyramid);

    // Get children of a parent bucket
    [[nodiscard]] std::span<const std::uint32_t> GetChildren(std::uint32_t parent_id) const;

    // Number of parent buckets
    [[nodiscard]] std::size_t NumParents() const noexcept {
        return offsets_.empty() ? 0 : offsets_.size() - 1;
    }

    // Total children across all parents
    [[nodiscard]] std::size_t TotalChildren() const noexcept {
        return children_.size();
    }

    // Check if index is valid
    [[nodiscard]] bool IsValid() const noexcept { return !offsets_.empty(); }

private:
    std::vector<std::uint32_t> offsets_;   // CSR offsets [num_parents + 1]
    std::vector<std::uint32_t> children_;  // Child bucket IDs
};

// Main refinement class
class Refinement {
public:
    explicit Refinement(const RefinementConfig& config = {});
    ~Refinement() = default;

    Refinement(const Refinement&) = default;
    Refinement& operator=(const Refinement&) = default;
    Refinement(Refinement&&) = default;
    Refinement& operator=(Refinement&&) = default;

    // Check if state is ready for refinement
    // Returns number of reads that would be refined
    [[nodiscard]] std::uint32_t CountRefinementCandidates(
        const WaveState& state,
        WaveLevel from_level) const;

    // Refine reads from L0 to L1
    // Modifies state in-place, returns stats
    RefinementStats RefineL0ToL1(
        WaveState& state,
        const ChildIndex& child_index);

    // Refine reads from L1 to L2
    RefinementStats RefineL1ToL2(
        WaveState& state,
        const ChildIndex& child_index);

    // Refine reads from specified level to next finer level
    // Uses pre-built child index
    RefinementStats Refine(
        WaveState& state,
        WaveLevel from_level,
        const ChildIndex& child_index);

    // Refine a single read
    // Returns number of children added (0 if read not refined)
    std::uint32_t RefineRead(
        WaveState& state,
        std::uint32_t read_idx,
        WaveLevel from_level,
        const ChildIndex& child_index);

    // Get/set configuration
    const RefinementConfig& Config() const { return config_; }
    void SetConfig(const RefinementConfig& config) { config_ = config; }

    // Set expansion threshold
    void SetExpansionThreshold(float threshold) {
        config_.expansion_threshold = threshold;
    }

    // Set AI prior weights for child distribution (optional)
    // Indexed by child bucket ID, used when weight_by_ai_prior=true
    void SetAiPrior(std::span<const float> ai_prior);

private:
    RefinementConfig config_;
    std::vector<float> ai_prior_;

    // Expand one parent bucket to its children
    // Returns vector of (child_id, probability) pairs
    std::vector<BucketProb> ExpandBucket(
        std::uint32_t parent_id,
        float parent_prob,
        const ChildIndex& child_index) const;
};

// Convenience: check if state should transition to next level
// Based on fraction of reads collapsed or above high probability
bool ShouldRefine(
    const WaveState& state,
    WaveLevel current_level,
    float collapsed_fraction_threshold = 0.8f,
    float high_prob_threshold = 0.9f);

// Get next finer level (returns L3 if already at L2)
WaveLevel NextLevel(WaveLevel level);

// Get coarser level (returns L0 if already at L0)
WaveLevel PrevLevel(WaveLevel level);

}  // namespace llmap
