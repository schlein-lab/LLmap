// LLmap — Member Propagation: propagate positions from representatives to members.
//
// After Stage 2 EM refinement converges for cluster representatives, this module
// propagates the refined positions back to all cluster members. Each member
// receives a position distribution based on:
//   1. The representative's converged distribution
//   2. The member's similarity to the representative
//   3. The cluster's overall cohesion
//
// Members with high similarity to their representative get more confident
// positions; outliers retain higher uncertainty.

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "core/wave_state.h"

namespace llmap::self_interference {
struct ClusteringResult;
struct ClusterRepResult;
struct SelfWaveCollapseResult;
class SimilarityGraph;
}  // namespace llmap::self_interference

namespace llmap {

// Configuration for member propagation
struct MemberPropagationConfig {
    float base_confidence_scaling = 0.8f;   // Member's confidence = rep_confidence * scaling
    float similarity_weight = 0.5f;         // How much similarity affects confidence
    float min_confidence = 0.1f;            // Minimum confidence for any member
    float max_confidence = 0.95f;           // Maximum confidence for members (< 1.0)

    float cohesion_weight = 0.3f;           // How much cluster cohesion affects propagation
    float outlier_penalty = 0.5f;           // Confidence reduction for low-similarity members

    bool propagate_top_candidates_only = true;   // Only propagate representative's top buckets
    std::uint32_t max_propagated_candidates = 10; // Max buckets to propagate per member
};

// Statistics from member propagation
struct MemberPropagationStats {
    std::uint32_t num_clusters = 0;
    std::uint32_t num_representatives = 0;
    std::uint32_t num_members_propagated = 0;
    std::uint32_t num_members_skipped = 0;    // Already collapsed or no rep

    float avg_member_confidence = 0.0f;
    float avg_similarity_to_rep = 0.0f;
    float avg_cluster_cohesion = 0.0f;

    std::uint32_t total_buckets_propagated = 0;
    float propagation_time_ms = 0.0f;
};

// Per-member propagation result
struct MemberPropagationEntry {
    std::uint32_t member_idx;               // Read index of member
    std::uint32_t cluster_id;               // Cluster this member belongs to
    std::uint32_t representative_idx;       // Read index of representative
    float similarity_to_rep;                // Member's similarity to representative
    float propagated_confidence;            // Confidence of propagated position
    std::uint32_t buckets_received;         // Number of buckets propagated
};

// Result of member propagation
struct MemberPropagationResult {
    std::vector<MemberPropagationEntry> entries;
    MemberPropagationStats stats;

    // Get entry for a specific member (or nullopt if not found)
    [[nodiscard]] std::optional<MemberPropagationEntry> GetEntry(std::uint32_t member_idx) const;

    // Get all members in a cluster
    [[nodiscard]] std::vector<std::uint32_t> GetClusterMembers(std::uint32_t cluster_id) const;
};

// Member propagation engine
//
// Stage 2 step 2.x: After representatives have converged positions from
// WaveCollapse EM, propagate those positions to all cluster members.
// This allows cluster members to inherit refined genomic positions without
// running full EM on each individual read.
class MemberPropagation {
public:
    MemberPropagation();
    explicit MemberPropagation(const MemberPropagationConfig& config);
    ~MemberPropagation();

    MemberPropagation(MemberPropagation&&) noexcept;
    MemberPropagation& operator=(MemberPropagation&&) noexcept;
    MemberPropagation(const MemberPropagation&) = delete;
    MemberPropagation& operator=(const MemberPropagation&) = delete;

    // Propagate positions from representatives to members
    //
    // Parameters:
    //   state: WaveState containing representative positions (will be modified)
    //   clustering: Cluster assignments for all reads
    //   rep_result: Representative selection result
    //   graph: Similarity graph for computing member-rep similarity
    //
    // Returns: statistics and per-member propagation info
    MemberPropagationResult Propagate(
        WaveState& state,
        const self_interference::ClusteringResult& clustering,
        const self_interference::ClusterRepResult& rep_result,
        const self_interference::SimilarityGraph& graph);

    // Propagate for a single cluster
    //
    // Parameters:
    //   state: WaveState to modify
    //   representative_idx: Read index of representative
    //   member_indices: Read indices of all cluster members (including representative)
    //   graph: Similarity graph
    //   cluster_id: Cluster ID for tracking
    //
    // Returns: entries for propagated members (excludes representative)
    std::vector<MemberPropagationEntry> PropagateCluster(
        WaveState& state,
        std::uint32_t representative_idx,
        std::span<const std::uint32_t> member_indices,
        const self_interference::SimilarityGraph& graph,
        std::uint32_t cluster_id);

    // Propagate to a single member
    //
    // Returns: number of buckets propagated (0 if skipped)
    std::uint32_t PropagateMember(
        WaveState& state,
        std::uint32_t representative_idx,
        std::uint32_t member_idx,
        float similarity_to_rep,
        float cluster_cohesion);

    // Compute cluster cohesion (avg pairwise similarity among members)
    [[nodiscard]] float ComputeClusterCohesion(
        const self_interference::SimilarityGraph& graph,
        std::span<const std::uint32_t> member_indices) const;

    // Get/set configuration
    [[nodiscard]] const MemberPropagationConfig& Config() const { return config_; }
    void SetConfig(const MemberPropagationConfig& config) { config_ = config; }

private:
    MemberPropagationConfig config_;

    struct InternalState;
    std::unique_ptr<InternalState> state_;

    // Compute confidence for a member based on similarity and cohesion
    [[nodiscard]] float ComputeMemberConfidence(
        float rep_confidence,
        float similarity_to_rep,
        float cluster_cohesion) const;

    // Scale probabilities by confidence factor
    void ScaleProbabilities(
        std::vector<BucketProb>& buckets,
        float confidence_factor) const;
};

// Convenience function: propagate with default config
MemberPropagationResult PropagateToMembers(
    WaveState& state,
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph);

}  // namespace llmap
