#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace llmap::self_interference {

class SimilarityGraph;
struct ClusteringResult;
struct SelfWaveCollapseResult;

// Configuration for cluster representative selection
struct ClusterRepConfig {
    // Selection method
    enum class Method {
        Medoid,           // Minimize sum of distances to all cluster members
        Centroid,         // Closest to centroid (requires embeddings)
        MaxDegree,        // Highest weighted degree within cluster
        MaxConfidence,    // Highest confidence from Self-WaveCollapse
    };

    Method method = Method::Medoid;

    // Medoid computation settings
    size_t max_medoid_candidates = 0;      // 0 = consider all members as candidates
    bool use_approximate_medoid = false;   // Use sampling for large clusters

    // Approximate medoid sampling
    size_t approx_sample_size = 100;       // Sample size for approximate medoid
    uint64_t seed = 42;                    // Random seed for sampling

    // When multiple methods available, use as tie-breaker
    bool use_confidence_tiebreaker = true;

    // Minimum cluster size to compute representative
    size_t min_cluster_size = 1;

    // Performance hints
    bool parallel = true;                  // Use parallel execution for large inputs
    size_t parallel_threshold = 1000;      // Use parallelism for >= this many clusters
};

// Statistics from representative selection
struct ClusterRepStats {
    size_t num_clusters = 0;
    size_t num_representatives = 0;
    size_t clusters_skipped = 0;           // Clusters too small

    // Per-method stats
    size_t exact_medoid_computed = 0;
    size_t approx_medoid_computed = 0;
    size_t max_degree_fallback = 0;

    // Quality metrics
    float avg_rep_confidence = 0.0f;       // Average representative confidence
    float avg_cluster_size = 0.0f;
    float total_within_cluster_distance = 0.0f;

    // Timing
    float total_time_ms = 0.0f;
    float medoid_time_ms = 0.0f;
};

// Information about a cluster representative
struct RepresentativeInfo {
    uint32_t cluster_id;                   // Cluster this represents
    uint32_t read_idx;                     // Global read index of representative
    float confidence;                      // Representative's confidence score
    size_t cluster_size;                   // Number of reads in cluster
    float avg_distance_to_members;         // Average distance to cluster members
    float centrality_score;                // Weighted degree within cluster
};

// Result of representative selection
struct ClusterRepResult {
    // Per-cluster representative info, indexed by cluster_id
    std::vector<RepresentativeInfo> representatives;

    // Quick lookup: representative read indices only
    std::vector<uint32_t> representative_reads;

    // Statistics
    ClusterRepStats stats;

    // Get representative for a specific cluster
    std::optional<RepresentativeInfo> GetRepresentative(uint32_t cluster_id) const;

    // Check if a read is a representative
    bool IsRepresentative(uint32_t read_idx) const;

    // Get all representative read indices as a sorted vector
    std::vector<uint32_t> GetRepresentativeReads() const;

    // Get cluster ID for a representative read
    std::optional<uint32_t> GetClusterForRepresentative(uint32_t read_idx) const;
};

// Cluster representative selector
//
// This is Stage 1 step 1.4 from the LLmap pipeline: after Self-WaveCollapse
// refines the clusters, select one representative read per cluster.
//
// The representative (medoid) minimizes the sum of distances to all other
// cluster members. This representative will be used in Stage 2 for projection
// to the reference genome.
//
// For efficiency:
// - Small clusters: exact medoid computation O(n²)
// - Large clusters: approximate via sampling
// - Very large clusters: fall back to max-degree heuristic
class ClusterRepSelector {
public:
    ClusterRepSelector();
    explicit ClusterRepSelector(const ClusterRepConfig& config);

    ~ClusterRepSelector();

    // Move-only
    ClusterRepSelector(ClusterRepSelector&&) noexcept;
    ClusterRepSelector& operator=(ClusterRepSelector&&) noexcept;
    ClusterRepSelector(const ClusterRepSelector&) = delete;
    ClusterRepSelector& operator=(const ClusterRepSelector&) = delete;

    // Select representatives from clustering result
    //
    // Uses similarity graph for distance computation
    std::unique_ptr<ClusterRepResult> Select(
        const SimilarityGraph& graph,
        const ClusteringResult& clustering);

    // Select representatives using Self-WaveCollapse results
    //
    // More accurate: uses confidence scores from EM refinement
    std::unique_ptr<ClusterRepResult> SelectWithConfidence(
        const SimilarityGraph& graph,
        const ClusteringResult& clustering,
        const SelfWaveCollapseResult& swc_result);

    // Select representative for a single cluster
    //
    // Returns read index of representative, or std::nullopt if cluster is empty
    std::optional<uint32_t> SelectForCluster(
        const SimilarityGraph& graph,
        std::span<const uint32_t> cluster_members);

    // Select representative for a single cluster with confidence scores
    std::optional<uint32_t> SelectForClusterWithConfidence(
        const SimilarityGraph& graph,
        std::span<const uint32_t> cluster_members,
        std::span<const float> confidence_scores);

    // Get/set configuration
    const ClusterRepConfig& GetConfig() const { return config_; }
    void SetConfig(const ClusterRepConfig& config) { config_ = config; }

private:
    ClusterRepConfig config_;

    // Internal state
    struct InternalState;
    std::unique_ptr<InternalState> state_;

    // Core medoid computation methods
    uint32_t ComputeExactMedoid(
        const SimilarityGraph& graph,
        std::span<const uint32_t> members);

    uint32_t ComputeApproxMedoid(
        const SimilarityGraph& graph,
        std::span<const uint32_t> members);

    uint32_t ComputeMaxDegree(
        const SimilarityGraph& graph,
        std::span<const uint32_t> members);

    uint32_t ComputeMaxConfidence(
        std::span<const uint32_t> members,
        std::span<const float> confidence_scores);

    // Distance computation
    float ComputeDistanceSum(
        const SimilarityGraph& graph,
        std::span<const uint32_t> members,
        uint32_t candidate);

    float ComputeCentrality(
        const SimilarityGraph& graph,
        std::span<const uint32_t> members,
        uint32_t candidate);
};

// Convenience function: select representatives with default config
std::unique_ptr<ClusterRepResult> SelectClusterRepresentatives(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering);

// Convenience function: select representatives using SWC results
std::unique_ptr<ClusterRepResult> SelectClusterRepresentatives(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering,
    const SelfWaveCollapseResult& swc_result);

}  // namespace llmap::self_interference
