#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace llmap::self_interference {

class SimilarityGraph;

// Configuration for Leiden clustering
struct LeidenConfig {
    // Resolution parameter: higher = more/smaller communities
    float resolution = 1.0f;

    // Convergence control
    size_t max_iterations = 100;
    float min_modularity_gain = 1e-7f;

    // Refinement settings
    bool enable_refinement = true;
    size_t refinement_iterations = 10;

    // Randomization
    uint64_t seed = 42;
    bool randomize_node_order = true;

    // Community constraints
    size_t min_community_size = 1;      // Drop communities smaller than this
    size_t max_community_size = 0;      // 0 = no limit

    // Memory/performance hints
    bool use_gpu = false;               // When CUDA available
    size_t batch_size = 10000;          // For GPU batching
};

// Statistics from a Leiden clustering run
struct LeidenStats {
    size_t num_nodes = 0;
    size_t num_communities = 0;
    size_t num_iterations = 0;
    float final_modularity = 0.0f;
    float time_ms = 0.0f;

    // Per-phase timing
    float local_move_time_ms = 0.0f;
    float refinement_time_ms = 0.0f;
    float aggregation_time_ms = 0.0f;

    // Quality metrics
    float avg_community_size = 0.0f;
    size_t min_community_size = 0;
    size_t max_community_size = 0;
    size_t singleton_count = 0;  // Communities with 1 node
};

// Result of community detection
struct ClusteringResult {
    // Community assignment per node: labels[i] = community ID for node i
    std::vector<uint32_t> labels;

    // Number of communities found
    size_t num_communities = 0;

    // Quality score (modularity)
    float modularity = 0.0f;

    // Statistics
    LeidenStats stats;

    // Get community members
    std::vector<std::vector<uint32_t>> GetCommunityMembers() const;

    // Get size of each community
    std::vector<size_t> GetCommunitySizes() const;

    // Filter: keep only communities meeting size criteria
    ClusteringResult FilterByCommunitySize(
        size_t min_size,
        size_t max_size = 0) const;

    // Relabel communities contiguously (0, 1, 2, ...)
    void RelabelContiguous();
};

// Leiden community detection algorithm
//
// Implements the Leiden algorithm (Traag et al. 2019) for finding communities
// in weighted graphs. Leiden improves on Louvain by guaranteeing well-connected
// communities through a refinement phase.
//
// Algorithm phases (each iteration):
// 1. Local move phase: nodes move to neighboring communities to increase modularity
// 2. Refinement phase: partition refined within each community
// 3. Aggregation phase: communities become nodes in a coarser graph
//
// References:
// - Traag, V.A., Waltman, L. & van Eck, N.J. From Louvain to Leiden:
//   guaranteeing well-connected communities. Sci Rep 9, 5233 (2019)
class LeidenClustering {
public:
    LeidenClustering();
    explicit LeidenClustering(const LeidenConfig& config);

    ~LeidenClustering();

    // Move-only (PIMPL idiom)
    LeidenClustering(LeidenClustering&&) noexcept;
    LeidenClustering& operator=(LeidenClustering&&) noexcept;
    LeidenClustering(const LeidenClustering&) = delete;
    LeidenClustering& operator=(const LeidenClustering&) = delete;

    // Perform community detection on a graph
    // Returns nullptr if clustering fails
    std::unique_ptr<ClusteringResult> Cluster(const SimilarityGraph& graph);

    // Cluster with pre-initialized labels (for hierarchical refinement)
    std::unique_ptr<ClusteringResult> ClusterWithInitialLabels(
        const SimilarityGraph& graph,
        std::span<const uint32_t> initial_labels);

    // Get/set configuration
    const LeidenConfig& GetConfig() const { return config_; }
    void SetConfig(const LeidenConfig& config) { config_ = config; }

    // Calculate modularity of a given partition
    static float CalculateModularity(
        const SimilarityGraph& graph,
        std::span<const uint32_t> labels,
        float resolution = 1.0f);

    // Check if communities are well-connected (Leiden guarantee)
    static bool VerifyWellConnected(
        const SimilarityGraph& graph,
        std::span<const uint32_t> labels);

private:
    LeidenConfig config_;

    // Internal state for iterative algorithm
    struct InternalState;
    std::unique_ptr<InternalState> state_;

    // Algorithm phases
    bool LocalMovePhase(InternalState& state);
    void RefinementPhase(InternalState& state);
    void AggregationPhase(InternalState& state);

    // Modularity change calculations
    float ComputeModularityGain(
        const InternalState& state,
        uint32_t node,
        uint32_t target_community) const;
};

// Convenience function: run Leiden with default config
std::unique_ptr<ClusteringResult> RunLeiden(
    const SimilarityGraph& graph,
    float resolution = 1.0f);

// Convenience function: run Leiden and get community labels directly
std::vector<uint32_t> GetCommunityLabels(
    const SimilarityGraph& graph,
    float resolution = 1.0f);

}  // namespace llmap::self_interference
