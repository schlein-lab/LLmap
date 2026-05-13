#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace llmap::self_interference {

class SimilarityGraph;
struct ClusteringResult;

// Configuration for Self-WaveCollapse intra-cluster EM
struct SelfWaveCollapseConfig {
    // EM convergence control
    size_t max_iterations = 20;
    float convergence_threshold = 1e-5f;  // Stop when max probability change < this
    float collapse_threshold = 0.95f;     // Read collapses when max P >= this

    // Damping factor: P_new = (1 - gamma) * P_old + gamma * update
    float gamma = 0.3f;

    // Similarity weighting
    float similarity_exponent = 2.0f;     // Higher = sharper similarity weighting

    // Cluster handling
    size_t min_cluster_size = 2;          // Skip EM for clusters smaller than this
    size_t max_cluster_size = 0;          // 0 = no limit; split larger clusters

    // Sub-clustering for large clusters
    bool enable_subclustering = true;
    size_t subcluster_threshold = 1000;   // Subcluster if cluster size > this

    // Memory/performance
    bool use_gpu = false;
    size_t batch_size = 10000;

    // Random seed for tie-breaking
    uint64_t seed = 42;
};

// Statistics from a Self-WaveCollapse run
struct SelfWaveCollapseStats {
    size_t num_reads = 0;
    size_t num_clusters = 0;

    // Convergence stats
    size_t total_iterations = 0;
    size_t clusters_converged = 0;
    size_t clusters_reached_max_iter = 0;
    float avg_iterations_per_cluster = 0.0f;

    // Read stats
    size_t reads_collapsed = 0;           // Reads that concentrated to single position
    size_t reads_split = 0;               // Reads that moved to different clusters
    size_t reads_unchanged = 0;           // Reads staying in original cluster

    // Timing
    float total_time_ms = 0.0f;
    float em_time_ms = 0.0f;
    float subclustering_time_ms = 0.0f;

    // Quality metrics
    float avg_final_entropy = 0.0f;       // Lower = more concentrated
    float avg_max_probability = 0.0f;     // Higher = more confident
};

// Per-read state after Self-WaveCollapse
struct ReadAssignment {
    uint32_t read_idx;                    // Original read index
    uint32_t cluster_id;                  // Final cluster assignment
    float confidence;                     // Max probability in final distribution
    bool collapsed;                       // True if P_max >= collapse_threshold
    uint32_t anchor_read;                 // Representative read in cluster (highest centrality)
};

// Result of Self-WaveCollapse
struct SelfWaveCollapseResult {
    // Per-read assignments (indexed by read_idx)
    std::vector<ReadAssignment> assignments;

    // Final cluster count (may differ from input if splits occurred)
    size_t num_clusters = 0;

    // Statistics
    SelfWaveCollapseStats stats;

    // Get reads in a specific cluster
    std::vector<uint32_t> GetClusterMembers(uint32_t cluster_id) const;

    // Get cluster sizes
    std::vector<size_t> GetClusterSizes() const;

    // Get anchor (representative) for each cluster
    std::vector<uint32_t> GetClusterAnchors() const;
};

// Self-WaveCollapse: intra-cluster EM refinement
//
// This is Stage 1 step 1.3 from the LLmap pipeline. Within each cluster
// produced by Leiden, reads inform each other's "position" through an
// EM-style iterative refinement.
//
// The algorithm treats cluster membership as a probability distribution:
// - P(r, c) = probability that read r belongs to sub-position c within cluster
// - Initially, c = other reads (each read is a potential "anchor")
// - Update: P(r, c) += similarity(r, c) * sum_j(P(j, c) * similarity(r, j))
// - Converged reads "collapse" to their most likely anchor
//
// This achieves read-to-read information sharing: reads with high mutual
// similarity reinforce each other, while outliers get pushed to the periphery.
class SelfWaveCollapse {
public:
    SelfWaveCollapse();
    explicit SelfWaveCollapse(const SelfWaveCollapseConfig& config);

    ~SelfWaveCollapse();

    // Move-only
    SelfWaveCollapse(SelfWaveCollapse&&) noexcept;
    SelfWaveCollapse& operator=(SelfWaveCollapse&&) noexcept;
    SelfWaveCollapse(const SelfWaveCollapse&) = delete;
    SelfWaveCollapse& operator=(const SelfWaveCollapse&) = delete;

    // Run Self-WaveCollapse on clustered reads
    //
    // Parameters:
    //   graph: similarity graph over reads (from FAISS k-NN)
    //   clustering: Leiden clustering result (defines initial clusters)
    //
    // Returns: refined assignments with per-read confidence and anchors
    std::unique_ptr<SelfWaveCollapseResult> Refine(
        const SimilarityGraph& graph,
        const ClusteringResult& clustering);

    // Run on a single cluster (useful for parallel processing)
    //
    // Parameters:
    //   graph: full similarity graph
    //   cluster_members: read indices in this cluster
    //   cluster_id: ID for this cluster in output
    //
    // Returns: assignments for reads in this cluster only
    std::vector<ReadAssignment> RefineCluster(
        const SimilarityGraph& graph,
        std::span<const uint32_t> cluster_members,
        uint32_t cluster_id);

    // Get/set configuration
    const SelfWaveCollapseConfig& GetConfig() const { return config_; }
    void SetConfig(const SelfWaveCollapseConfig& config) { config_ = config; }

private:
    SelfWaveCollapseConfig config_;

    // Internal state
    struct InternalState;
    std::unique_ptr<InternalState> state_;

    // Core EM iteration for a single cluster
    void RunEM(
        InternalState& state,
        const SimilarityGraph& graph,
        std::span<const uint32_t> members);

    // Compute similarity weight between two reads
    float ComputeSimilarityWeight(
        const SimilarityGraph& graph,
        uint32_t read_a,
        uint32_t read_b) const;

    // Find anchor (most central) read in cluster
    uint32_t FindAnchor(
        const SimilarityGraph& graph,
        std::span<const uint32_t> members) const;

    // Check convergence of probability matrix
    bool CheckConvergence(
        std::span<const float> old_probs,
        std::span<const float> new_probs) const;

    // Compute entropy of probability distribution
    static float ComputeEntropy(std::span<const float> probs);
};

// Convenience function: run Self-WaveCollapse with default config
std::unique_ptr<SelfWaveCollapseResult> RunSelfWaveCollapse(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering);

}  // namespace llmap::self_interference
