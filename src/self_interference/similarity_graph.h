#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace llmap::self_interference {

// Forward declaration
struct BatchKnnResult;

// Edge in the similarity graph
struct Edge {
    uint32_t source;      // Source node index
    uint32_t target;      // Target node index
    float weight;         // Similarity weight (1 - normalized_distance, or raw similarity)

    bool operator==(const Edge& other) const {
        return source == other.source && target == other.target;
    }
};

// Statistics about the similarity graph
struct SimilarityGraphStats {
    size_t num_nodes = 0;
    size_t num_edges = 0;
    size_t num_self_loops_removed = 0;
    size_t num_duplicates_merged = 0;
    float avg_degree = 0.0f;
    float max_weight = 0.0f;
    float min_weight = 0.0f;
    float build_time_ms = 0.0f;
};

// Configuration for building the similarity graph
struct SimilarityGraphConfig {
    // Distance to weight conversion
    bool convert_distance_to_similarity = true;  // weight = 1 / (1 + distance)
    float distance_scale = 1.0f;                 // Scale factor for distance normalization

    // Edge filtering
    float min_weight_threshold = 0.0f;           // Drop edges below this weight
    float max_distance_threshold = -1.0f;        // Drop edges above this distance (-1 = no limit)
    size_t max_edges_per_node = 0;               // 0 = keep all from k-NN

    // Graph structure
    bool make_symmetric = true;                  // Add reverse edges for undirected graph
    bool remove_self_loops = true;               // Remove edges where source == target
    bool merge_duplicates = true;                // Merge duplicate edges (keep max weight)

    // Memory hints
    size_t expected_num_nodes = 0;               // For pre-allocation
};

// Sparse similarity graph in CSR (Compressed Sparse Row) format
// Optimized for Leiden clustering and graph algorithms
//
// CSR Layout:
//   row_offsets[i] = start index in col_indices/weights for node i
//   row_offsets[num_nodes] = total edges (sentinel)
//   col_indices[j] = target node for edge j
//   weights[j] = weight for edge j
class SimilarityGraph {
public:
    SimilarityGraph() = default;

    // Build from batch k-NN results
    // indices: [num_queries * k] neighbor indices from FAISS
    // distances: [num_queries * k] distances from FAISS
    // num_queries: number of query nodes
    // k: number of neighbors per query
    static std::unique_ptr<SimilarityGraph> BuildFromKnn(
        std::span<const int64_t> indices,
        std::span<const float> distances,
        size_t num_queries,
        size_t k,
        const SimilarityGraphConfig& config = {});

    // Build from edge list
    static std::unique_ptr<SimilarityGraph> BuildFromEdgeList(
        std::span<const Edge> edges,
        size_t num_nodes,
        const SimilarityGraphConfig& config = {});

    // Build from BatchKnnResult directly
    static std::unique_ptr<SimilarityGraph> BuildFromBatchKnnResult(
        const BatchKnnResult& result,
        const SimilarityGraphConfig& config = {});

    ~SimilarityGraph() = default;

    // Non-copyable, movable
    SimilarityGraph(const SimilarityGraph&) = delete;
    SimilarityGraph& operator=(const SimilarityGraph&) = delete;
    SimilarityGraph(SimilarityGraph&&) noexcept = default;
    SimilarityGraph& operator=(SimilarityGraph&&) noexcept = default;

    // ========== Accessors ==========

    // Number of nodes in the graph
    size_t NumNodes() const noexcept { return num_nodes_; }

    // Number of edges (directed count)
    size_t NumEdges() const noexcept { return col_indices_.size(); }

    // Get neighbors of a node
    std::span<const uint32_t> Neighbors(size_t node) const;
    std::span<const float> NeighborWeights(size_t node) const;

    // Get degree of a node
    size_t Degree(size_t node) const;

    // Check if edge exists
    bool HasEdge(uint32_t source, uint32_t target) const;

    // Get edge weight (returns 0 if edge doesn't exist)
    float GetEdgeWeight(uint32_t source, uint32_t target) const;

    // ========== Raw CSR Access ==========

    std::span<const uint32_t> RowOffsets() const noexcept {
        return row_offsets_;
    }

    std::span<const uint32_t> ColIndices() const noexcept {
        return col_indices_;
    }

    std::span<const float> Weights() const noexcept {
        return weights_;
    }

    // ========== Graph Properties ==========

    bool IsSymmetric() const noexcept { return is_symmetric_; }
    bool IsEmpty() const noexcept { return num_nodes_ == 0 || col_indices_.empty(); }

    // Statistics
    SimilarityGraphStats GetStats() const;
    float AverageDegree() const;
    size_t MaxDegree() const;

    // ========== Export ==========

    // Convert to edge list (for debugging or other algorithms)
    std::vector<Edge> ToEdgeList() const;

    // Export to simple text format (source target weight per line)
    std::string ToEdgeListString() const;

    // ========== Serialization ==========

    bool Save(const std::string& path) const;
    static std::unique_ptr<SimilarityGraph> Load(const std::string& path);

private:
    size_t num_nodes_ = 0;
    bool is_symmetric_ = false;

    // CSR arrays
    std::vector<uint32_t> row_offsets_;  // size = num_nodes + 1
    std::vector<uint32_t> col_indices_;  // size = num_edges
    std::vector<float> weights_;         // size = num_edges

    // Build stats (for debugging/profiling)
    SimilarityGraphStats build_stats_;
};

// ========== Utility Functions ==========

// Convert L2 distance to similarity weight
// Default: weight = 1 / (1 + distance * scale)
inline float DistanceToSimilarity(float distance, float scale = 1.0f) {
    return 1.0f / (1.0f + distance * scale);
}

// Convert inner product similarity to distance
// IP distance = 1 - IP (for normalized vectors, IP in [0,1])
inline float InnerProductToDistance(float ip) {
    return 1.0f - ip;
}

}  // namespace llmap::self_interference
