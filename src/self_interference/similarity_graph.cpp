#include "self_interference/similarity_graph.h"
#include "self_interference/faiss_wrapper.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace llmap::self_interference {

namespace {

// Hash for Edge (used in duplicate detection)
struct EdgeHash {
    size_t operator()(const Edge& e) const {
        return std::hash<uint64_t>{}(
            (static_cast<uint64_t>(e.source) << 32) | e.target);
    }
};

// Build adjacency list from edges, then convert to CSR
struct AdjacencyBuilder {
    size_t num_nodes = 0;
    std::vector<std::vector<std::pair<uint32_t, float>>> adj;

    explicit AdjacencyBuilder(size_t n) : num_nodes(n), adj(n) {}

    void AddEdge(uint32_t src, uint32_t dst, float weight) {
        if (src < num_nodes && dst < num_nodes) {
            adj[src].emplace_back(dst, weight);
        }
    }

    void SortAndDeduplicate(bool keep_max_weight) {
        for (auto& neighbors : adj) {
            if (neighbors.empty()) continue;

            // Sort by target
            std::sort(neighbors.begin(), neighbors.end(),
                [](const auto& a, const auto& b) {
                    return a.first < b.first;
                });

            // Deduplicate (keep max weight if specified)
            auto write = neighbors.begin();
            for (auto read = neighbors.begin(); read != neighbors.end(); ) {
                uint32_t target = read->first;
                float best_weight = read->second;

                ++read;
                while (read != neighbors.end() && read->first == target) {
                    if (keep_max_weight) {
                        best_weight = std::max(best_weight, read->second);
                    }
                    ++read;
                }

                *write = {target, best_weight};
                ++write;
            }
            neighbors.erase(write, neighbors.end());
        }
    }

    void ToCsr(
        std::vector<uint32_t>& row_offsets,
        std::vector<uint32_t>& col_indices,
        std::vector<float>& weights) const
    {
        row_offsets.clear();
        col_indices.clear();
        weights.clear();

        row_offsets.reserve(num_nodes + 1);

        size_t total_edges = 0;
        for (const auto& neighbors : adj) {
            total_edges += neighbors.size();
        }
        col_indices.reserve(total_edges);
        weights.reserve(total_edges);

        row_offsets.push_back(0);
        for (const auto& neighbors : adj) {
            for (const auto& [target, weight] : neighbors) {
                col_indices.push_back(target);
                weights.push_back(weight);
            }
            row_offsets.push_back(static_cast<uint32_t>(col_indices.size()));
        }
    }
};

}  // namespace

std::unique_ptr<SimilarityGraph> SimilarityGraph::BuildFromKnn(
    std::span<const int64_t> indices,
    std::span<const float> distances,
    size_t num_queries,
    size_t k,
    const SimilarityGraphConfig& config)
{
    auto start = std::chrono::high_resolution_clock::now();

    if (indices.size() < num_queries * k || distances.size() < num_queries * k) {
        return nullptr;
    }

    auto graph = std::unique_ptr<SimilarityGraph>(new SimilarityGraph());

    // Determine number of nodes (max index + 1, or from config)
    size_t num_nodes = config.expected_num_nodes;
    if (num_nodes == 0) {
        int64_t max_idx = -1;
        for (size_t i = 0; i < num_queries * k; ++i) {
            if (indices[i] >= 0) {
                max_idx = std::max(max_idx, indices[i]);
            }
        }
        // Also consider query indices
        num_nodes = std::max(static_cast<size_t>(max_idx + 1), num_queries);
    }

    graph->num_nodes_ = num_nodes;

    // Build adjacency list
    AdjacencyBuilder builder(num_nodes);

    size_t self_loops_removed = 0;
    size_t edges_added = 0;

    for (size_t q = 0; q < num_queries; ++q) {
        const size_t base = q * k;
        size_t edges_for_node = 0;

        for (size_t i = 0; i < k; ++i) {
            const int64_t neighbor_idx = indices[base + i];
            const float dist = distances[base + i];

            // Skip invalid indices
            if (neighbor_idx < 0) continue;

            const auto source = static_cast<uint32_t>(q);
            const auto target = static_cast<uint32_t>(neighbor_idx);

            // Self-loop removal
            if (config.remove_self_loops && source == target) {
                ++self_loops_removed;
                continue;
            }

            // Distance threshold
            if (config.max_distance_threshold >= 0 && dist > config.max_distance_threshold) {
                continue;
            }

            // Convert distance to weight
            float weight = dist;
            if (config.convert_distance_to_similarity) {
                weight = DistanceToSimilarity(dist, config.distance_scale);
            }

            // Weight threshold
            if (weight < config.min_weight_threshold) {
                continue;
            }

            // Max edges per node
            if (config.max_edges_per_node > 0 && edges_for_node >= config.max_edges_per_node) {
                continue;
            }

            builder.AddEdge(source, target, weight);
            ++edges_for_node;
            ++edges_added;

            // Add reverse edge for symmetric graph
            if (config.make_symmetric && source != target) {
                builder.AddEdge(target, source, weight);
            }
        }
    }

    // Sort and deduplicate
    builder.SortAndDeduplicate(config.merge_duplicates);

    // Convert to CSR
    builder.ToCsr(graph->row_offsets_, graph->col_indices_, graph->weights_);

    graph->is_symmetric_ = config.make_symmetric;

    auto end = std::chrono::high_resolution_clock::now();

    // Build stats
    graph->build_stats_.num_nodes = graph->num_nodes_;
    graph->build_stats_.num_edges = graph->col_indices_.size();
    graph->build_stats_.num_self_loops_removed = self_loops_removed;
    graph->build_stats_.build_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    if (!graph->weights_.empty()) {
        graph->build_stats_.max_weight = *std::max_element(
            graph->weights_.begin(), graph->weights_.end());
        graph->build_stats_.min_weight = *std::min_element(
            graph->weights_.begin(), graph->weights_.end());
    }

    if (graph->num_nodes_ > 0) {
        graph->build_stats_.avg_degree =
            static_cast<float>(graph->col_indices_.size()) /
            static_cast<float>(graph->num_nodes_);
    }

    return graph;
}

std::unique_ptr<SimilarityGraph> SimilarityGraph::BuildFromEdgeList(
    std::span<const Edge> edges,
    size_t num_nodes,
    const SimilarityGraphConfig& config)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto graph = std::unique_ptr<SimilarityGraph>(new SimilarityGraph());
    graph->num_nodes_ = num_nodes;

    AdjacencyBuilder builder(num_nodes);

    size_t self_loops_removed = 0;

    for (const auto& edge : edges) {
        // Self-loop removal
        if (config.remove_self_loops && edge.source == edge.target) {
            ++self_loops_removed;
            continue;
        }

        // Weight threshold
        if (edge.weight < config.min_weight_threshold) {
            continue;
        }

        builder.AddEdge(edge.source, edge.target, edge.weight);

        // Add reverse edge for symmetric graph
        if (config.make_symmetric && edge.source != edge.target) {
            builder.AddEdge(edge.target, edge.source, edge.weight);
        }
    }

    builder.SortAndDeduplicate(config.merge_duplicates);
    builder.ToCsr(graph->row_offsets_, graph->col_indices_, graph->weights_);

    graph->is_symmetric_ = config.make_symmetric;

    auto end = std::chrono::high_resolution_clock::now();

    graph->build_stats_.num_nodes = graph->num_nodes_;
    graph->build_stats_.num_edges = graph->col_indices_.size();
    graph->build_stats_.num_self_loops_removed = self_loops_removed;
    graph->build_stats_.build_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    if (!graph->weights_.empty()) {
        graph->build_stats_.max_weight = *std::max_element(
            graph->weights_.begin(), graph->weights_.end());
        graph->build_stats_.min_weight = *std::min_element(
            graph->weights_.begin(), graph->weights_.end());
    }

    if (graph->num_nodes_ > 0) {
        graph->build_stats_.avg_degree =
            static_cast<float>(graph->col_indices_.size()) /
            static_cast<float>(graph->num_nodes_);
    }

    return graph;
}

std::unique_ptr<SimilarityGraph> SimilarityGraph::BuildFromBatchKnnResult(
    const BatchKnnResult& result,
    const SimilarityGraphConfig& config)
{
    return BuildFromKnn(
        result.indices,
        result.distances,
        result.num_queries,
        result.k,
        config);
}

std::span<const uint32_t> SimilarityGraph::Neighbors(size_t node) const {
    if (node >= num_nodes_) {
        return {};
    }
    const size_t start = row_offsets_[node];
    const size_t end = row_offsets_[node + 1];
    return std::span<const uint32_t>(col_indices_.data() + start, end - start);
}

std::span<const float> SimilarityGraph::NeighborWeights(size_t node) const {
    if (node >= num_nodes_) {
        return {};
    }
    const size_t start = row_offsets_[node];
    const size_t end = row_offsets_[node + 1];
    return std::span<const float>(weights_.data() + start, end - start);
}

size_t SimilarityGraph::Degree(size_t node) const {
    if (node >= num_nodes_) {
        return 0;
    }
    return row_offsets_[node + 1] - row_offsets_[node];
}

bool SimilarityGraph::HasEdge(uint32_t source, uint32_t target) const {
    if (source >= num_nodes_) {
        return false;
    }
    auto neighbors = Neighbors(source);
    return std::binary_search(neighbors.begin(), neighbors.end(), target);
}

float SimilarityGraph::GetEdgeWeight(uint32_t source, uint32_t target) const {
    if (source >= num_nodes_) {
        return 0.0f;
    }
    const size_t start = row_offsets_[source];
    const size_t end = row_offsets_[source + 1];

    auto it = std::lower_bound(
        col_indices_.begin() + start,
        col_indices_.begin() + end,
        target);

    if (it != col_indices_.begin() + end && *it == target) {
        return weights_[it - col_indices_.begin()];
    }
    return 0.0f;
}

SimilarityGraphStats SimilarityGraph::GetStats() const {
    return build_stats_;
}

float SimilarityGraph::AverageDegree() const {
    if (num_nodes_ == 0) return 0.0f;
    return static_cast<float>(col_indices_.size()) / static_cast<float>(num_nodes_);
}

size_t SimilarityGraph::MaxDegree() const {
    if (num_nodes_ == 0) return 0;
    size_t max_deg = 0;
    for (size_t i = 0; i < num_nodes_; ++i) {
        const size_t deg = row_offsets_[i + 1] - row_offsets_[i];
        max_deg = std::max(max_deg, deg);
    }
    return max_deg;
}

std::vector<Edge> SimilarityGraph::ToEdgeList() const {
    std::vector<Edge> edges;
    edges.reserve(col_indices_.size());

    for (size_t src = 0; src < num_nodes_; ++src) {
        const size_t start = row_offsets_[src];
        const size_t end = row_offsets_[src + 1];
        for (size_t i = start; i < end; ++i) {
            edges.push_back({
                static_cast<uint32_t>(src),
                col_indices_[i],
                weights_[i]
            });
        }
    }

    return edges;
}

std::string SimilarityGraph::ToEdgeListString() const {
    std::ostringstream oss;
    oss << "# SimilarityGraph: " << num_nodes_ << " nodes, "
        << col_indices_.size() << " edges\n";
    oss << "# source target weight\n";

    for (size_t src = 0; src < num_nodes_; ++src) {
        const size_t start = row_offsets_[src];
        const size_t end = row_offsets_[src + 1];
        for (size_t i = start; i < end; ++i) {
            oss << src << " " << col_indices_[i] << " " << weights_[i] << "\n";
        }
    }

    return oss.str();
}

bool SimilarityGraph::Save(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    // Header
    const char magic[] = "LLMG";  // LLMap Graph
    file.write(magic, 4);

    const uint32_t version = 1;
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));

    const auto n = static_cast<uint64_t>(num_nodes_);
    const auto e = static_cast<uint64_t>(col_indices_.size());
    file.write(reinterpret_cast<const char*>(&n), sizeof(n));
    file.write(reinterpret_cast<const char*>(&e), sizeof(e));

    const uint8_t sym = is_symmetric_ ? 1 : 0;
    file.write(reinterpret_cast<const char*>(&sym), sizeof(sym));

    // Arrays
    file.write(reinterpret_cast<const char*>(row_offsets_.data()),
               row_offsets_.size() * sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(col_indices_.data()),
               col_indices_.size() * sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(weights_.data()),
               weights_.size() * sizeof(float));

    return file.good();
}

std::unique_ptr<SimilarityGraph> SimilarityGraph::Load(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return nullptr;

    // Header
    char magic[4];
    file.read(magic, 4);
    if (std::memcmp(magic, "LLMG", 4) != 0) return nullptr;

    uint32_t version;
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 1) return nullptr;

    uint64_t n, e;
    file.read(reinterpret_cast<char*>(&n), sizeof(n));
    file.read(reinterpret_cast<char*>(&e), sizeof(e));

    uint8_t sym;
    file.read(reinterpret_cast<char*>(&sym), sizeof(sym));

    auto graph = std::unique_ptr<SimilarityGraph>(new SimilarityGraph());
    graph->num_nodes_ = static_cast<size_t>(n);
    graph->is_symmetric_ = (sym != 0);

    // Arrays
    graph->row_offsets_.resize(n + 1);
    graph->col_indices_.resize(e);
    graph->weights_.resize(e);

    file.read(reinterpret_cast<char*>(graph->row_offsets_.data()),
              graph->row_offsets_.size() * sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(graph->col_indices_.data()),
              graph->col_indices_.size() * sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(graph->weights_.data()),
              graph->weights_.size() * sizeof(float));

    if (!file.good()) return nullptr;

    // Rebuild stats
    graph->build_stats_.num_nodes = graph->num_nodes_;
    graph->build_stats_.num_edges = graph->col_indices_.size();
    if (!graph->weights_.empty()) {
        graph->build_stats_.max_weight = *std::max_element(
            graph->weights_.begin(), graph->weights_.end());
        graph->build_stats_.min_weight = *std::min_element(
            graph->weights_.begin(), graph->weights_.end());
    }
    if (graph->num_nodes_ > 0) {
        graph->build_stats_.avg_degree =
            static_cast<float>(graph->col_indices_.size()) /
            static_cast<float>(graph->num_nodes_);
    }

    return graph;
}

}  // namespace llmap::self_interference
