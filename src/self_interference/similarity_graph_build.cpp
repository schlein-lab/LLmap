#include "self_interference/similarity_graph.h"
#include "self_interference/similarity_graph_internal.h"
#include "self_interference/faiss_wrapper.h"

#include <algorithm>
#include <chrono>

namespace llmap::self_interference {

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

    size_t num_nodes = config.expected_num_nodes;
    if (num_nodes == 0) {
        int64_t max_idx = -1;
        for (size_t i = 0; i < num_queries * k; ++i) {
            if (indices[i] >= 0) {
                max_idx = std::max(max_idx, indices[i]);
            }
        }
        num_nodes = std::max(static_cast<size_t>(max_idx + 1), num_queries);
    }

    graph->num_nodes_ = num_nodes;

    internal::AdjacencyBuilder builder(num_nodes);

    size_t self_loops_removed = 0;
    size_t edges_added = 0;

    for (size_t q = 0; q < num_queries; ++q) {
        const size_t base = q * k;
        size_t edges_for_node = 0;

        for (size_t i = 0; i < k; ++i) {
            const int64_t neighbor_idx = indices[base + i];
            const float dist = distances[base + i];

            if (neighbor_idx < 0) continue;

            const auto source = static_cast<uint32_t>(q);
            const auto target = static_cast<uint32_t>(neighbor_idx);

            if (config.remove_self_loops && source == target) {
                ++self_loops_removed;
                continue;
            }

            if (config.max_distance_threshold >= 0 && dist > config.max_distance_threshold) {
                continue;
            }

            float weight = dist;
            if (config.convert_distance_to_similarity) {
                weight = DistanceToSimilarity(dist, config.distance_scale);
            }

            if (weight < config.min_weight_threshold) {
                continue;
            }

            if (config.max_edges_per_node > 0 && edges_for_node >= config.max_edges_per_node) {
                continue;
            }

            builder.AddEdge(source, target, weight);
            ++edges_for_node;
            ++edges_added;

            if (config.make_symmetric && source != target) {
                builder.AddEdge(target, source, weight);
            }
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

std::unique_ptr<SimilarityGraph> SimilarityGraph::BuildFromEdgeList(
    std::span<const Edge> edges,
    size_t num_nodes,
    const SimilarityGraphConfig& config)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto graph = std::unique_ptr<SimilarityGraph>(new SimilarityGraph());
    graph->num_nodes_ = num_nodes;

    internal::AdjacencyBuilder builder(num_nodes);

    size_t self_loops_removed = 0;

    for (const auto& edge : edges) {
        if (config.remove_self_loops && edge.source == edge.target) {
            ++self_loops_removed;
            continue;
        }

        if (edge.weight < config.min_weight_threshold) {
            continue;
        }

        builder.AddEdge(edge.source, edge.target, edge.weight);

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

}  // namespace llmap::self_interference
