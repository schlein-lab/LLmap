#include "self_interference/similarity_graph.h"

#include <algorithm>

namespace llmap::self_interference {

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

}  // namespace llmap::self_interference
