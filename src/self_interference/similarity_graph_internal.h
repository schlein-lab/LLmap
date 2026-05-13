#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace llmap::self_interference::internal {

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

            std::sort(neighbors.begin(), neighbors.end(),
                [](const auto& a, const auto& b) {
                    return a.first < b.first;
                });

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

}  // namespace llmap::self_interference::internal
