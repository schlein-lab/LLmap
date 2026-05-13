#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace llmap::self_interference {

float LeidenClustering::CalculateModularity(
    const SimilarityGraph& graph,
    std::span<const uint32_t> labels,
    float resolution)
{
    const size_t n = graph.NumNodes();
    if (n == 0 || labels.size() != n) return 0.0f;

    // Compute total weight
    double total_weight = 0.0;
    std::vector<float> node_weights(n, 0.0f);

    for (size_t i = 0; i < n; ++i) {
        auto weights = graph.NeighborWeights(i);
        for (float w : weights) {
            node_weights[i] += w;
        }
        total_weight += node_weights[i];
    }
    total_weight /= 2.0;  // Each edge counted twice

    if (total_weight == 0.0) return 0.0f;

    // Compute community totals
    uint32_t max_label = *std::max_element(labels.begin(), labels.end());
    std::vector<double> community_weights(max_label + 1, 0.0);
    std::vector<double> community_internal(max_label + 1, 0.0);

    for (size_t i = 0; i < n; ++i) {
        community_weights[labels[i]] += node_weights[i];

        auto neighbors = graph.Neighbors(i);
        auto weights = graph.NeighborWeights(i);
        for (size_t j = 0; j < neighbors.size(); ++j) {
            if (labels[neighbors[j]] == labels[i]) {
                community_internal[labels[i]] += weights[j];
            }
        }
    }

    // Internal edges counted twice, so divide
    for (auto& ic : community_internal) {
        ic /= 2.0;
    }

    // Q = Sigma_c [e_c / m - gamma(Sigma_c / 2m)^2]
    double Q = 0.0;
    for (size_t c = 0; c <= max_label; ++c) {
        if (community_weights[c] == 0.0) continue;

        double e_c = community_internal[c] / total_weight;
        double a_c = community_weights[c] / (2.0 * total_weight);
        Q += e_c - resolution * a_c * a_c;
    }

    return static_cast<float>(Q);
}

bool LeidenClustering::VerifyWellConnected(
    const SimilarityGraph& graph,
    std::span<const uint32_t> labels)
{
    const size_t n = graph.NumNodes();
    if (n == 0) return true;

    // Group nodes by community
    std::unordered_map<uint32_t, std::vector<uint32_t>> communities;
    for (uint32_t i = 0; i < n; ++i) {
        communities[labels[i]].push_back(i);
    }

    // For each community, verify it's connected (all nodes reachable)
    for (const auto& [comm_id, members] : communities) {
        if (members.size() <= 1) continue;

        // BFS to check connectivity within community
        std::unordered_set<uint32_t> member_set(members.begin(), members.end());
        std::unordered_set<uint32_t> visited;
        std::vector<uint32_t> queue;

        queue.push_back(members[0]);
        visited.insert(members[0]);

        while (!queue.empty()) {
            uint32_t node = queue.back();
            queue.pop_back();

            auto neighbors = graph.Neighbors(node);
            for (uint32_t neighbor : neighbors) {
                if (member_set.count(neighbor) && !visited.count(neighbor)) {
                    visited.insert(neighbor);
                    queue.push_back(neighbor);
                }
            }
        }

        // If not all members visited, community is disconnected
        if (visited.size() != members.size()) {
            return false;
        }
    }

    return true;
}

// ClusteringResult methods

std::vector<std::vector<uint32_t>> ClusteringResult::GetCommunityMembers() const {
    if (num_communities == 0) return {};

    std::vector<std::vector<uint32_t>> members(num_communities);
    for (uint32_t i = 0; i < labels.size(); ++i) {
        if (labels[i] < num_communities) {
            members[labels[i]].push_back(i);
        }
    }
    return members;
}

std::vector<size_t> ClusteringResult::GetCommunitySizes() const {
    if (num_communities == 0) return {};

    std::vector<size_t> sizes(num_communities, 0);
    for (uint32_t label : labels) {
        if (label < num_communities) {
            ++sizes[label];
        }
    }
    return sizes;
}

ClusteringResult ClusteringResult::FilterByCommunitySize(
    size_t min_size,
    size_t max_size) const
{
    auto sizes = GetCommunitySizes();

    // Find which communities to keep
    std::vector<bool> keep(num_communities, true);
    size_t new_num_communities = 0;

    for (size_t c = 0; c < num_communities; ++c) {
        if (sizes[c] < min_size || (max_size > 0 && sizes[c] > max_size)) {
            keep[c] = false;
        } else {
            ++new_num_communities;
        }
    }

    // Build remapping
    std::vector<uint32_t> remap(num_communities, UINT32_MAX);
    uint32_t next_id = 0;
    for (size_t c = 0; c < num_communities; ++c) {
        if (keep[c]) {
            remap[c] = next_id++;
        }
    }

    // Create filtered result
    ClusteringResult result;
    result.labels.resize(labels.size());
    result.num_communities = new_num_communities;
    result.modularity = modularity;  // Approximate; should recompute
    result.stats = stats;
    result.stats.num_communities = new_num_communities;

    // Special community ID for filtered-out nodes
    const uint32_t FILTERED = new_num_communities;
    result.num_communities = new_num_communities + 1;  // +1 for filtered

    for (size_t i = 0; i < labels.size(); ++i) {
        if (remap[labels[i]] != UINT32_MAX) {
            result.labels[i] = remap[labels[i]];
        } else {
            result.labels[i] = FILTERED;
        }
    }

    return result;
}

void ClusteringResult::RelabelContiguous() {
    if (labels.empty()) return;

    // Find all unique labels and create mapping to contiguous IDs
    std::unordered_map<uint32_t, uint32_t> label_map;
    uint32_t next_id = 0;

    for (uint32_t& label : labels) {
        auto it = label_map.find(label);
        if (it == label_map.end()) {
            label_map[label] = next_id;
            label = next_id++;
        } else {
            label = it->second;
        }
    }

    num_communities = next_id;
}

// Convenience functions

std::unique_ptr<ClusteringResult> RunLeiden(
    const SimilarityGraph& graph,
    float resolution)
{
    LeidenConfig config;
    config.resolution = resolution;

    LeidenClustering leiden(config);
    return leiden.Cluster(graph);
}

std::vector<uint32_t> GetCommunityLabels(
    const SimilarityGraph& graph,
    float resolution)
{
    auto result = RunLeiden(graph, resolution);
    if (result) {
        return std::move(result->labels);
    }
    return {};
}

}  // namespace llmap::self_interference
