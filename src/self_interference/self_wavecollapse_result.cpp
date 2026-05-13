#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"
#include "self_interference/leiden_clustering.h"

#include <unordered_map>

namespace llmap::self_interference {

std::vector<uint32_t> SelfWaveCollapseResult::GetClusterMembers(uint32_t cluster_id) const {
    std::vector<uint32_t> members;
    for (const auto& a : assignments) {
        if (a.cluster_id == cluster_id) {
            members.push_back(a.read_idx);
        }
    }
    return members;
}

std::vector<size_t> SelfWaveCollapseResult::GetClusterSizes() const {
    if (num_clusters == 0) return {};

    std::vector<size_t> sizes(num_clusters, 0);
    for (const auto& a : assignments) {
        if (a.cluster_id < num_clusters) {
            ++sizes[a.cluster_id];
        }
    }
    return sizes;
}

std::vector<uint32_t> SelfWaveCollapseResult::GetClusterAnchors() const {
    if (num_clusters == 0) return {};

    // For each cluster, find the most common anchor (by vote)
    // or the anchor with highest total confidence
    std::vector<std::unordered_map<uint32_t, float>> anchor_scores(num_clusters);

    for (const auto& a : assignments) {
        if (a.cluster_id < num_clusters) {
            anchor_scores[a.cluster_id][a.anchor_read] += a.confidence;
        }
    }

    std::vector<uint32_t> anchors(num_clusters, 0);
    for (size_t c = 0; c < num_clusters; ++c) {
        float best_score = 0.0f;
        for (const auto& [anchor, score] : anchor_scores[c]) {
            if (score > best_score) {
                best_score = score;
                anchors[c] = anchor;
            }
        }
    }

    return anchors;
}

// Convenience function: run Self-WaveCollapse with default config
std::unique_ptr<SelfWaveCollapseResult> RunSelfWaveCollapse(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering)
{
    SelfWaveCollapse swc;
    return swc.Refine(graph, clustering);
}

}  // namespace llmap::self_interference
