#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>

namespace llmap::self_interference {

std::optional<RepresentativeInfo> ClusterRepResult::GetRepresentative(
    uint32_t cluster_id) const
{
    for (const auto& rep : representatives) {
        if (rep.cluster_id == cluster_id) {
            return rep;
        }
    }
    return std::nullopt;
}

bool ClusterRepResult::IsRepresentative(uint32_t read_idx) const {
    return std::find(representative_reads.begin(),
                     representative_reads.end(),
                     read_idx) != representative_reads.end();
}

std::vector<uint32_t> ClusterRepResult::GetRepresentativeReads() const {
    std::vector<uint32_t> sorted_reps = representative_reads;
    std::sort(sorted_reps.begin(), sorted_reps.end());
    return sorted_reps;
}

std::optional<uint32_t> ClusterRepResult::GetClusterForRepresentative(
    uint32_t read_idx) const
{
    for (const auto& rep : representatives) {
        if (rep.read_idx == read_idx) {
            return rep.cluster_id;
        }
    }
    return std::nullopt;
}

std::unique_ptr<ClusterRepResult> SelectClusterRepresentatives(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering)
{
    ClusterRepSelector selector;
    return selector.Select(graph, clustering);
}

std::unique_ptr<ClusterRepResult> SelectClusterRepresentatives(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering,
    const SelfWaveCollapseResult& swc_result)
{
    ClusterRepSelector selector;
    return selector.SelectWithConfidence(graph, clustering, swc_result);
}

}  // namespace llmap::self_interference
