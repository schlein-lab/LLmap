#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"

#include <chrono>
#include <unordered_set>

namespace llmap::self_interference {

namespace {

class FastRng {
public:
    explicit FastRng(uint64_t seed) : state_(seed) {}

    uint64_t Next() {
        state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
        return state_;
    }

    uint32_t NextInRange(uint32_t max_exclusive) {
        if (max_exclusive == 0) return 0;
        return static_cast<uint32_t>(Next() % max_exclusive);
    }

private:
    uint64_t state_;
};

}  // namespace

struct ClusterRepSelector::InternalState {
    FastRng rng{42};
    std::vector<float> distance_sums;
    std::unordered_set<uint32_t> rep_set;
};

ClusterRepSelector::ClusterRepSelector()
    : state_(std::make_unique<InternalState>()) {}

ClusterRepSelector::ClusterRepSelector(const ClusterRepConfig& config)
    : config_(config), state_(std::make_unique<InternalState>()) {}

ClusterRepSelector::~ClusterRepSelector() = default;

ClusterRepSelector::ClusterRepSelector(ClusterRepSelector&&) noexcept = default;
ClusterRepSelector& ClusterRepSelector::operator=(ClusterRepSelector&&) noexcept = default;

std::unique_ptr<ClusterRepResult> ClusterRepSelector::Select(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto result = std::make_unique<ClusterRepResult>();
    const size_t n_reads = clustering.labels.size();

    if (n_reads == 0 || clustering.num_communities == 0) {
        return result;
    }

    auto community_members = clustering.GetCommunityMembers();

    result->representatives.reserve(clustering.num_communities);
    result->representative_reads.reserve(clustering.num_communities);
    state_->rep_set.clear();
    state_->rng = FastRng(config_.seed);

    float total_confidence = 0.0f;
    size_t total_cluster_size = 0;

    for (uint32_t cluster_id = 0; cluster_id < community_members.size(); ++cluster_id) {
        const auto& members = community_members[cluster_id];

        if (members.size() < config_.min_cluster_size) {
            ++result->stats.clusters_skipped;
            continue;
        }

        auto rep_opt = SelectForCluster(graph, members);
        if (!rep_opt.has_value()) {
            ++result->stats.clusters_skipped;
            continue;
        }

        uint32_t rep_read = *rep_opt;

        RepresentativeInfo info;
        info.cluster_id = cluster_id;
        info.read_idx = rep_read;
        info.confidence = 1.0f;
        info.cluster_size = members.size();
        info.avg_distance_to_members = 0.0f;
        info.centrality_score = ComputeCentrality(graph, members, rep_read);

        if (members.size() > 1) {
            float dist_sum = ComputeDistanceSum(graph, members, rep_read);
            info.avg_distance_to_members = dist_sum / static_cast<float>(members.size() - 1);
            result->stats.total_within_cluster_distance += dist_sum;
        }

        result->representatives.push_back(info);
        result->representative_reads.push_back(rep_read);
        state_->rep_set.insert(rep_read);

        total_confidence += info.confidence;
        total_cluster_size += members.size();
    }

    result->stats.num_clusters = clustering.num_communities;
    result->stats.num_representatives = result->representatives.size();

    if (!result->representatives.empty()) {
        result->stats.avg_rep_confidence =
            total_confidence / static_cast<float>(result->representatives.size());
        result->stats.avg_cluster_size =
            static_cast<float>(total_cluster_size) / static_cast<float>(result->representatives.size());
    }

    auto end = std::chrono::high_resolution_clock::now();
    result->stats.total_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return result;
}

std::unique_ptr<ClusterRepResult> ClusterRepSelector::SelectWithConfidence(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering,
    const SelfWaveCollapseResult& swc_result)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto result = std::make_unique<ClusterRepResult>();
    const size_t n_reads = clustering.labels.size();

    if (n_reads == 0 || clustering.num_communities == 0) {
        return result;
    }

    std::vector<float> confidence_by_read(n_reads, 0.0f);
    for (const auto& assignment : swc_result.assignments) {
        if (assignment.read_idx < n_reads) {
            confidence_by_read[assignment.read_idx] = assignment.confidence;
        }
    }

    auto community_members = clustering.GetCommunityMembers();

    result->representatives.reserve(clustering.num_communities);
    result->representative_reads.reserve(clustering.num_communities);
    state_->rep_set.clear();
    state_->rng = FastRng(config_.seed);

    float total_confidence = 0.0f;
    size_t total_cluster_size = 0;

    for (uint32_t cluster_id = 0; cluster_id < community_members.size(); ++cluster_id) {
        const auto& members = community_members[cluster_id];

        if (members.size() < config_.min_cluster_size) {
            ++result->stats.clusters_skipped;
            continue;
        }

        std::vector<float> member_confidences;
        member_confidences.reserve(members.size());
        for (uint32_t read_idx : members) {
            member_confidences.push_back(confidence_by_read[read_idx]);
        }

        auto rep_opt = SelectForClusterWithConfidence(
            graph, members, member_confidences);

        if (!rep_opt.has_value()) {
            ++result->stats.clusters_skipped;
            continue;
        }

        uint32_t rep_read = *rep_opt;

        RepresentativeInfo info;
        info.cluster_id = cluster_id;
        info.read_idx = rep_read;
        info.confidence = confidence_by_read[rep_read];
        info.cluster_size = members.size();
        info.avg_distance_to_members = 0.0f;
        info.centrality_score = ComputeCentrality(graph, members, rep_read);

        if (members.size() > 1) {
            float dist_sum = ComputeDistanceSum(graph, members, rep_read);
            info.avg_distance_to_members = dist_sum / static_cast<float>(members.size() - 1);
            result->stats.total_within_cluster_distance += dist_sum;
        }

        result->representatives.push_back(info);
        result->representative_reads.push_back(rep_read);
        state_->rep_set.insert(rep_read);

        total_confidence += info.confidence;
        total_cluster_size += members.size();
    }

    result->stats.num_clusters = clustering.num_communities;
    result->stats.num_representatives = result->representatives.size();

    if (!result->representatives.empty()) {
        result->stats.avg_rep_confidence =
            total_confidence / static_cast<float>(result->representatives.size());
        result->stats.avg_cluster_size =
            static_cast<float>(total_cluster_size) / static_cast<float>(result->representatives.size());
    }

    auto end = std::chrono::high_resolution_clock::now();
    result->stats.total_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return result;
}

}  // namespace llmap::self_interference
