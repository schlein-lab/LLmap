#include "self_interference/cluster_rep.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

FastRng& GetRng(uint64_t seed) {
    thread_local FastRng rng(42);
    rng = FastRng(seed);
    return rng;
}

}  // namespace

std::optional<uint32_t> ClusterRepSelector::SelectForCluster(
    const SimilarityGraph& graph,
    std::span<const uint32_t> cluster_members)
{
    if (cluster_members.empty()) {
        return std::nullopt;
    }

    if (cluster_members.size() == 1) {
        return cluster_members[0];
    }

    switch (config_.method) {
        case ClusterRepConfig::Method::Medoid:
            if (config_.use_approximate_medoid &&
                cluster_members.size() > config_.approx_sample_size) {
                return ComputeApproxMedoid(graph, cluster_members);
            }
            return ComputeExactMedoid(graph, cluster_members);

        case ClusterRepConfig::Method::MaxDegree:
            return ComputeMaxDegree(graph, cluster_members);

        case ClusterRepConfig::Method::MaxConfidence:
            return ComputeMaxDegree(graph, cluster_members);

        case ClusterRepConfig::Method::Centroid:
            return ComputeExactMedoid(graph, cluster_members);
    }

    return ComputeExactMedoid(graph, cluster_members);
}

std::optional<uint32_t> ClusterRepSelector::SelectForClusterWithConfidence(
    const SimilarityGraph& graph,
    std::span<const uint32_t> cluster_members,
    std::span<const float> confidence_scores)
{
    if (cluster_members.empty()) {
        return std::nullopt;
    }

    if (cluster_members.size() == 1) {
        return cluster_members[0];
    }

    if (confidence_scores.size() != cluster_members.size()) {
        return SelectForCluster(graph, cluster_members);
    }

    switch (config_.method) {
        case ClusterRepConfig::Method::MaxConfidence:
            return ComputeMaxConfidence(cluster_members, confidence_scores);

        case ClusterRepConfig::Method::Medoid: {
            uint32_t medoid = ComputeExactMedoid(graph, cluster_members);

            if (config_.use_confidence_tiebreaker) {
                size_t medoid_idx = 0;
                for (size_t i = 0; i < cluster_members.size(); ++i) {
                    if (cluster_members[i] == medoid) {
                        medoid_idx = i;
                        break;
                    }
                }

                float medoid_dist = ComputeDistanceSum(graph, cluster_members, medoid);

                for (size_t i = 0; i < cluster_members.size(); ++i) {
                    if (cluster_members[i] == medoid) continue;

                    float dist = ComputeDistanceSum(graph, cluster_members, cluster_members[i]);
                    if (std::abs(dist - medoid_dist) < 1e-6f) {
                        if (confidence_scores[i] > confidence_scores[medoid_idx]) {
                            medoid = cluster_members[i];
                            medoid_idx = i;
                        }
                    }
                }
            }

            return medoid;
        }

        default:
            return SelectForCluster(graph, cluster_members);
    }
}

uint32_t ClusterRepSelector::ComputeExactMedoid(
    const SimilarityGraph& graph,
    std::span<const uint32_t> members)
{
    if (members.empty()) return 0;
    if (members.size() == 1) return members[0];

    uint32_t best_medoid = members[0];
    float best_dist_sum = std::numeric_limits<float>::max();

    std::span<const uint32_t> candidates = members;
    std::vector<uint32_t> limited_candidates;

    if (config_.max_medoid_candidates > 0 &&
        members.size() > config_.max_medoid_candidates) {
        limited_candidates.reserve(config_.max_medoid_candidates);
        std::vector<uint32_t> all_members(members.begin(), members.end());
        auto& rng = GetRng(config_.seed);

        for (size_t i = 0; i < config_.max_medoid_candidates; ++i) {
            size_t swap_idx = i + rng.NextInRange(
                static_cast<uint32_t>(all_members.size() - i));
            std::swap(all_members[i], all_members[swap_idx]);
            limited_candidates.push_back(all_members[i]);
        }
        candidates = limited_candidates;
    }

    for (uint32_t candidate : candidates) {
        float dist_sum = ComputeDistanceSum(graph, members, candidate);

        if (dist_sum < best_dist_sum) {
            best_dist_sum = dist_sum;
            best_medoid = candidate;
        }
    }

    return best_medoid;
}

uint32_t ClusterRepSelector::ComputeApproxMedoid(
    const SimilarityGraph& graph,
    std::span<const uint32_t> members)
{
    if (members.empty()) return 0;
    if (members.size() == 1) return members[0];

    std::vector<uint32_t> sample;
    sample.reserve(config_.approx_sample_size);

    if (members.size() <= config_.approx_sample_size) {
        sample.assign(members.begin(), members.end());
    } else {
        std::vector<uint32_t> all_members(members.begin(), members.end());
        auto& rng = GetRng(config_.seed);
        for (size_t i = 0; i < config_.approx_sample_size; ++i) {
            size_t swap_idx = i + rng.NextInRange(
                static_cast<uint32_t>(all_members.size() - i));
            std::swap(all_members[i], all_members[swap_idx]);
            sample.push_back(all_members[i]);
        }
    }

    uint32_t best_medoid = sample[0];
    float best_dist_sum = std::numeric_limits<float>::max();

    for (uint32_t candidate : sample) {
        float dist_sum = 0.0f;
        for (uint32_t other : sample) {
            if (other != candidate) {
                float weight = graph.GetEdgeWeight(candidate, other);
                dist_sum += (weight > 0.0f) ? (1.0f - weight) : 1.0f;
            }
        }

        if (dist_sum < best_dist_sum) {
            best_dist_sum = dist_sum;
            best_medoid = candidate;
        }
    }

    return best_medoid;
}

uint32_t ClusterRepSelector::ComputeMaxDegree(
    const SimilarityGraph& graph,
    std::span<const uint32_t> members)
{
    if (members.empty()) return 0;
    if (members.size() == 1) return members[0];

    std::unordered_set<uint32_t> member_set(members.begin(), members.end());

    uint32_t best_member = members[0];
    float best_degree = 0.0f;

    for (uint32_t read_idx : members) {
        float weighted_degree = 0.0f;

        auto neighbors = graph.Neighbors(read_idx);
        auto weights = graph.NeighborWeights(read_idx);

        for (size_t j = 0; j < neighbors.size(); ++j) {
            if (member_set.count(neighbors[j])) {
                weighted_degree += weights[j];
            }
        }

        if (weighted_degree > best_degree) {
            best_degree = weighted_degree;
            best_member = read_idx;
        }
    }

    return best_member;
}

uint32_t ClusterRepSelector::ComputeMaxConfidence(
    std::span<const uint32_t> members,
    std::span<const float> confidence_scores)
{
    if (members.empty()) return 0;
    if (members.size() == 1) return members[0];

    uint32_t best_member = members[0];
    float best_confidence = confidence_scores[0];

    for (size_t i = 1; i < members.size(); ++i) {
        if (confidence_scores[i] > best_confidence) {
            best_confidence = confidence_scores[i];
            best_member = members[i];
        }
    }

    return best_member;
}

float ClusterRepSelector::ComputeDistanceSum(
    const SimilarityGraph& graph,
    std::span<const uint32_t> members,
    uint32_t candidate)
{
    float dist_sum = 0.0f;

    for (uint32_t other : members) {
        if (other == candidate) continue;

        float weight = graph.GetEdgeWeight(candidate, other);
        dist_sum += (weight > 0.0f) ? (1.0f - weight) : 1.0f;
    }

    return dist_sum;
}

float ClusterRepSelector::ComputeCentrality(
    const SimilarityGraph& graph,
    std::span<const uint32_t> members,
    uint32_t candidate)
{
    std::unordered_set<uint32_t> member_set(members.begin(), members.end());
    float centrality = 0.0f;

    auto neighbors = graph.Neighbors(candidate);
    auto weights = graph.NeighborWeights(candidate);

    for (size_t j = 0; j < neighbors.size(); ++j) {
        if (member_set.count(neighbors[j])) {
            centrality += weights[j];
        }
    }

    return centrality;
}

}  // namespace llmap::self_interference
