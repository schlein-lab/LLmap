#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace llmap::self_interference {

namespace {

// Fast RNG for sampling
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

// Internal state for ClusterRepSelector
struct ClusterRepSelector::InternalState {
    FastRng rng{42};

    // Cached distance sums for current cluster
    std::vector<float> distance_sums;

    // For tracking representative read indices quickly
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

    // Group reads by cluster
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
        info.confidence = 1.0f;  // No confidence info available
        info.cluster_size = members.size();
        info.avg_distance_to_members = 0.0f;
        info.centrality_score = ComputeCentrality(graph, members, rep_read);

        // Compute average distance to members
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

    // Compute stats
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

    // Build confidence lookup by read index
    std::vector<float> confidence_by_read(n_reads, 0.0f);
    for (const auto& assignment : swc_result.assignments) {
        if (assignment.read_idx < n_reads) {
            confidence_by_read[assignment.read_idx] = assignment.confidence;
        }
    }

    // Group reads by cluster
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

        // Extract confidence scores for cluster members
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

        // Compute average distance to members
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

    // Compute stats
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
            // Without confidence scores, fall back to max degree
            return ComputeMaxDegree(graph, cluster_members);

        case ClusterRepConfig::Method::Centroid:
            // Without embeddings, fall back to medoid
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
        // Fallback to regular selection
        return SelectForCluster(graph, cluster_members);
    }

    // For medoid with confidence: use confidence as a tiebreaker
    switch (config_.method) {
        case ClusterRepConfig::Method::MaxConfidence:
            return ComputeMaxConfidence(cluster_members, confidence_scores);

        case ClusterRepConfig::Method::Medoid: {
            // Compute medoid, but use confidence for ties
            uint32_t medoid = ComputeExactMedoid(graph, cluster_members);

            if (config_.use_confidence_tiebreaker) {
                // Find the medoid's position to get its confidence
                size_t medoid_idx = 0;
                for (size_t i = 0; i < cluster_members.size(); ++i) {
                    if (cluster_members[i] == medoid) {
                        medoid_idx = i;
                        break;
                    }
                }

                float medoid_dist = ComputeDistanceSum(graph, cluster_members, medoid);

                // Check if any other member has same distance but higher confidence
                for (size_t i = 0; i < cluster_members.size(); ++i) {
                    if (cluster_members[i] == medoid) continue;

                    float dist = ComputeDistanceSum(graph, cluster_members, cluster_members[i]);
                    // Allow small tolerance for distance equality
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

    // Determine candidate set
    std::span<const uint32_t> candidates = members;
    std::vector<uint32_t> limited_candidates;

    if (config_.max_medoid_candidates > 0 &&
        members.size() > config_.max_medoid_candidates) {
        // Sample candidates
        limited_candidates.reserve(config_.max_medoid_candidates);
        std::vector<uint32_t> all_members(members.begin(), members.end());

        for (size_t i = 0; i < config_.max_medoid_candidates; ++i) {
            size_t swap_idx = i + state_->rng.NextInRange(
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

    // Sample a subset of members
    std::vector<uint32_t> sample;
    sample.reserve(config_.approx_sample_size);

    if (members.size() <= config_.approx_sample_size) {
        sample.assign(members.begin(), members.end());
    } else {
        std::vector<uint32_t> all_members(members.begin(), members.end());
        for (size_t i = 0; i < config_.approx_sample_size; ++i) {
            size_t swap_idx = i + state_->rng.NextInRange(
                static_cast<uint32_t>(all_members.size() - i));
            std::swap(all_members[i], all_members[swap_idx]);
            sample.push_back(all_members[i]);
        }
    }

    // Find medoid within sample
    uint32_t best_medoid = sample[0];
    float best_dist_sum = std::numeric_limits<float>::max();

    for (uint32_t candidate : sample) {
        // Compute distance to all sample members (not all cluster members)
        float dist_sum = 0.0f;
        for (uint32_t other : sample) {
            if (other != candidate) {
                float weight = graph.GetEdgeWeight(candidate, other);
                // Distance = 1 - similarity
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
        // Distance = 1 - similarity (weight is already similarity)
        // If no edge exists, assume maximum distance (1.0)
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

// ClusterRepResult methods

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

// Convenience functions

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
