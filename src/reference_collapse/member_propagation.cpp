#include "reference_collapse/member_propagation.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <unordered_set>

namespace llmap {

struct MemberPropagation::InternalState {
    std::vector<BucketProb> temp_buckets;
};

MemberPropagation::MemberPropagation()
    : state_(std::make_unique<InternalState>()) {}

MemberPropagation::MemberPropagation(const MemberPropagationConfig& config)
    : config_(config), state_(std::make_unique<InternalState>()) {}

MemberPropagation::~MemberPropagation() = default;

MemberPropagation::MemberPropagation(MemberPropagation&&) noexcept = default;
MemberPropagation& MemberPropagation::operator=(MemberPropagation&&) noexcept = default;

MemberPropagationResult MemberPropagation::Propagate(
    WaveState& state,
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph)
{
    auto start = std::chrono::high_resolution_clock::now();

    MemberPropagationResult result;
    result.stats.num_clusters = static_cast<std::uint32_t>(clustering.num_communities);
    result.stats.num_representatives = static_cast<std::uint32_t>(rep_result.representatives.size());

    if (rep_result.representatives.empty()) {
        return result;
    }

    auto community_members = clustering.GetCommunityMembers();

    float total_similarity = 0.0f;
    float total_cohesion = 0.0f;
    float total_confidence = 0.0f;
    std::uint32_t cohesion_count = 0;

    for (const auto& rep_info : rep_result.representatives) {
        if (rep_info.cluster_id >= community_members.size()) {
            continue;
        }

        const auto& members = community_members[rep_info.cluster_id];
        if (members.empty()) {
            continue;
        }

        auto entries = PropagateCluster(
            state,
            rep_info.read_idx,
            members,
            graph,
            rep_info.cluster_id);

        for (auto& entry : entries) {
            total_similarity += entry.similarity_to_rep;
            total_confidence += entry.propagated_confidence;
            result.entries.push_back(std::move(entry));
        }

        float cohesion = ComputeClusterCohesion(graph, members);
        total_cohesion += cohesion;
        ++cohesion_count;
    }

    result.stats.num_members_propagated = static_cast<std::uint32_t>(result.entries.size());

    if (!result.entries.empty()) {
        result.stats.avg_similarity_to_rep =
            total_similarity / static_cast<float>(result.entries.size());
        result.stats.avg_member_confidence =
            total_confidence / static_cast<float>(result.entries.size());
    }

    if (cohesion_count > 0) {
        result.stats.avg_cluster_cohesion = total_cohesion / static_cast<float>(cohesion_count);
    }

    for (const auto& entry : result.entries) {
        result.stats.total_buckets_propagated += entry.buckets_received;
    }

    auto end = std::chrono::high_resolution_clock::now();
    result.stats.propagation_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return result;
}

std::vector<MemberPropagationEntry> MemberPropagation::PropagateCluster(
    WaveState& state,
    std::uint32_t representative_idx,
    std::span<const std::uint32_t> member_indices,
    const self_interference::SimilarityGraph& graph,
    std::uint32_t cluster_id)
{
    std::vector<MemberPropagationEntry> entries;

    if (member_indices.empty()) {
        return entries;
    }

    if (representative_idx >= state.n_reads()) {
        return entries;
    }

    if (state.is_collapsed(representative_idx)) {
        return entries;
    }

    auto rep_buckets = state.bucket_indices_for_read(representative_idx);
    auto rep_probs = state.probabilities_for_read(representative_idx);

    if (rep_buckets.empty()) {
        return entries;
    }

    float cluster_cohesion = ComputeClusterCohesion(graph, member_indices);

    entries.reserve(member_indices.size() - 1);

    for (std::uint32_t member_idx : member_indices) {
        if (member_idx == representative_idx) {
            continue;
        }

        if (member_idx >= state.n_reads()) {
            continue;
        }

        if (state.is_collapsed(member_idx)) {
            continue;
        }

        float similarity = graph.GetEdgeWeight(representative_idx, member_idx);
        if (similarity <= 0.0f) {
            similarity = 0.01f;
        }

        std::uint32_t buckets_propagated = PropagateMember(
            state,
            representative_idx,
            member_idx,
            similarity,
            cluster_cohesion);

        MemberPropagationEntry entry;
        entry.member_idx = member_idx;
        entry.cluster_id = cluster_id;
        entry.representative_idx = representative_idx;
        entry.similarity_to_rep = similarity;
        entry.buckets_received = buckets_propagated;

        float rep_max_prob = 0.0f;
        for (float p : rep_probs) {
            rep_max_prob = std::max(rep_max_prob, p);
        }
        entry.propagated_confidence = ComputeMemberConfidence(
            rep_max_prob, similarity, cluster_cohesion);

        entries.push_back(entry);
    }

    return entries;
}

std::uint32_t MemberPropagation::PropagateMember(
    WaveState& state,
    std::uint32_t representative_idx,
    std::uint32_t member_idx,
    float similarity_to_rep,
    float cluster_cohesion)
{
    auto rep_buckets = state.bucket_indices_for_read(representative_idx);
    auto rep_probs = state.probabilities_for_read(representative_idx);

    if (rep_buckets.empty()) {
        return 0;
    }

    float rep_max_prob = 0.0f;
    for (float p : rep_probs) {
        rep_max_prob = std::max(rep_max_prob, p);
    }

    float member_confidence = ComputeMemberConfidence(
        rep_max_prob, similarity_to_rep, cluster_cohesion);

    state_->temp_buckets.clear();
    state_->temp_buckets.reserve(rep_buckets.size());

    for (std::size_t i = 0; i < rep_buckets.size(); ++i) {
        BucketProb bp;
        bp.bucket_id = rep_buckets[i];
        bp.probability = rep_probs[i] * member_confidence;
        state_->temp_buckets.push_back(bp);
    }

    std::sort(state_->temp_buckets.begin(), state_->temp_buckets.end(),
              [](const BucketProb& a, const BucketProb& b) {
                  return a.probability > b.probability;
              });

    if (config_.propagate_top_candidates_only &&
        state_->temp_buckets.size() > config_.max_propagated_candidates) {
        state_->temp_buckets.resize(config_.max_propagated_candidates);
    }

    std::sort(state_->temp_buckets.begin(), state_->temp_buckets.end(),
              [](const BucketProb& a, const BucketProb& b) {
                  return a.bucket_id < b.bucket_id;
              });

    float prob_sum = 0.0f;
    for (const auto& bp : state_->temp_buckets) {
        prob_sum += bp.probability;
    }
    if (prob_sum > 0.0f) {
        for (auto& bp : state_->temp_buckets) {
            bp.probability /= prob_sum;
        }
    }

    state.set_read_candidates(member_idx, state_->temp_buckets);
    state.set_level(member_idx, state.get_level(representative_idx));

    return static_cast<std::uint32_t>(state_->temp_buckets.size());
}

float MemberPropagation::ComputeClusterCohesion(
    const self_interference::SimilarityGraph& graph,
    std::span<const std::uint32_t> member_indices) const
{
    if (member_indices.size() <= 1) {
        return 1.0f;
    }

    float total_similarity = 0.0f;
    std::size_t pair_count = 0;

    constexpr std::size_t MAX_PAIRS_TO_SAMPLE = 100;
    std::size_t total_pairs = member_indices.size() * (member_indices.size() - 1) / 2;

    if (total_pairs <= MAX_PAIRS_TO_SAMPLE) {
        for (std::size_t i = 0; i < member_indices.size(); ++i) {
            for (std::size_t j = i + 1; j < member_indices.size(); ++j) {
                float sim = graph.GetEdgeWeight(member_indices[i], member_indices[j]);
                total_similarity += sim;
                ++pair_count;
            }
        }
    } else {
        std::size_t step = std::max<std::size_t>(1, member_indices.size() / 10);
        for (std::size_t i = 0; i < member_indices.size(); i += step) {
            for (std::size_t j = i + 1; j < member_indices.size(); j += step) {
                float sim = graph.GetEdgeWeight(member_indices[i], member_indices[j]);
                total_similarity += sim;
                ++pair_count;
            }
        }
    }

    return (pair_count > 0) ? (total_similarity / static_cast<float>(pair_count)) : 0.0f;
}

float MemberPropagation::ComputeMemberConfidence(
    float rep_confidence,
    float similarity_to_rep,
    float cluster_cohesion) const
{
    float base = rep_confidence * config_.base_confidence_scaling;

    float sim_factor = 1.0f - config_.similarity_weight + config_.similarity_weight * similarity_to_rep;

    float cohesion_factor = 1.0f - config_.cohesion_weight + config_.cohesion_weight * cluster_cohesion;

    float confidence = base * sim_factor * cohesion_factor;

    if (similarity_to_rep < 0.3f) {
        confidence *= config_.outlier_penalty;
    }

    confidence = std::clamp(confidence, config_.min_confidence, config_.max_confidence);

    return confidence;
}

void MemberPropagation::ScaleProbabilities(
    std::vector<BucketProb>& buckets,
    float confidence_factor) const
{
    for (auto& bp : buckets) {
        bp.probability *= confidence_factor;
    }
}

std::optional<MemberPropagationEntry> MemberPropagationResult::GetEntry(
    std::uint32_t member_idx) const
{
    for (const auto& entry : entries) {
        if (entry.member_idx == member_idx) {
            return entry;
        }
    }
    return std::nullopt;
}

std::vector<std::uint32_t> MemberPropagationResult::GetClusterMembers(
    std::uint32_t cluster_id) const
{
    std::vector<std::uint32_t> members;
    for (const auto& entry : entries) {
        if (entry.cluster_id == cluster_id) {
            members.push_back(entry.member_idx);
        }
    }
    return members;
}

MemberPropagationResult PropagateToMembers(
    WaveState& state,
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph)
{
    MemberPropagation propagation;
    return propagation.Propagate(state, clustering, rep_result, graph);
}

}  // namespace llmap
