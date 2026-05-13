#include "self_interference/self_wavecollapse.h"
#include "self_interference/self_wavecollapse_impl.h"

#include <algorithm>
#include <chrono>
#include <unordered_set>

namespace llmap::self_interference {

SelfWaveCollapse::SelfWaveCollapse() = default;

SelfWaveCollapse::SelfWaveCollapse(const SelfWaveCollapseConfig& config)
    : config_(config) {}

SelfWaveCollapse::~SelfWaveCollapse() = default;

SelfWaveCollapse::SelfWaveCollapse(SelfWaveCollapse&&) noexcept = default;
SelfWaveCollapse& SelfWaveCollapse::operator=(SelfWaveCollapse&&) noexcept = default;

std::unique_ptr<SelfWaveCollapseResult> SelfWaveCollapse::Refine(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto result = std::make_unique<SelfWaveCollapseResult>();
    const size_t n_reads = clustering.labels.size();

    if (n_reads == 0) {
        return result;
    }

    // Initialize assignments
    result->assignments.resize(n_reads);
    for (uint32_t i = 0; i < n_reads; ++i) {
        result->assignments[i].read_idx = i;
        result->assignments[i].cluster_id = clustering.labels[i];
        result->assignments[i].confidence = 0.0f;
        result->assignments[i].collapsed = false;
        result->assignments[i].anchor_read = i;  // Self initially
    }

    // Group reads by cluster
    auto community_members = clustering.GetCommunityMembers();

    result->stats.num_reads = n_reads;
    result->stats.num_clusters = clustering.num_communities;

    float total_max_prob = 0.0f;
    size_t processed_reads = 0;

    // Process each cluster
    for (uint32_t cluster_id = 0; cluster_id < community_members.size(); ++cluster_id) {
        const auto& members = community_members[cluster_id];

        if (members.size() < config_.min_cluster_size) {
            // Small cluster: mark all as collapsed to self
            for (uint32_t read_idx : members) {
                result->assignments[read_idx].confidence = 1.0f;
                result->assignments[read_idx].collapsed = true;
                result->assignments[read_idx].anchor_read = read_idx;
            }
            result->stats.reads_collapsed += members.size();
            continue;
        }

        // Run EM on this cluster
        auto cluster_assignments = RefineCluster(graph, members, cluster_id);

        // Copy results
        for (const auto& assignment : cluster_assignments) {
            result->assignments[assignment.read_idx] = assignment;

            if (assignment.collapsed) {
                ++result->stats.reads_collapsed;
            }

            // Track quality metrics
            total_max_prob += assignment.confidence;
            ++processed_reads;
        }
    }

    result->num_clusters = clustering.num_communities;

    // Compute final stats
    auto end = std::chrono::high_resolution_clock::now();
    result->stats.total_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    if (processed_reads > 0) {
        result->stats.avg_max_probability = total_max_prob / static_cast<float>(processed_reads);
    }

    return result;
}

std::vector<ReadAssignment> SelfWaveCollapse::RefineCluster(
    const SimilarityGraph& graph,
    std::span<const uint32_t> cluster_members,
    uint32_t cluster_id)
{
    std::vector<ReadAssignment> assignments;
    const size_t n = cluster_members.size();

    if (n == 0) {
        return assignments;
    }

    // Single-read cluster: trivially collapsed
    if (n == 1) {
        ReadAssignment a;
        a.read_idx = cluster_members[0];
        a.cluster_id = cluster_id;
        a.confidence = 1.0f;
        a.collapsed = true;
        a.anchor_read = cluster_members[0];
        assignments.push_back(a);
        return assignments;
    }

    // Initialize state
    state_ = std::make_unique<InternalState>();
    state_->rng = FastRng(config_.seed + cluster_id);

    // Build local <-> global mappings
    state_->local_to_global.assign(cluster_members.begin(), cluster_members.end());
    for (size_t i = 0; i < n; ++i) {
        state_->global_to_local[cluster_members[i]] = static_cast<uint32_t>(i);
    }

    // Allocate probability matrices
    const size_t matrix_size = n * n;
    state_->prob_matrix.resize(matrix_size);
    state_->prob_matrix_new.resize(matrix_size);
    state_->similarity_cache.resize(matrix_size, -1.0f);  // -1 = not computed
    state_->collapsed.resize(n, false);
    state_->collapsed_to.resize(n, 0);

    // Initialize probability matrix: uniform distribution
    const float uniform_prob = 1.0f / static_cast<float>(n);
    std::fill(state_->prob_matrix.begin(), state_->prob_matrix.end(), uniform_prob);

    // Pre-compute similarity cache for cluster members
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            float sim = ComputeSimilarityWeight(
                graph,
                state_->local_to_global[i],
                state_->local_to_global[j]);

            state_->similarity_cache[i * n + j] = sim;
            state_->similarity_cache[j * n + i] = sim;
        }
    }

    // Run EM iterations
    RunEM(*state_, graph, cluster_members);

    // Extract results
    assignments.resize(n);

    // Find cluster anchor (most central read)
    uint32_t global_anchor = FindAnchor(graph, cluster_members);

    for (size_t i = 0; i < n; ++i) {
        ReadAssignment& a = assignments[i];
        a.read_idx = state_->local_to_global[i];
        a.cluster_id = cluster_id;

        // Find max probability and its anchor
        float max_prob = 0.0f;
        uint32_t max_anchor_local = 0;

        for (size_t j = 0; j < n; ++j) {
            float p = state_->prob_matrix[i * n + j];
            if (p > max_prob) {
                max_prob = p;
                max_anchor_local = static_cast<uint32_t>(j);
            }
        }

        a.confidence = max_prob;
        a.collapsed = (max_prob >= config_.collapse_threshold) || state_->collapsed[i];
        a.anchor_read = state_->local_to_global[max_anchor_local];

        // If collapsed, could also use the cluster anchor
        if (a.collapsed && max_anchor_local == i) {
            a.anchor_read = global_anchor;
        }
    }

    state_.reset();
    return assignments;
}

uint32_t SelfWaveCollapse::FindAnchor(
    const SimilarityGraph& graph,
    std::span<const uint32_t> members) const
{
    if (members.empty()) {
        return 0;
    }

    if (members.size() == 1) {
        return members[0];
    }

    // Find read with highest weighted degree within cluster
    std::unordered_set<uint32_t> member_set(members.begin(), members.end());

    uint32_t best_anchor = members[0];
    float best_centrality = 0.0f;

    for (uint32_t read_idx : members) {
        float centrality = 0.0f;

        auto neighbors = graph.Neighbors(read_idx);
        auto weights = graph.NeighborWeights(read_idx);

        for (size_t j = 0; j < neighbors.size(); ++j) {
            if (member_set.count(neighbors[j])) {
                centrality += weights[j];
            }
        }

        if (centrality > best_centrality) {
            best_centrality = centrality;
            best_anchor = read_idx;
        }
    }

    return best_anchor;
}

}  // namespace llmap::self_interference
