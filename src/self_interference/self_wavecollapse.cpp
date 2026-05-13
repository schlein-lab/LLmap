#include "self_interference/self_wavecollapse.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace llmap::self_interference {

namespace {

// Fast RNG for tie-breaking
class FastRng {
public:
    explicit FastRng(uint64_t seed) : state_(seed) {}

    uint64_t Next() {
        state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
        return state_;
    }

    float NextFloat() {
        return static_cast<float>(Next() & 0xFFFFFF) / static_cast<float>(0x1000000);
    }

private:
    uint64_t state_;
};

}  // namespace

// Internal state for Self-WaveCollapse
struct SelfWaveCollapse::InternalState {
    // Working arrays for a single cluster
    std::vector<float> prob_matrix;       // [n_members * n_members] probabilities
    std::vector<float> prob_matrix_new;   // Double buffer for updates
    std::vector<float> similarity_cache;  // [n_members * n_members] cached similarities

    // Per-read state
    std::vector<bool> collapsed;          // Has this read collapsed?
    std::vector<uint32_t> collapsed_to;   // Which anchor did it collapse to?

    // Mapping from cluster-local index to global read index
    std::vector<uint32_t> local_to_global;
    std::unordered_map<uint32_t, uint32_t> global_to_local;

    // Statistics tracking
    size_t iteration_count = 0;
    bool converged = false;

    FastRng rng{42};
};

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
    // P(read i, anchor j) = 1/n initially
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

void SelfWaveCollapse::RunEM(
    InternalState& state,
    const SimilarityGraph& graph,
    std::span<const uint32_t> members)
{
    const size_t n = members.size();
    if (n == 0) return;

    state.iteration_count = 0;
    state.converged = false;

    for (size_t iter = 0; iter < config_.max_iterations; ++iter) {
        ++state.iteration_count;

        // E-step + M-step combined:
        // P_new(i, j) = (1-gamma) * P_old(i, j) + gamma * normalized_update(i, j)
        //
        // update(i, j) = S(i, j) * sum_k[P(k, j) * S(i, k)]
        // where S(i, j) = similarity_cache[i * n + j]

        for (size_t i = 0; i < n; ++i) {
            if (state.collapsed[i]) {
                // Already collapsed: keep distribution frozen
                for (size_t j = 0; j < n; ++j) {
                    state.prob_matrix_new[i * n + j] = state.prob_matrix[i * n + j];
                }
                continue;
            }

            // Compute update for read i
            float sum_update = 0.0f;
            std::vector<float> updates(n, 0.0f);

            for (size_t j = 0; j < n; ++j) {
                float s_ij = state.similarity_cache[i * n + j];

                // Sum over all reads k: P(k, j) * S(i, k)
                float support = 0.0f;
                for (size_t k = 0; k < n; ++k) {
                    float p_kj = state.prob_matrix[k * n + j];
                    float s_ik = state.similarity_cache[i * n + k];
                    support += p_kj * s_ik;
                }

                updates[j] = std::pow(s_ij, config_.similarity_exponent) * support;
                sum_update += updates[j];
            }

            // Normalize and apply damping
            if (sum_update > 0.0f) {
                for (size_t j = 0; j < n; ++j) {
                    float old_p = state.prob_matrix[i * n + j];
                    float new_p = updates[j] / sum_update;
                    state.prob_matrix_new[i * n + j] =
                        (1.0f - config_.gamma) * old_p + config_.gamma * new_p;
                }
            } else {
                // No update possible; keep old probabilities
                for (size_t j = 0; j < n; ++j) {
                    state.prob_matrix_new[i * n + j] = state.prob_matrix[i * n + j];
                }
            }

            // Renormalize to ensure sum = 1
            float row_sum = 0.0f;
            for (size_t j = 0; j < n; ++j) {
                row_sum += state.prob_matrix_new[i * n + j];
            }
            if (row_sum > 0.0f) {
                for (size_t j = 0; j < n; ++j) {
                    state.prob_matrix_new[i * n + j] /= row_sum;
                }
            }

            // Check for collapse
            float max_p = 0.0f;
            uint32_t max_j = 0;
            for (size_t j = 0; j < n; ++j) {
                if (state.prob_matrix_new[i * n + j] > max_p) {
                    max_p = state.prob_matrix_new[i * n + j];
                    max_j = static_cast<uint32_t>(j);
                }
            }

            if (max_p >= config_.collapse_threshold) {
                state.collapsed[i] = true;
                state.collapsed_to[i] = max_j;
            }
        }

        // Check convergence
        if (CheckConvergence(state.prob_matrix, state.prob_matrix_new)) {
            state.converged = true;
            std::swap(state.prob_matrix, state.prob_matrix_new);
            break;
        }

        // Swap buffers
        std::swap(state.prob_matrix, state.prob_matrix_new);
    }
}

float SelfWaveCollapse::ComputeSimilarityWeight(
    const SimilarityGraph& graph,
    uint32_t read_a,
    uint32_t read_b) const
{
    if (read_a == read_b) {
        return 1.0f;
    }

    // Get edge weight from similarity graph
    float edge_weight = graph.GetEdgeWeight(read_a, read_b);

    if (edge_weight > 0.0f) {
        return edge_weight;
    }

    // No direct edge: use 0 (no similarity contribution)
    return 0.0f;
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

bool SelfWaveCollapse::CheckConvergence(
    std::span<const float> old_probs,
    std::span<const float> new_probs) const
{
    if (old_probs.size() != new_probs.size()) {
        return false;
    }

    float max_diff = 0.0f;
    for (size_t i = 0; i < old_probs.size(); ++i) {
        float diff = std::abs(old_probs[i] - new_probs[i]);
        max_diff = std::max(max_diff, diff);
    }

    return max_diff < config_.convergence_threshold;
}

float SelfWaveCollapse::ComputeEntropy(std::span<const float> probs) {
    float entropy = 0.0f;
    for (float p : probs) {
        if (p > 0.0f) {
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

// SelfWaveCollapseResult methods

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

// Convenience function

std::unique_ptr<SelfWaveCollapseResult> RunSelfWaveCollapse(
    const SimilarityGraph& graph,
    const ClusteringResult& clustering)
{
    SelfWaveCollapse swc;
    return swc.Refine(graph, clustering);
}

}  // namespace llmap::self_interference
