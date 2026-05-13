#include "self_interference/leiden_clustering.h"
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

// Fast linear congruential generator for shuffling
class FastRng {
public:
    explicit FastRng(uint64_t seed) : state_(seed) {}

    uint64_t Next() {
        state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
        return state_;
    }

    size_t NextInRange(size_t max) {
        return Next() % max;
    }

private:
    uint64_t state_;
};

}  // namespace

// Internal state for the Leiden algorithm
struct LeidenClustering::InternalState {
    // Graph reference
    const SimilarityGraph* graph = nullptr;

    // Node-level state
    std::vector<uint32_t> labels;           // Current community assignment
    std::vector<float> node_weights;        // Sum of incident edge weights per node

    // Community-level state
    std::vector<float> community_weights;   // Sum of node weights in community
    std::vector<float> community_internal;  // Sum of internal edge weights
    size_t num_communities = 0;

    // Total graph weight (sum of all edge weights)
    double total_weight = 0.0;

    // Resolution parameter
    float resolution = 1.0f;

    // Random number generator
    FastRng rng{42};

    // Iteration tracking
    size_t iteration = 0;
    float current_modularity = 0.0f;

    // Timing
    float local_move_time_ms = 0.0f;
    float refinement_time_ms = 0.0f;
    float aggregation_time_ms = 0.0f;
};

LeidenClustering::LeidenClustering() = default;

LeidenClustering::LeidenClustering(const LeidenConfig& config)
    : config_(config) {}

LeidenClustering::~LeidenClustering() = default;

LeidenClustering::LeidenClustering(LeidenClustering&&) noexcept = default;
LeidenClustering& LeidenClustering::operator=(LeidenClustering&&) noexcept = default;

std::unique_ptr<ClusteringResult> LeidenClustering::Cluster(
    const SimilarityGraph& graph)
{
    // Initialize each node in its own community
    std::vector<uint32_t> initial_labels(graph.NumNodes());
    std::iota(initial_labels.begin(), initial_labels.end(), 0u);

    return ClusterWithInitialLabels(graph, initial_labels);
}

std::unique_ptr<ClusteringResult> LeidenClustering::ClusterWithInitialLabels(
    const SimilarityGraph& graph,
    std::span<const uint32_t> initial_labels)
{
    auto start = std::chrono::high_resolution_clock::now();

    const size_t n = graph.NumNodes();
    if (n == 0) {
        auto result = std::make_unique<ClusteringResult>();
        result->num_communities = 0;
        result->modularity = 0.0f;
        return result;
    }

    if (initial_labels.size() != n) {
        return nullptr;
    }

    // Initialize internal state
    state_ = std::make_unique<InternalState>();
    state_->graph = &graph;
    state_->rng = FastRng(config_.seed);
    state_->resolution = config_.resolution;

    // Copy initial labels
    state_->labels.assign(initial_labels.begin(), initial_labels.end());

    // Compute node weights (degree in weighted sense)
    state_->node_weights.resize(n, 0.0f);
    state_->total_weight = 0.0;

    for (size_t i = 0; i < n; ++i) {
        auto weights = graph.NeighborWeights(i);
        float sum = 0.0f;
        for (float w : weights) {
            sum += w;
        }
        state_->node_weights[i] = sum;
        state_->total_weight += sum;
    }
    state_->total_weight /= 2.0;  // Each edge counted twice

    // Initialize community weights
    uint32_t max_label = *std::max_element(state_->labels.begin(), state_->labels.end());
    state_->community_weights.resize(max_label + 1, 0.0f);
    state_->community_internal.resize(max_label + 1, 0.0f);

    for (size_t i = 0; i < n; ++i) {
        state_->community_weights[state_->labels[i]] += state_->node_weights[i];
    }

    // Compute initial internal edge weights
    for (size_t i = 0; i < n; ++i) {
        auto neighbors = graph.Neighbors(i);
        auto weights = graph.NeighborWeights(i);
        uint32_t my_comm = state_->labels[i];

        for (size_t j = 0; j < neighbors.size(); ++j) {
            if (state_->labels[neighbors[j]] == my_comm && neighbors[j] > i) {
                state_->community_internal[my_comm] += weights[j];
            }
        }
    }

    // Count initial communities
    std::unordered_set<uint32_t> active_communities(
        state_->labels.begin(), state_->labels.end());
    state_->num_communities = active_communities.size();

    // Main iteration loop
    bool improved = true;
    state_->iteration = 0;

    while (improved && state_->iteration < config_.max_iterations) {
        ++state_->iteration;

        // Phase 1: Local move
        improved = LocalMovePhase(*state_);

        // Phase 2: Refinement
        if (config_.enable_refinement) {
            RefinementPhase(*state_);
        }

        // Note: Full aggregation phase is skipped for CPU implementation
        // as it requires graph reconstruction. The local move already
        // achieves good quality for most use cases.
    }

    // Compute final modularity
    state_->current_modularity = CalculateModularity(
        graph, state_->labels, config_.resolution);

    // Build result
    auto result = std::make_unique<ClusteringResult>();
    result->labels = std::move(state_->labels);
    result->modularity = state_->current_modularity;

    // Relabel contiguously
    result->RelabelContiguous();
    result->num_communities = 0;
    if (!result->labels.empty()) {
        result->num_communities = *std::max_element(
            result->labels.begin(), result->labels.end()) + 1;
    }

    // Apply community size filter if configured
    if (config_.min_community_size > 1 || config_.max_community_size > 0) {
        *result = result->FilterByCommunitySize(
            config_.min_community_size,
            config_.max_community_size);
    }

    auto end = std::chrono::high_resolution_clock::now();

    // Fill stats
    result->stats.num_nodes = n;
    result->stats.num_communities = result->num_communities;
    result->stats.num_iterations = state_->iteration;
    result->stats.final_modularity = result->modularity;
    result->stats.time_ms = std::chrono::duration<float, std::milli>(end - start).count();
    result->stats.local_move_time_ms = state_->local_move_time_ms;
    result->stats.refinement_time_ms = state_->refinement_time_ms;

    // Community size stats
    auto sizes = result->GetCommunitySizes();
    if (!sizes.empty()) {
        result->stats.avg_community_size =
            static_cast<float>(n) / static_cast<float>(result->num_communities);
        result->stats.min_community_size = *std::min_element(sizes.begin(), sizes.end());
        result->stats.max_community_size = *std::max_element(sizes.begin(), sizes.end());
        result->stats.singleton_count = std::count(sizes.begin(), sizes.end(), 1UL);
    }

    state_.reset();
    return result;
}

bool LeidenClustering::LocalMovePhase(InternalState& state) {
    auto start = std::chrono::high_resolution_clock::now();

    const size_t n = state.graph->NumNodes();
    bool any_moved = false;

    // Create shuffled node order
    std::vector<size_t> node_order(n);
    std::iota(node_order.begin(), node_order.end(), 0);

    if (config_.randomize_node_order) {
        for (size_t i = n - 1; i > 0; --i) {
            size_t j = state.rng.NextInRange(i + 1);
            std::swap(node_order[i], node_order[j]);
        }
    }

    // Local move: try to move each node to a neighboring community
    bool improved_this_pass = true;

    while (improved_this_pass) {
        improved_this_pass = false;

        for (size_t idx : node_order) {
            const auto i = static_cast<uint32_t>(idx);
            const uint32_t current_comm = state.labels[i];
            auto neighbors = state.graph->Neighbors(i);
            auto weights = state.graph->NeighborWeights(i);

            if (neighbors.empty()) continue;

            // Find neighboring communities and their connection weights
            std::unordered_map<uint32_t, float> neighbor_comm_weights;
            for (size_t j = 0; j < neighbors.size(); ++j) {
                uint32_t neighbor_comm = state.labels[neighbors[j]];
                neighbor_comm_weights[neighbor_comm] += weights[j];
            }

            // Find best community to move to
            uint32_t best_comm = current_comm;
            float best_gain = 0.0f;

            for (const auto& [comm, edge_weight_to_comm] : neighbor_comm_weights) {
                if (comm == current_comm) continue;

                // Compute modularity gain for moving node i to community comm
                float gain = ComputeModularityGain(state, i, comm);

                if (gain > best_gain + config_.min_modularity_gain) {
                    best_gain = gain;
                    best_comm = comm;
                }
            }

            // Move node if beneficial
            if (best_comm != current_comm) {
                // Update community weights
                state.community_weights[current_comm] -= state.node_weights[i];
                state.community_weights[best_comm] += state.node_weights[i];

                // Update internal weights
                float edges_to_old = neighbor_comm_weights.count(current_comm)
                    ? neighbor_comm_weights[current_comm] : 0.0f;
                float edges_to_new = neighbor_comm_weights[best_comm];

                state.community_internal[current_comm] -= edges_to_old;
                state.community_internal[best_comm] += edges_to_new;

                // Update label
                state.labels[i] = best_comm;

                improved_this_pass = true;
                any_moved = true;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    state.local_move_time_ms += std::chrono::duration<float, std::milli>(end - start).count();

    return any_moved;
}

void LeidenClustering::RefinementPhase(InternalState& state) {
    auto start = std::chrono::high_resolution_clock::now();

    const size_t n = state.graph->NumNodes();

    // Refinement: within each community, check if nodes should form sub-communities
    // This is a simplified version - full Leiden does recursive partitioning

    std::unordered_map<uint32_t, std::vector<uint32_t>> community_members;
    for (uint32_t i = 0; i < n; ++i) {
        community_members[state.labels[i]].push_back(i);
    }

    // For each community with multiple members, check connectivity
    for (auto& [comm_id, members] : community_members) {
        if (members.size() < 2) continue;

        // Simple refinement: if a node has no internal edges, consider singleton
        for (uint32_t node : members) {
            auto neighbors = state.graph->Neighbors(node);
            bool has_internal = false;

            for (uint32_t neighbor : neighbors) {
                if (state.labels[neighbor] == comm_id && neighbor != node) {
                    has_internal = true;
                    break;
                }
            }

            // Node has no internal edges - might be poorly connected
            // (In full Leiden, this triggers sub-partition; here we just note it)
            (void)has_internal;  // Placeholder for more sophisticated refinement
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    state.refinement_time_ms += std::chrono::duration<float, std::milli>(end - start).count();
}

void LeidenClustering::AggregationPhase(InternalState& state) {
    // Aggregation: build a coarser graph where each community becomes a node
    // This is computationally expensive and skipped in the basic CPU implementation.
    // The local move phase achieves good results without full aggregation.
    (void)state;
}

float LeidenClustering::ComputeModularityGain(
    const InternalState& state,
    uint32_t node,
    uint32_t target_community) const
{
    // Modularity gain formula for moving node from current community to target
    // ΔQ = [e_in / m - γ(Σ_tot + k_i)² / (2m)²] - [e_in' / m - γ(Σ_tot')² / (2m)²]
    //
    // Simplified: ΔQ ≈ (k_i,in - k_i,out) / m - γ * k_i * (Σ_target - Σ_current) / (2m²)

    const double m = state.total_weight;
    if (m == 0) return 0.0f;

    const uint32_t current_comm = state.labels[node];
    const float k_i = state.node_weights[node];

    // Compute edges to current and target communities
    auto neighbors = state.graph->Neighbors(node);
    auto weights = state.graph->NeighborWeights(node);

    float k_i_current = 0.0f;
    float k_i_target = 0.0f;

    for (size_t j = 0; j < neighbors.size(); ++j) {
        uint32_t neighbor_comm = state.labels[neighbors[j]];
        if (neighbor_comm == current_comm) {
            k_i_current += weights[j];
        } else if (neighbor_comm == target_community) {
            k_i_target += weights[j];
        }
    }

    // Community totals (excluding node i from current)
    float sigma_current = state.community_weights[current_comm] - k_i;
    float sigma_target = state.community_weights[target_community];

    // Modularity gain
    double gain = (k_i_target - k_i_current) / m;
    gain -= state.resolution * k_i * (sigma_target - sigma_current) / (2.0 * m * m);

    return static_cast<float>(gain);
}

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

    // Q = Σ_c [e_c / m - γ(Σ_c / 2m)²]
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
