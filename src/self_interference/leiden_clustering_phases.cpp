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

// Fast linear congruential generator for shuffling (duplicated for linkage)
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

// Internal state mirror for phase methods
struct LeidenClustering::InternalState {
    const SimilarityGraph* graph = nullptr;
    std::vector<uint32_t> labels;
    std::vector<float> node_weights;
    std::vector<float> community_weights;
    std::vector<float> community_internal;
    size_t num_communities = 0;
    double total_weight = 0.0;
    float resolution = 1.0f;
    FastRng rng{42};
    size_t iteration = 0;
    float current_modularity = 0.0f;
    float local_move_time_ms = 0.0f;
    float refinement_time_ms = 0.0f;
    float aggregation_time_ms = 0.0f;
};

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
    // Delta_Q = [e_in / m - gamma(Sigma_tot + k_i)^2 / (2m)^2] - [e_in' / m - gamma(Sigma_tot')^2 / (2m)^2]
    //
    // Simplified: Delta_Q approx (k_i,in - k_i,out) / m - gamma * k_i * (Sigma_target - Sigma_current) / (2m^2)

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

}  // namespace llmap::self_interference
