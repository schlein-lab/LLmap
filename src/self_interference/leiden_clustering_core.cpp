#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <chrono>
#include <numeric>
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

}  // namespace llmap::self_interference
