#pragma once

// Internal implementation header for self_wavecollapse split files.
// Not for external use.

#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"
#include "self_interference/leiden_clustering.h"

#include <unordered_map>

namespace llmap::self_interference {

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

}  // namespace llmap::self_interference
