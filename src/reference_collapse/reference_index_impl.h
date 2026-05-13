// LLmap — Reference Index PIMPL internal header.
// Shared by reference_index*.cpp files.

#pragma once

#include "reference_collapse/reference_index.h"

#include <unordered_map>
#include <vector>

namespace llmap {

// PIMPL implementation - shared across reference_index_*.cpp
class ReferenceIndexImpl {
public:
    ReferenceIndexConfig config;
    ReferenceIndexStats stats;
    BucketPyramid pyramid;
    std::vector<ReferenceTarget> targets;
    std::vector<float> l0_embeddings;
    std::vector<float> l1_embeddings;
    std::vector<float> l2_embeddings;

    // Spatial index for fast bucket lookup
    std::unordered_map<std::string, std::vector<std::uint32_t>> target_to_l2_buckets;

    void BuildSpatialIndex() {
        target_to_l2_buckets.clear();
        auto buckets = pyramid.l2_buckets();
        for (std::size_t i = 0; i < buckets.size(); ++i) {
            target_to_l2_buckets[buckets[i].target_id].push_back(
                static_cast<std::uint32_t>(i));
        }
    }
};

}  // namespace llmap
