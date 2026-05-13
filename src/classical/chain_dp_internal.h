#pragma once

#include "classical/chain.h"

#include <span>
#include <vector>

namespace llmap::classical {
namespace internal {

// DP state for a single anchor (legacy non-scratch version)
struct DpState {
    int32_t score = 0;
    int32_t pred = -1;
};

// Run colinear DP for a single ref_id and strand (legacy version)
void ChainSingleRefStrand(
    std::span<const Anchor> anchors,
    std::span<const size_t> sorted_indices,
    std::vector<DpState>& dp,
    const ChainConfig& config);

// Backtrack from a chain endpoint (legacy version)
Chain BacktrackChain(
    std::span<const Anchor> anchors,
    const std::vector<DpState>& dp,
    size_t end_idx);

// Run colinear DP for a single ref_id and strand (scratch version)
void ChainSingleRefStrandScratch(
    std::span<const Anchor> anchors,
    std::span<const size_t> sorted_indices,
    core::ScratchBuffer<int32_t>& dp_score,
    core::ScratchBuffer<int32_t>& dp_pred,
    const ChainConfig& config);

// Backtrack using scratch buffer arrays
Chain BacktrackChainScratch(
    std::span<const Anchor> anchors,
    const core::ScratchBuffer<int32_t>& dp_score,
    const core::ScratchBuffer<int32_t>& dp_pred,
    size_t end_idx,
    core::ScratchBuffer<uint32_t>& backtrack);

}  // namespace internal
}  // namespace llmap::classical
