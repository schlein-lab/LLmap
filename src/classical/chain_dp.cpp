#include "classical/chain.h"
#include "classical/minimizer_index.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <numeric>

namespace llmap::classical {

namespace {

// Run colinear DP for a single ref_id and strand (uses separate score/pred arrays)
void ChainSingleRefStrandScratch(
    std::span<const Anchor> anchors,
    std::span<const size_t> sorted_indices,
    core::ScratchBuffer<int32_t>& dp_score,
    core::ScratchBuffer<int32_t>& dp_pred,
    const ChainConfig& config)
{
    const size_t n = sorted_indices.size();
    if (n == 0) return;

    // Initialize DP for this group
    for (size_t idx : sorted_indices) {
        dp_score[idx] = config.match_score;
        dp_pred[idx] = -1;
    }

    // Main DP loop
    for (size_t i = 1; i < n; ++i) {
        const size_t curr_idx = sorted_indices[i];
        const Anchor& curr = anchors[curr_idx];

        size_t skip_count = 0;
        for (size_t j = i; j > 0 && skip_count < config.max_skip; ) {
            --j;
            const size_t prev_idx = sorted_indices[j];
            const Anchor& prev = anchors[prev_idx];

            if (curr.ref_pos - prev.ref_pos > config.max_gap_ref) {
                break;
            }

            int32_t pair_score = AnchorPairScore(prev, curr, config);
            if (pair_score == std::numeric_limits<int32_t>::min()) {
                ++skip_count;
                continue;
            }

            int32_t new_score = dp_score[prev_idx] + pair_score;
            if (new_score > dp_score[curr_idx]) {
                dp_score[curr_idx] = new_score;
                dp_pred[curr_idx] = static_cast<int32_t>(prev_idx);
            }
        }
    }
}

// Backtrack using scratch buffer arrays
Chain BacktrackChainScratch(
    std::span<const Anchor> anchors,
    const core::ScratchBuffer<int32_t>& dp_score,
    const core::ScratchBuffer<int32_t>& dp_pred,
    size_t end_idx,
    core::ScratchBuffer<uint32_t>& backtrack)
{
    Chain chain;
    chain.score = dp_score[end_idx];

    backtrack.clear();
    int32_t idx = static_cast<int32_t>(end_idx);
    while (idx >= 0) {
        backtrack.push_back(static_cast<uint32_t>(idx));
        idx = dp_pred[idx];
    }

    // Reverse to start-to-end order
    chain.anchors.resize(backtrack.size());
    for (size_t i = 0; i < backtrack.size(); ++i) {
        chain.anchors[i] = backtrack[backtrack.size() - 1 - i];
    }

    const Anchor& first = anchors[chain.anchors.front()];
    const Anchor& last = anchors[chain.anchors.back()];

    chain.ref_id = first.ref_id;
    chain.is_forward = first.same_strand;
    chain.ref_start = first.ref_pos;
    chain.ref_end = last.ref_pos;

    if (first.same_strand) {
        chain.query_start = first.query_pos;
        chain.query_end = last.query_pos;
    } else {
        chain.query_start = last.query_pos;
        chain.query_end = first.query_pos;
    }

    return chain;
}

// DP state for a single anchor (legacy, used by non-scratch version)
struct DpState {
    int32_t score = 0;
    int32_t pred = -1;
};

// Run colinear DP for a single ref_id and strand
void ChainSingleRefStrand(
    std::span<const Anchor> anchors,
    std::span<const size_t> sorted_indices,
    std::vector<DpState>& dp,
    const ChainConfig& config)
{
    const size_t n = sorted_indices.size();
    if (n == 0) return;

    for (size_t idx : sorted_indices) {
        dp[idx].score = config.match_score;
        dp[idx].pred = -1;
    }

    for (size_t i = 1; i < n; ++i) {
        const size_t curr_idx = sorted_indices[i];
        const Anchor& curr = anchors[curr_idx];

        size_t skip_count = 0;
        for (size_t j = i; j > 0 && skip_count < config.max_skip; ) {
            --j;
            const size_t prev_idx = sorted_indices[j];
            const Anchor& prev = anchors[prev_idx];

            if (curr.ref_pos - prev.ref_pos > config.max_gap_ref) {
                break;
            }

            int32_t pair_score = AnchorPairScore(prev, curr, config);
            if (pair_score == std::numeric_limits<int32_t>::min()) {
                ++skip_count;
                continue;
            }

            int32_t new_score = dp[prev_idx].score + pair_score;
            if (new_score > dp[curr_idx].score) {
                dp[curr_idx].score = new_score;
                dp[curr_idx].pred = static_cast<int32_t>(prev_idx);
            }
        }
    }
}

// Backtrack from a chain endpoint to extract the full chain
Chain BacktrackChain(
    std::span<const Anchor> anchors,
    const std::vector<DpState>& dp,
    size_t end_idx)
{
    Chain chain;
    chain.score = dp[end_idx].score;

    // Collect anchor indices by backtracking
    std::vector<uint32_t> indices;
    int32_t idx = static_cast<int32_t>(end_idx);
    while (idx >= 0) {
        indices.push_back(static_cast<uint32_t>(idx));
        idx = dp[idx].pred;
    }

    // Reverse to get start-to-end order
    std::reverse(indices.begin(), indices.end());
    chain.anchors = std::move(indices);

    // Fill in chain metadata
    const Anchor& first = anchors[chain.anchors.front()];
    const Anchor& last = anchors[chain.anchors.back()];

    chain.ref_id = first.ref_id;
    chain.is_forward = first.same_strand;
    chain.ref_start = first.ref_pos;
    chain.ref_end = last.ref_pos;

    if (first.same_strand) {
        chain.query_start = first.query_pos;
        chain.query_end = last.query_pos;
    } else {
        chain.query_start = last.query_pos;
        chain.query_end = first.query_pos;
    }

    return chain;
}

}  // namespace

ChainResult ExtractChainsFromAnchors(
    std::span<const Anchor> anchors,
    uint32_t query_len,
    const ChainConfig& config)
{
    auto start = std::chrono::high_resolution_clock::now();

    ChainResult result;
    result.total_anchors = anchors.size();

    if (anchors.empty()) {
        return result;
    }

    // Group anchors by (ref_id, strand)
    // Key: (ref_id << 1) | (same_strand ? 1 : 0)
    std::vector<std::pair<uint64_t, size_t>> keyed_indices;
    keyed_indices.reserve(anchors.size());
    for (size_t i = 0; i < anchors.size(); ++i) {
        uint64_t key = (static_cast<uint64_t>(anchors[i].ref_id) << 1) |
                       (anchors[i].same_strand ? 1 : 0);
        keyed_indices.emplace_back(key, i);
    }

    // Sort by (key, ref_pos)
    std::sort(keyed_indices.begin(), keyed_indices.end(),
        [&](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first < b.first;
            return anchors[a.second].ref_pos < anchors[b.second].ref_pos;
        });

    // Allocate DP states
    std::vector<DpState> dp(anchors.size());
    std::vector<bool> used(anchors.size(), false);

    // Process each (ref_id, strand) group
    size_t group_start = 0;
    while (group_start < keyed_indices.size()) {
        uint64_t group_key = keyed_indices[group_start].first;
        size_t group_end = group_start + 1;
        while (group_end < keyed_indices.size() &&
               keyed_indices[group_end].first == group_key) {
            ++group_end;
        }

        // Extract indices for this group
        std::vector<size_t> group_indices;
        group_indices.reserve(group_end - group_start);
        for (size_t i = group_start; i < group_end; ++i) {
            group_indices.push_back(keyed_indices[i].second);
        }

        // Run DP
        ChainSingleRefStrand(anchors, group_indices, dp, config);

        group_start = group_end;
    }

    // Extract chains by finding high-scoring endpoints
    // Sort all anchors by DP score descending
    std::vector<size_t> score_order(anchors.size());
    std::iota(score_order.begin(), score_order.end(), 0);
    std::sort(score_order.begin(), score_order.end(),
        [&](size_t a, size_t b) { return dp[a].score > dp[b].score; });

    std::vector<Chain> all_chains;

    for (size_t end_idx : score_order) {
        if (used[end_idx]) continue;
        if (dp[end_idx].score < config.min_chain_score) break;

        // Backtrack to get chain
        Chain chain = BacktrackChain(anchors, dp, end_idx);

        // Skip if too few anchors
        if (chain.NumAnchors() < config.min_chain_anchors) {
            continue;
        }

        // Mark all anchors in this chain as used
        for (uint32_t idx : chain.anchors) {
            used[idx] = true;
        }

        all_chains.push_back(std::move(chain));
    }

    // Find best score
    if (!all_chains.empty()) {
        result.best_score = all_chains[0].score;  // Already sorted by score
    }

    // Filter by min_score_fraction
    int32_t score_threshold = static_cast<int32_t>(
        config.min_score_fraction * result.best_score);

    for (auto& chain : all_chains) {
        if (chain.score >= score_threshold) {
            result.chains.push_back(std::move(chain));
        } else {
            ++result.filtered_chains;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    result.chain_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return result;
}

ChainResult ExtractChains(
    std::span<const MinimizerHit> hits,
    uint32_t query_len,
    const ChainConfig& config)
{
    std::vector<Anchor> anchors;
    anchors.reserve(hits.size());

    for (const auto& hit : hits) {
        anchors.push_back({
            hit.ref_id,
            hit.ref_pos,
            hit.query_pos,
            hit.same_strand
        });
    }

    return ExtractChainsFromAnchors(anchors, query_len, config);
}

ChainResult ExtractChainsFromAnchorsWithScratch(
    std::span<const Anchor> anchors,
    uint32_t query_len,
    const ChainConfig& config,
    ChainScratch& scratch)
{
    auto start = std::chrono::high_resolution_clock::now();

    ChainResult result;
    result.total_anchors = anchors.size();

    if (anchors.empty()) {
        return result;
    }

    const size_t n = anchors.size();

    // Prepare scratch buffers (resize, don't reallocate if capacity sufficient)
    scratch.dp_score.resize(n);
    scratch.dp_pred.resize(n);
    scratch.used.resize_zero(n);  // Zero-initialize

    // Build keyed indices for grouping
    scratch.keyed_indices.clear();
    scratch.keyed_indices.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        uint64_t key = (static_cast<uint64_t>(anchors[i].ref_id) << 1) |
                       (anchors[i].same_strand ? 1 : 0);
        scratch.keyed_indices.push_back({key, i});
    }

    // Sort by (key, ref_pos)
    std::sort(scratch.keyed_indices.begin(), scratch.keyed_indices.end(),
        [&](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first < b.first;
            return anchors[a.second].ref_pos < anchors[b.second].ref_pos;
        });

    // Process each (ref_id, strand) group
    size_t group_start = 0;
    while (group_start < scratch.keyed_indices.size()) {
        uint64_t group_key = scratch.keyed_indices[group_start].first;
        size_t group_end = group_start + 1;
        while (group_end < scratch.keyed_indices.size() &&
               scratch.keyed_indices[group_end].first == group_key) {
            ++group_end;
        }

        scratch.group_indices.clear();
        for (size_t i = group_start; i < group_end; ++i) {
            scratch.group_indices.push_back(scratch.keyed_indices[i].second);
        }

        ChainSingleRefStrandScratch(
            anchors, scratch.group_indices.span(),
            scratch.dp_score, scratch.dp_pred, config);

        group_start = group_end;
    }

    // Build score order for chain extraction
    scratch.score_order.resize(n);
    for (size_t i = 0; i < n; ++i) {
        scratch.score_order[i] = i;
    }
    std::sort(scratch.score_order.begin(), scratch.score_order.end(),
        [&](size_t a, size_t b) {
            return scratch.dp_score[a] > scratch.dp_score[b];
        });

    std::vector<Chain> all_chains;

    for (size_t i = 0; i < scratch.score_order.size(); ++i) {
        size_t end_idx = scratch.score_order[i];
        if (scratch.used[end_idx]) continue;
        if (scratch.dp_score[end_idx] < config.min_chain_score) break;

        Chain chain = BacktrackChainScratch(
            anchors, scratch.dp_score, scratch.dp_pred,
            end_idx, scratch.backtrack);

        if (chain.NumAnchors() < config.min_chain_anchors) {
            continue;
        }

        for (uint32_t idx : chain.anchors) {
            scratch.used[idx] = true;
        }

        all_chains.push_back(std::move(chain));
    }

    if (!all_chains.empty()) {
        result.best_score = all_chains[0].score;
    }

    int32_t score_threshold = static_cast<int32_t>(
        config.min_score_fraction * result.best_score);

    for (auto& chain : all_chains) {
        if (chain.score >= score_threshold) {
            result.chains.push_back(std::move(chain));
        } else {
            ++result.filtered_chains;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    result.chain_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return result;
}

}  // namespace llmap::classical
