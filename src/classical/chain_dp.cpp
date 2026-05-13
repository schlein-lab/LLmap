#include "classical/chain.h"
#include "classical/minimizer_index.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <numeric>

namespace llmap::classical {

namespace {

// DP state for a single anchor
struct DpState {
    int32_t score = 0;       // Best chain score ending here
    int32_t pred = -1;       // Predecessor index (-1 = chain start)
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

    // Initialize DP
    for (size_t idx : sorted_indices) {
        dp[idx].score = config.match_score;  // Single anchor score
        dp[idx].pred = -1;
    }

    // Main DP loop
    for (size_t i = 1; i < n; ++i) {
        const size_t curr_idx = sorted_indices[i];
        const Anchor& curr = anchors[curr_idx];

        // Look back at previous anchors
        size_t skip_count = 0;
        for (size_t j = i; j > 0 && skip_count < config.max_skip; ) {
            --j;
            const size_t prev_idx = sorted_indices[j];
            const Anchor& prev = anchors[prev_idx];

            // Gap check in reference coordinates
            if (curr.ref_pos - prev.ref_pos > config.max_gap_ref) {
                break;  // Too far back, no more valid predecessors
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
    // Convert hits to anchors
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

}  // namespace llmap::classical
