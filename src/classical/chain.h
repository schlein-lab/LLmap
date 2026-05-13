#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace llmap::classical {

// Forward declarations
struct MinimizerHit;

// Configuration for chaining algorithm
struct ChainConfig {
    // Scoring parameters (matches minimap2 defaults)
    int32_t match_score = 1;
    int32_t gap_penalty = 1;      // Per-base gap cost (applied to |gap| - 1)

    // Chain length constraints
    uint32_t max_gap_ref = 5000;      // Max gap in reference coordinates
    uint32_t max_gap_query = 5000;    // Max gap in query coordinates
    uint32_t min_chain_anchors = 3;   // Min anchors to form valid chain

    // Score filtering
    float min_score_fraction = 0.9f;  // Keep chains ≥ this × best_score
    int32_t min_chain_score = 20;     // Absolute minimum chain score

    // DP optimization
    uint32_t max_skip = 25;           // Max predecessors to skip in DP
    uint32_t bandwidth = 500;         // DP bandwidth for colinear extension
};

// An anchor for chaining (derived from MinimizerHit)
struct Anchor {
    uint32_t ref_id;      // Reference sequence ID
    uint32_t ref_pos;     // Position in reference
    uint32_t query_pos;   // Position in query
    bool same_strand;     // True if same strand as query

    bool operator<(const Anchor& o) const {
        if (ref_id != o.ref_id) return ref_id < o.ref_id;
        return ref_pos < o.ref_pos;
    }
};

// A chain of colinear anchors
struct Chain {
    uint32_t ref_id;                 // Reference sequence containing this chain
    bool is_forward;                 // Chain orientation (true = forward strand)
    int32_t score;                   // Chain score from DP
    uint32_t ref_start;              // Start position in reference
    uint32_t ref_end;                // End position in reference
    uint32_t query_start;            // Start position in query
    uint32_t query_end;              // End position in query
    std::vector<uint32_t> anchors;   // Indices into original anchor array

    uint32_t NumAnchors() const { return static_cast<uint32_t>(anchors.size()); }
    uint32_t RefSpan() const { return ref_end - ref_start; }
    uint32_t QuerySpan() const { return query_end - query_start; }
};

// Result of chain extraction
struct ChainResult {
    std::vector<Chain> chains;       // All chains passing filters, sorted by score desc
    int32_t best_score = 0;          // Score of best chain
    size_t total_anchors = 0;        // Total input anchors
    size_t filtered_chains = 0;      // Chains removed by score filter
    float chain_time_ms = 0.0f;      // Time for chaining
};

// Extract chains from minimizer hits
//
// This implements the colinear chaining algorithm from minimap2:
// 1. Sort anchors by (ref_id, ref_pos)
// 2. For each ref_id and strand, run sparse DP
// 3. Backtrack to extract chains
// 4. Filter chains by score relative to best
//
// Parameters:
//   hits: minimizer hits from index query
//   query_len: length of query sequence
//   config: chaining parameters
//
// Returns: ChainResult with all chains passing filters
ChainResult ExtractChains(
    std::span<const MinimizerHit> hits,
    uint32_t query_len,
    const ChainConfig& config = {});

// Extract chains from raw anchors (internal/testing)
ChainResult ExtractChainsFromAnchors(
    std::span<const Anchor> anchors,
    uint32_t query_len,
    const ChainConfig& config = {});

// Compute chain score between two anchors
// Returns negative if anchors are not colinear
int32_t AnchorPairScore(
    const Anchor& a,
    const Anchor& b,
    const ChainConfig& config);

// Check if two anchors are colinear (compatible for chaining)
bool IsColinear(
    const Anchor& a,
    const Anchor& b,
    const ChainConfig& config);

}  // namespace llmap::classical
