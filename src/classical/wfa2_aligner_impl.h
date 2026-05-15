#pragma once

#include "classical/wfa2_aligner.h"

#include <algorithm>
#include <limits>
#include <vector>

namespace llmap::classical {

// PIMPL implementation for WFA2Aligner
// Contains either WFA2-lib state or fallback DP state
class WFA2AlignerImpl {
public:
    explicit WFA2AlignerImpl(const WFA2Config& config) : config_(config) {}

#ifdef LLMAP_HAS_WFA2
    // WFA2-lib state would be here
    void* wfa_aligner = nullptr;
#endif

    WFA2Config config_;
};

// Fallback gap-affine DP alignment (Gotoh algorithm)
// Used when WFA2-lib is not available
class FallbackAligner {
public:
    explicit FallbackAligner(const WFA2Config& config) : config_(config) {}

    std::optional<WFA2Result> Align(
        std::string_view query,
        std::string_view reference) const {

        if (query.empty() || reference.empty()) {
            return std::nullopt;
        }

        if (query.size() > config_.max_pattern_length ||
            reference.size() > config_.max_text_length) {
            return std::nullopt;
        }

        const int32_t m = static_cast<int32_t>(query.size());
        const int32_t n = static_cast<int32_t>(reference.size());

        // DP matrices: M[i][j] = best score ending with match/mismatch
        //              I[i][j] = best score ending with insertion (gap in ref)
        //              D[i][j] = best score ending with deletion (gap in query)
        // Using 1D arrays with row-major indexing for cache efficiency
        const int32_t INF = std::numeric_limits<int32_t>::max() / 2;
        const size_t stride = n + 1;

        std::vector<int32_t> M((m + 1) * stride, INF);
        std::vector<int32_t> I((m + 1) * stride, INF);
        std::vector<int32_t> D((m + 1) * stride, INF);

        // Initialize: semi-global alignment allows free gaps at reference ends.
        // First column: query gaps are penalized (must align entire query).
        // First row: reference gaps are FREE (can start alignment anywhere in ref).
        M[0] = 0;
        for (int32_t j = 1; j <= n; ++j) {
            D[j] = 0;  // Free gap at reference start (semi-global)
            M[j] = INF;
        }
        for (int32_t i = 1; i <= m; ++i) {
            size_t idx = i * stride;
            I[idx] = config_.gap_open + (i - 1) * config_.gap_extend;
            M[idx] = INF;
        }

        // Fill DP matrices
        for (int32_t i = 1; i <= m; ++i) {
            for (int32_t j = 1; j <= n; ++j) {
                size_t idx = i * stride + j;
                size_t diag = (i - 1) * stride + (j - 1);
                size_t up = (i - 1) * stride + j;
                size_t left = i * stride + (j - 1);

                // Match/mismatch
                int32_t sub_cost = (query[i - 1] == reference[j - 1])
                    ? config_.match_score
                    : config_.mismatch_penalty;

                int32_t m_score = std::min({M[diag], I[diag], D[diag]}) + sub_cost;

                // Insertion (gap in reference)
                int32_t i_extend = I[up] + config_.gap_extend;
                int32_t i_open = M[up] + config_.gap_open;
                I[idx] = std::min(i_extend, i_open);

                // Deletion (gap in query)
                int32_t d_extend = D[left] + config_.gap_extend;
                int32_t d_open = M[left] + config_.gap_open;
                D[idx] = std::min(d_extend, d_open);

                M[idx] = m_score;
            }
        }

        // Semi-global: find best endpoint in last row (allow free gaps at ref end).
        // The alignment must consume entire query but can end anywhere in reference.
        int32_t best_score = INF;
        int32_t best_j = n;
        for (int32_t jj = 1; jj <= n; ++jj) {
            size_t idx = m * stride + jj;
            int32_t score = std::min({M[idx], I[idx], D[idx]});
            if (score < best_score) {
                best_score = score;
                best_j = jj;
            }
        }

        if (best_score > config_.max_alignment_score) {
            return std::nullopt;
        }

        // Traceback to build CIGAR
        std::vector<CigarElement> cigar;
        int32_t i = m, j = best_j;

        // Determine which matrix we're starting from
        size_t end_idx = m * stride + best_j;
        enum class State { M_STATE, I_STATE, D_STATE };
        State state;
        if (M[end_idx] <= I[end_idx] && M[end_idx] <= D[end_idx]) {
            state = State::M_STATE;
        } else if (I[end_idx] <= D[end_idx]) {
            state = State::I_STATE;
        } else {
            state = State::D_STATE;
        }

        // Count statistics
        size_t matches = 0, mismatches = 0, insertions = 0, deletions = 0;
        int32_t ref_start_pos = 0;  // Track where alignment starts in reference

        while (i > 0 || j > 0) {
            size_t idx = i * stride + j;

            if (i == 0) {
                // Semi-global: remaining ref bases are free gaps at start.
                // Record where alignment actually starts.
                ref_start_pos = j;
                break;
            }
            if (j == 0) {
                // Must insert into reference (query not fully aligned yet)
                AddCigarOp(cigar, CigarOp::Insertion, i);
                insertions += i;
                break;
            }

            // Check if we reached a free-gap entry point (D[j] was initialized to 0)
            // If D[idx] == 0 and we came from D state at first row, we can stop.
            if (state == State::D_STATE && D[idx] == 0 && i == 0) {
                ref_start_pos = j;
                break;
            }

            size_t diag = (i - 1) * stride + (j - 1);
            size_t up = (i - 1) * stride + j;
            size_t left = i * stride + (j - 1);

            switch (state) {
                case State::M_STATE: {
                    bool match = (query[i - 1] == reference[j - 1]);
                    if (match) {
                        AddCigarOp(cigar, CigarOp::Equal, 1);
                        ++matches;
                    } else {
                        AddCigarOp(cigar, CigarOp::Diff, 1);
                        ++mismatches;
                    }

                    // Determine predecessor
                    int32_t sub_cost = match ? config_.match_score : config_.mismatch_penalty;
                    int32_t expected = M[idx] - sub_cost;

                    if (M[diag] == expected) {
                        state = State::M_STATE;
                    } else if (I[diag] == expected) {
                        state = State::I_STATE;
                    } else {
                        state = State::D_STATE;
                    }
                    --i; --j;
                    break;
                }
                case State::I_STATE: {
                    AddCigarOp(cigar, CigarOp::Insertion, 1);
                    ++insertions;

                    // Check if we extended or opened
                    if (I[up] + config_.gap_extend == I[idx]) {
                        state = State::I_STATE;
                    } else {
                        state = State::M_STATE;
                    }
                    --i;
                    break;
                }
                case State::D_STATE: {
                    AddCigarOp(cigar, CigarOp::Deletion, 1);
                    ++deletions;

                    if (D[left] + config_.gap_extend == D[idx]) {
                        state = State::D_STATE;
                    } else {
                        state = State::M_STATE;
                    }
                    --j;
                    break;
                }
            }
        }

        // Reverse CIGAR (built backwards)
        std::reverse(cigar.begin(), cigar.end());

        // Merge adjacent same-type operations
        MergeCigar(cigar);

        WFA2Result result;
        result.cigar = std::move(cigar);
        result.score = best_score;
        result.query_start = 0;
        result.query_end = m;
        result.ref_start = ref_start_pos;  // Semi-global: actual start from traceback
        result.ref_end = best_j;  // Semi-global: actual end where alignment ends
        result.num_matches = matches;
        result.num_mismatches = mismatches;
        result.num_insertions = insertions;
        result.num_deletions = deletions;

        size_t aligned = matches + mismatches + insertions + deletions;
        result.identity = aligned > 0
            ? static_cast<float>(matches) / static_cast<float>(aligned)
            : 0.0f;

        return result;
    }

private:
    WFA2Config config_;

    static void AddCigarOp(std::vector<CigarElement>& cigar, CigarOp op, uint32_t len) {
        if (!cigar.empty() && cigar.back().op == op) {
            cigar.back().length += len;
        } else {
            cigar.push_back({op, len});
        }
    }

    static void MergeCigar(std::vector<CigarElement>& cigar) {
        if (cigar.size() <= 1) return;

        std::vector<CigarElement> merged;
        merged.push_back(cigar[0]);

        for (size_t i = 1; i < cigar.size(); ++i) {
            if (cigar[i].op == merged.back().op) {
                merged.back().length += cigar[i].length;
            } else {
                merged.push_back(cigar[i]);
            }
        }

        cigar = std::move(merged);
    }
};

}  // namespace llmap::classical
