#include "self_interference/self_wavecollapse.h"
#include "self_interference/self_wavecollapse_impl.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace llmap::self_interference {

void SelfWaveCollapse::RunEM(
    InternalState& state,
    const SimilarityGraph& graph,
    std::span<const uint32_t> members)
{
    const size_t n = members.size();
    if (n == 0) return;

    state.iteration_count = 0;
    state.converged = false;

    for (size_t iter = 0; iter < config_.max_iterations; ++iter) {
        ++state.iteration_count;

        // E-step + M-step combined:
        // P_new(i, j) = (1-gamma) * P_old(i, j) + gamma * normalized_update(i, j)
        //
        // update(i, j) = S(i, j) * sum_k[P(k, j) * S(i, k)]
        // where S(i, j) = similarity_cache[i * n + j]

        for (size_t i = 0; i < n; ++i) {
            if (state.collapsed[i]) {
                // Already collapsed: keep distribution frozen
                for (size_t j = 0; j < n; ++j) {
                    state.prob_matrix_new[i * n + j] = state.prob_matrix[i * n + j];
                }
                continue;
            }

            // Compute update for read i
            float sum_update = 0.0f;
            std::vector<float> updates(n, 0.0f);

            for (size_t j = 0; j < n; ++j) {
                float s_ij = state.similarity_cache[i * n + j];

                // Sum over all reads k: P(k, j) * S(i, k)
                float support = 0.0f;
                for (size_t k = 0; k < n; ++k) {
                    float p_kj = state.prob_matrix[k * n + j];
                    float s_ik = state.similarity_cache[i * n + k];
                    support += p_kj * s_ik;
                }

                updates[j] = std::pow(s_ij, config_.similarity_exponent) * support;
                sum_update += updates[j];
            }

            // Normalize and apply damping
            if (sum_update > 0.0f) {
                for (size_t j = 0; j < n; ++j) {
                    float old_p = state.prob_matrix[i * n + j];
                    float new_p = updates[j] / sum_update;
                    state.prob_matrix_new[i * n + j] =
                        (1.0f - config_.gamma) * old_p + config_.gamma * new_p;
                }
            } else {
                // No update possible; keep old probabilities
                for (size_t j = 0; j < n; ++j) {
                    state.prob_matrix_new[i * n + j] = state.prob_matrix[i * n + j];
                }
            }

            // Renormalize to ensure sum = 1
            float row_sum = 0.0f;
            for (size_t j = 0; j < n; ++j) {
                row_sum += state.prob_matrix_new[i * n + j];
            }
            if (row_sum > 0.0f) {
                for (size_t j = 0; j < n; ++j) {
                    state.prob_matrix_new[i * n + j] /= row_sum;
                }
            }

            // Check for collapse
            float max_p = 0.0f;
            uint32_t max_j = 0;
            for (size_t j = 0; j < n; ++j) {
                if (state.prob_matrix_new[i * n + j] > max_p) {
                    max_p = state.prob_matrix_new[i * n + j];
                    max_j = static_cast<uint32_t>(j);
                }
            }

            if (max_p >= config_.collapse_threshold) {
                state.collapsed[i] = true;
                state.collapsed_to[i] = max_j;
            }
        }

        // Check convergence
        if (CheckConvergence(state.prob_matrix, state.prob_matrix_new)) {
            state.converged = true;
            std::swap(state.prob_matrix, state.prob_matrix_new);
            break;
        }

        // Swap buffers
        std::swap(state.prob_matrix, state.prob_matrix_new);
    }
}

float SelfWaveCollapse::ComputeSimilarityWeight(
    const SimilarityGraph& graph,
    uint32_t read_a,
    uint32_t read_b) const
{
    if (read_a == read_b) {
        return 1.0f;
    }

    // Get edge weight from similarity graph
    float edge_weight = graph.GetEdgeWeight(read_a, read_b);

    if (edge_weight > 0.0f) {
        return edge_weight;
    }

    // No direct edge: use 0 (no similarity contribution)
    return 0.0f;
}

bool SelfWaveCollapse::CheckConvergence(
    std::span<const float> old_probs,
    std::span<const float> new_probs) const
{
    if (old_probs.size() != new_probs.size()) {
        return false;
    }

    float max_diff = 0.0f;
    for (size_t i = 0; i < old_probs.size(); ++i) {
        float diff = std::abs(old_probs[i] - new_probs[i]);
        max_diff = std::max(max_diff, diff);
    }

    return max_diff < config_.convergence_threshold;
}

float SelfWaveCollapse::ComputeEntropy(std::span<const float> probs) {
    float entropy = 0.0f;
    for (float p : probs) {
        if (p > 0.0f) {
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

}  // namespace llmap::self_interference
