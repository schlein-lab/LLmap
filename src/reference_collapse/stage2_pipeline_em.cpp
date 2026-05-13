// LLmap — Stage2Pipeline EM iteration and extraction.

#include "reference_collapse/stage2_pipeline.h"
#include "reference_collapse/stage2_pipeline_impl.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <chrono>

namespace llmap {

bool Stage2Pipeline::RunEMLevel(WaveLevel level) {
    if (!state_->wave_state) {
        state_->last_error = "WaveState not initialized";
        return false;
    }

    auto start = std::chrono::steady_clock::now();

    LevelStats* level_stats = nullptr;
    std::uint32_t num_buckets = 0;

    switch (level) {
        case WaveLevel::L0:
            level_stats = &state_->stats.l0_stats;
            num_buckets = static_cast<std::uint32_t>(state_->ref_index->Pyramid().l0_count());
            break;
        case WaveLevel::L1:
            level_stats = &state_->stats.l1_stats;
            num_buckets = static_cast<std::uint32_t>(state_->ref_index->Pyramid().l1_count());
            break;
        case WaveLevel::L2:
            level_stats = &state_->stats.l2_stats;
            num_buckets = static_cast<std::uint32_t>(state_->ref_index->Pyramid().l2_count());
            break;
        default:
            state_->last_error = "Invalid wave level";
            return false;
    }

    // EM iteration loop
    for (std::uint32_t iter = 0; iter < config_.max_iterations_per_level; ++iter) {
        auto iter_stats = state_->em_iterator.Step(*state_->wave_state, num_buckets);
        auto collapse_stats = state_->collapse_checker.Check(*state_->wave_state);

        level_stats->iterations = iter + 1;
        level_stats->reads_processed = iter_stats.reads_processed;
        level_stats->reads_collapsed = collapse_stats.already_collapsed + collapse_stats.newly_collapsed;
        level_stats->avg_entropy = iter_stats.mean_entropy;

        // Check for convergence
        if (iter_stats.max_prob_delta < config_.convergence_delta) {
            break;
        }

        // Check if we should refine to next level
        float collapse_rate = collapse_stats.collapse_rate();
        if (collapse_rate >= config_.refine_trigger_rate && level != WaveLevel::L2) {
            break;
        }
    }

    auto end = std::chrono::steady_clock::now();
    level_stats->time_ms = std::chrono::duration<float, std::milli>(end - start).count();
    state_->stats.total_iterations += level_stats->iterations;

    return true;
}

bool Stage2Pipeline::RefineLevel(WaveLevel from_level) {
    if (!state_->wave_state) {
        state_->last_error = "WaveState not initialized";
        return false;
    }

    RefinementStats refine_stats;

    if (from_level == WaveLevel::L0) {
        refine_stats = state_->refinement.RefineL0ToL1(
            *state_->wave_state, state_->l0_to_l1_index);
        state_->stats.l0_stats.reads_refined = refine_stats.reads_refined;
    } else if (from_level == WaveLevel::L1) {
        refine_stats = state_->refinement.RefineL1ToL2(
            *state_->wave_state, state_->l1_to_l2_index);
        state_->stats.l1_stats.reads_refined = refine_stats.reads_refined;
    } else {
        state_->last_error = "Cannot refine from L2 level";
        return false;
    }

    return true;
}

bool Stage2Pipeline::PropagateToMembers(
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph) {

    auto start = std::chrono::steady_clock::now();

    auto prop_result = state_->member_propagation.Propagate(
        *state_->wave_state, clustering, rep_result, graph);

    state_->stats.members_propagated = prop_result.stats.num_members_propagated;
    state_->stats.num_clusters = prop_result.stats.num_clusters;
    state_->stats.total_reads = prop_result.stats.num_members_propagated +
                                prop_result.stats.num_representatives;

    auto end = std::chrono::steady_clock::now();
    state_->stats.propagation_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return true;
}

bool Stage2Pipeline::ExtractAlignments() {
    if (!state_->wave_state || !state_->ref_index) {
        state_->last_error = "State not properly initialized";
        return false;
    }

    auto& result = state_->result;
    auto& ws = *state_->wave_state;
    auto l2_buckets = state_->ref_index->Pyramid().l2_buckets();

    std::uint32_t collapsed_count = 0;

    for (std::uint32_t read_idx = 0; read_idx < ws.n_reads(); ++read_idx) {
        if (ws.is_collapsed(read_idx)) {
            auto bucket_id = ws.collapsed_bucket(read_idx);
            if (bucket_id < l2_buckets.size()) {
                const auto& bucket = l2_buckets[bucket_id];

                // Find max probability from probabilities for this read
                auto probs = ws.probabilities_for_read(read_idx);
                float max_prob = probs.empty() ? 0.0f : *std::max_element(probs.begin(), probs.end());

                ReadAlignment alignment;
                alignment.read_idx = read_idx;
                alignment.target_name = bucket.target_id;
                alignment.position = bucket.start;
                alignment.confidence = max_prob;
                alignment.is_representative = (read_idx < state_->stats.num_representatives);
                alignment.cluster_id = 0;  // Would need to track this
                alignment.l2_bucket_id = bucket_id;

                result->alignments.push_back(alignment);
                ++collapsed_count;
            }
        }
    }

    state_->stats.collapse_rate =
        static_cast<float>(collapsed_count) / static_cast<float>(ws.n_reads());

    return true;
}

}  // namespace llmap
