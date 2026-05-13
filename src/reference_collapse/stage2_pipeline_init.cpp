// LLmap — Stage2Pipeline initialization routines.

#include "reference_collapse/stage2_pipeline.h"
#include "reference_collapse/stage2_pipeline_impl.h"
#include "self_interference/cluster_rep.h"

namespace llmap {

bool Stage2Pipeline::LoadReferenceIndex() {
    if (state_->ref_index) {
        return true;  // Already loaded via SetReferenceIndex
    }

    if (config_.index_path.empty()) {
        state_->last_error = "No reference index path specified";
        return false;
    }

    state_->ref_index = ReferenceIndex::Load(config_.index_path);
    if (!state_->ref_index) {
        state_->last_error = "Failed to load reference index: " + config_.index_path.string();
        return false;
    }

    // Build child indices for refinement
    state_->l0_to_l1_index = ChildIndex::BuildL0ToL1(state_->ref_index->Pyramid());
    state_->l1_to_l2_index = ChildIndex::BuildL1ToL2(state_->ref_index->Pyramid());

    state_->stats.ref_targets = state_->ref_index->NumTargets();
    state_->stats.ref_l2_buckets = state_->ref_index->Pyramid().l2_count();

    return true;
}

bool Stage2Pipeline::InitializeWaveState(
    const self_interference::ClusterRepResult& rep_result) {

    if (!state_->ref_index) {
        state_->last_error = "Reference index not loaded";
        return false;
    }

    auto num_reps = static_cast<std::uint32_t>(rep_result.representatives.size());
    state_->stats.num_representatives = num_reps;

    // Initialize WaveState for all representatives at L0 level
    auto num_l0 = state_->ref_index->Pyramid().l0_count();
    state_->wave_state = std::make_unique<WaveState>(num_reps, WaveLevel::L0);

    // Initialize uniform distribution over all L0 buckets for each rep
    auto prob = 1.0f / static_cast<float>(num_l0);
    std::vector<BucketProb> candidates;
    candidates.reserve(num_l0);
    for (std::uint32_t b = 0; b < num_l0; ++b) {
        candidates.push_back({b, prob});
    }

    for (std::uint32_t i = 0; i < num_reps; ++i) {
        state_->wave_state->set_read_candidates(i, candidates);
    }

    // Configure EM iterator
    EmIteratorConfig em_config;
    em_config.gamma = config_.gamma;
    em_config.tau_collapse = config_.tau_collapse;
    em_config.apply_smoothing = config_.apply_smoothing;
    state_->em_iterator.SetConfig(em_config);

    // Configure collapse checker
    CollapseCheckConfig collapse_config;
    collapse_config.tau_collapse = config_.tau_collapse;
    state_->collapse_checker.SetConfig(collapse_config);

    // Configure refinement
    RefinementConfig refine_config;
    refine_config.expansion_threshold = config_.expansion_threshold;
    refine_config.max_candidates = config_.max_candidates;
    state_->refinement.SetConfig(refine_config);

    // Configure member propagation
    MemberPropagationConfig prop_config;
    prop_config.base_confidence_scaling = config_.propagation_confidence;
    prop_config.similarity_weight = config_.similarity_weight;
    state_->member_propagation.SetConfig(prop_config);

    return true;
}

}  // namespace llmap
