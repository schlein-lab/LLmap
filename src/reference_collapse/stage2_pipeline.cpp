// LLmap — Stage2Pipeline implementation.

#include "reference_collapse/stage2_pipeline.h"
#include "reference_collapse/collapse_check.h"
#include "reference_collapse/em_iterator.h"
#include "reference_collapse/member_propagation.h"
#include "reference_collapse/reference_index.h"
#include "reference_collapse/refinement.h"
#include "core/wave_state.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <chrono>

namespace llmap {

// Internal state for Stage 2 pipeline
struct Stage2Pipeline::InternalState {
    std::unique_ptr<ReferenceIndex> ref_index;
    std::unique_ptr<WaveState> wave_state;
    std::unique_ptr<Stage2Result> result;

    EmIterator em_iterator;
    CollapseChecker collapse_checker;
    Refinement refinement;
    MemberPropagation member_propagation;

    ChildIndex l0_to_l1_index;
    ChildIndex l1_to_l2_index;

    Stage2Stats stats;
    std::string last_error;

    std::chrono::steady_clock::time_point start_time;
};

Stage2Pipeline::Stage2Pipeline() : Stage2Pipeline(Stage2Config{}) {}

Stage2Pipeline::Stage2Pipeline(const Stage2Config& config)
    : config_(config), state_(std::make_unique<InternalState>()) {}

Stage2Pipeline::~Stage2Pipeline() = default;

Stage2Pipeline::Stage2Pipeline(Stage2Pipeline&&) noexcept = default;
Stage2Pipeline& Stage2Pipeline::operator=(Stage2Pipeline&&) noexcept = default;

std::unique_ptr<Stage2Result> Stage2Pipeline::Run(
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph) {
    return Run(clustering, rep_result, graph, nullptr);
}

std::unique_ptr<Stage2Result> Stage2Pipeline::Run(
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph,
    Stage2ProgressCallback callback) {

    progress_callback_ = std::move(callback);
    state_->start_time = std::chrono::steady_clock::now();
    state_->result = std::make_unique<Stage2Result>();
    state_->stats = {};

    // Load reference index
    ReportProgress("init", 0, 6, "Loading reference index");
    if (!LoadReferenceIndex()) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    // Initialize WaveState from representatives
    ReportProgress("init", 1, 6, "Initializing WaveState");
    if (!InitializeWaveState(rep_result)) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    // Run EM at L0 level
    ReportProgress("em_l0", 2, 6, "Running EM at L0 level");
    if (!RunEMLevel(WaveLevel::L0)) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    // Refine L0 → L1 and run EM at L1
    ReportProgress("em_l1", 3, 6, "Refining to L1 and running EM");
    if (!RefineLevel(WaveLevel::L0) || !RunEMLevel(WaveLevel::L1)) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    // Refine L1 → L2 and run EM at L2
    ReportProgress("em_l2", 4, 6, "Refining to L2 and running EM");
    if (!RefineLevel(WaveLevel::L1) || !RunEMLevel(WaveLevel::L2)) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    // Propagate positions to cluster members
    ReportProgress("propagate", 5, 6, "Propagating to cluster members");
    if (!PropagateToMembers(clustering, rep_result, graph)) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    // Extract final alignments
    ReportProgress("finalize", 6, 6, "Extracting alignments");
    if (!ExtractAlignments()) {
        state_->result->success = false;
        state_->result->error_message = state_->last_error;
        return std::move(state_->result);
    }

    auto end = std::chrono::steady_clock::now();
    state_->stats.total_time_ms =
        std::chrono::duration<float, std::milli>(end - state_->start_time).count();

    state_->result->stats = state_->stats;
    state_->result->success = true;
    return std::move(state_->result);
}

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

void Stage2Pipeline::SetReferenceIndex(std::unique_ptr<ReferenceIndex> index) {
    state_->ref_index = std::move(index);
}

std::string Stage2Pipeline::LastError() const {
    return state_->last_error;
}

const WaveState* Stage2Pipeline::GetWaveState() const {
    return state_->wave_state.get();
}

void Stage2Pipeline::ReportProgress(const std::string& stage, std::size_t current,
                                    std::size_t total, const std::string& message) {
    if (progress_callback_) {
        progress_callback_(stage, current, total, message);
    }
}

// --- Result accessors ---

std::vector<ReadAlignment> Stage2Result::GetAlignmentsForRead(std::uint32_t read_idx) const {
    std::vector<ReadAlignment> result;
    for (const auto& a : alignments) {
        if (a.read_idx == read_idx) {
            result.push_back(a);
        }
    }
    return result;
}

std::vector<ReadAlignment> Stage2Result::GetHighConfidenceAlignments(float min_confidence) const {
    std::vector<ReadAlignment> result;
    for (const auto& a : alignments) {
        if (a.confidence >= min_confidence) {
            result.push_back(a);
        }
    }
    return result;
}

// --- Convenience function ---

std::unique_ptr<Stage2Result> RunStage2Pipeline(
    const std::filesystem::path& index_path,
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph) {

    Stage2Config config;
    config.index_path = index_path;

    Stage2Pipeline pipeline(config);
    return pipeline.Run(clustering, rep_result, graph);
}

}  // namespace llmap
