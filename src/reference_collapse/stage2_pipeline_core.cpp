// LLmap — Stage2Pipeline core orchestration and accessors.

#include "reference_collapse/stage2_pipeline.h"
#include "reference_collapse/stage2_pipeline_impl.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <chrono>

namespace llmap {

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
