// LLmap — AllpairPipeline core implementation.
//
// Contains: AllpairResult, constructor/destructor, Run() orchestration.

#include "self_interference/allpair_pipeline.h"
#include "self_interference/allpair_pipeline_internal.h"

#include <algorithm>
#include <chrono>

namespace llmap::self_interference {

// ========== AllpairResult ==========

std::vector<uint32_t> AllpairResult::GetRepresentatives() const {
    std::vector<uint32_t> reps;
    for (size_t i = 0; i < is_representative.size(); ++i) {
        if (is_representative[i]) {
            reps.push_back(static_cast<uint32_t>(i));
        }
    }
    return reps;
}

// ========== AllpairPipeline ==========

AllpairPipeline::AllpairPipeline()
    : config_(), state_(std::make_unique<InternalState>()) {}

AllpairPipeline::AllpairPipeline(const AllpairConfig& config)
    : config_(config), state_(std::make_unique<InternalState>()) {}

AllpairPipeline::~AllpairPipeline() = default;

AllpairPipeline::AllpairPipeline(AllpairPipeline&&) noexcept = default;
AllpairPipeline& AllpairPipeline::operator=(AllpairPipeline&&) noexcept = default;

std::string AllpairPipeline::LastError() const {
    return state_ ? state_->last_error : "";
}

void AllpairPipeline::ReportProgress(const std::string& stage, size_t current,
                                      size_t total, const std::string& message) {
    if (progress_callback_) {
        progress_callback_(stage, current, total, message);
    }
    if (config_.verbose) {
        std::fprintf(stderr, "[%s] %zu/%zu: %s\n",
                     stage.c_str(), current, total, message.c_str());
    }
}

std::unique_ptr<AllpairResult> AllpairPipeline::Run() {
    return Run(nullptr);
}

std::unique_ptr<AllpairResult> AllpairPipeline::Run(ProgressCallback callback) {
    progress_callback_ = std::move(callback);
    state_->start_time = std::chrono::steady_clock::now();

    auto result = std::make_unique<AllpairResult>();

    // Stage 1: Load reads
    ReportProgress("load", 0, 6, "Loading reads from FASTQ");
    if (!LoadReads()) {
        result->success = false;
        result->error_message = state_->last_error;
        return result;
    }

    // Stage 2: Compute embeddings
    ReportProgress("embed", 1, 6, "Computing read embeddings");
    if (!ComputeEmbeddings()) {
        result->success = false;
        result->error_message = state_->last_error;
        return result;
    }

    // Stage 3: Build similarity graph
    ReportProgress("graph", 2, 6, "Building similarity graph");
    if (!BuildSimilarityGraph()) {
        result->success = false;
        result->error_message = state_->last_error;
        return result;
    }

    // Stage 4: Leiden clustering
    ReportProgress("cluster", 3, 6, "Running Leiden clustering");
    if (!RunClustering()) {
        result->success = false;
        result->error_message = state_->last_error;
        return result;
    }

    // Stage 5: Self-WaveCollapse refinement
    ReportProgress("refine", 4, 6, "Running Self-WaveCollapse refinement");
    if (!RunSelfWaveCollapse()) {
        result->success = false;
        result->error_message = state_->last_error;
        return result;
    }

    // Stage 6: Select representatives
    if (config_.select_representatives) {
        ReportProgress("reps", 5, 6, "Selecting cluster representatives");
        if (!SelectRepresentatives()) {
            result->success = false;
            result->error_message = state_->last_error;
            return result;
        }
    }

    // Stage 7: Write output
    ReportProgress("output", 6, 6, "Writing output");
    if (!WriteOutput()) {
        result->success = false;
        result->error_message = state_->last_error;
        return result;
    }

    // Populate result
    const size_t num_reads = state_->reads.size();
    result->cluster_ids.resize(num_reads);
    result->confidences.resize(num_reads);
    result->is_representative.resize(num_reads, false);
    result->read_ids.reserve(num_reads);

    for (size_t i = 0; i < num_reads; ++i) {
        result->read_ids.push_back(state_->reads[i].id);
    }

    // Copy from SWC result
    if (state_->swc_result) {
        for (const auto& assignment : state_->swc_result->assignments) {
            if (assignment.read_idx < num_reads) {
                result->cluster_ids[assignment.read_idx] = assignment.cluster_id;
                result->confidences[assignment.read_idx] = assignment.confidence;
            }
        }
        result->cluster_sizes = state_->swc_result->GetClusterSizes();
    }

    // Mark representatives
    if (state_->rep_result) {
        for (const auto& rep : state_->rep_result->representatives) {
            if (rep.read_idx < num_reads) {
                result->is_representative[rep.read_idx] = true;
            }
        }
        result->stats.num_representatives = state_->rep_result->representatives.size();
    }

    // Calculate final statistics
    auto end_time = std::chrono::steady_clock::now();
    result->stats.total_time_ms =
        std::chrono::duration<float, std::milli>(end_time - state_->start_time).count();
    result->stats.input_reads = num_reads;
    result->stats.reads_per_second = num_reads * 1000.0f / result->stats.total_time_ms;

    if (state_->graph) {
        result->stats.graph_nodes = state_->graph->NumNodes();
        result->stats.graph_edges = state_->graph->NumEdges();
    }

    if (state_->clustering) {
        result->stats.num_clusters = state_->clustering->num_communities;
        result->stats.modularity = state_->clustering->modularity;
    }

    if (state_->swc_result) {
        result->stats.reads_collapsed = state_->swc_result->stats.reads_collapsed;
    }

    result->success = true;
    return result;
}

// ========== Convenience function ==========

std::unique_ptr<AllpairResult> RunAllpairPipeline(
    const std::filesystem::path& reads,
    const std::filesystem::path& output,
    const std::filesystem::path& model_path) {

    AllpairConfig config;
    config.reads_a = reads;
    config.output = output;
    config.model_path = model_path;

    AllpairPipeline pipeline(config);
    return pipeline.Run();
}

}  // namespace llmap::self_interference
