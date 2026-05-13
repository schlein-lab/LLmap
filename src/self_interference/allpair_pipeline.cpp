// LLmap — AllpairPipeline implementation.

#include "self_interference/allpair_pipeline.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>

#include "ai/foundation_embedder.h"
#include "io/fastq_reader.h"
#include "output/cluster_writer.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/faiss_wrapper.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"

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

// ========== Internal State ==========

struct AllpairPipeline::InternalState {
    // Loaded reads
    std::vector<io::FastqRecord> reads;

    // Computed embeddings (flattened: [num_reads * embedding_dim])
    std::vector<float> embeddings;

    // FAISS index
    std::unique_ptr<FaissIndex> faiss_index;

    // Similarity graph
    std::unique_ptr<SimilarityGraph> graph;

    // Clustering result
    std::unique_ptr<ClusteringResult> clustering;

    // Self-WaveCollapse result
    std::unique_ptr<SelfWaveCollapseResult> swc_result;

    // Representative selection result
    std::unique_ptr<ClusterRepResult> rep_result;

    // Timing
    std::chrono::steady_clock::time_point start_time;

    // Error message
    std::string last_error;
};

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

bool AllpairPipeline::LoadReads() {
    io::FastqReaderConfig reader_config;
    reader_config.min_length = config_.min_read_length;
    reader_config.max_length = config_.max_read_length;
    reader_config.max_records = config_.max_reads;

    auto reader = io::FastqReader::Open(config_.reads_a, reader_config);
    if (!reader) {
        state_->last_error = "Failed to open FASTQ file: " + config_.reads_a.string();
        return false;
    }

    state_->reads = reader->ReadAll();

    if (state_->reads.empty()) {
        state_->last_error = "No reads loaded from FASTQ file";
        return false;
    }

    auto stats = reader->GetStats();
    if (config_.verbose) {
        std::fprintf(stderr, "Loaded %zu reads, avg length %.1f bp\n",
                     stats.total_records, stats.avg_length);
    }

    return true;
}

bool AllpairPipeline::ComputeEmbeddings() {
    // If no model path provided, use placeholder embeddings for testing
    if (config_.model_path.empty() || !std::filesystem::exists(config_.model_path)) {
        // Generate deterministic pseudo-embeddings based on sequence content
        const size_t num_reads = state_->reads.size();
        const size_t dim = config_.embedding_dim;
        state_->embeddings.resize(num_reads * dim);

        for (size_t i = 0; i < num_reads; ++i) {
            const auto& seq = state_->reads[i].sequence;

            // Simple hash-based embedding for testing
            std::hash<std::string> hasher;
            size_t h = hasher(seq);

            for (size_t d = 0; d < dim; ++d) {
                h = h * 2654435761u + d;
                state_->embeddings[i * dim + d] =
                    static_cast<float>(h % 1000) / 500.0f - 1.0f;
            }

            // Normalize
            float norm = 0.0f;
            for (size_t d = 0; d < dim; ++d) {
                norm += state_->embeddings[i * dim + d] * state_->embeddings[i * dim + d];
            }
            norm = std::sqrt(norm);
            if (norm > 0) {
                for (size_t d = 0; d < dim; ++d) {
                    state_->embeddings[i * dim + d] /= norm;
                }
            }
        }

        return true;
    }

    // Use real embedder
    ai::EmbedderConfig embedder_config;
    embedder_config.model_path = config_.model_path;
    embedder_config.embedding_dim = config_.embedding_dim;
    embedder_config.batch_size = config_.embedding_batch_size;

    if (config_.use_gpu_embedder) {
        embedder_config.provider = ai::ExecutionProvider::CUDA;
        embedder_config.device_id = config_.gpu_device_id;
    }

    auto embedder = ai::FoundationEmbedder::Create(embedder_config);
    if (!embedder) {
        state_->last_error = "Failed to create embedder";
        return false;
    }

    const size_t num_reads = state_->reads.size();
    const size_t dim = embedder->EmbeddingDim();
    state_->embeddings.resize(num_reads * dim);

    // Embed in batches
    std::vector<std::string_view> batch_seqs;
    batch_seqs.reserve(config_.embedding_batch_size);

    for (size_t i = 0; i < num_reads; i += config_.embedding_batch_size) {
        batch_seqs.clear();
        const size_t batch_end = std::min(i + config_.embedding_batch_size, num_reads);

        for (size_t j = i; j < batch_end; ++j) {
            batch_seqs.push_back(state_->reads[j].sequence);
        }

        std::span<float> output_span(
            state_->embeddings.data() + i * dim,
            (batch_end - i) * dim);

        if (!embedder->EmbedBatchInto(batch_seqs, output_span)) {
            state_->last_error = "Embedding failed";
            return false;
        }
    }

    return true;
}

bool AllpairPipeline::BuildSimilarityGraph() {
    const size_t num_reads = state_->reads.size();

    if (num_reads < 2) {
        state_->last_error = "Need at least 2 reads to build similarity graph";
        return false;
    }

    // Create FAISS index
    FaissIndexConfig faiss_config;
    faiss_config.embedding_dim = config_.embedding_dim;
    faiss_config.nlist = std::min(config_.faiss_nlist, num_reads / 10 + 1);
    faiss_config.nprobe = config_.faiss_nprobe;
    faiss_config.default_k = config_.faiss_k;
    faiss_config.normalize_vectors = true;

    if (config_.use_gpu_faiss) {
        faiss_config.provider = FaissProvider::GPU;
        faiss_config.gpu_device_id = config_.gpu_device_id;
    }

    // For small datasets, use exact search
    if (num_reads < 10000) {
        faiss_config.index_type = FaissIndexType::FlatL2;
    } else {
        faiss_config.index_type = FaissIndexType::IVFFlat;
    }

    state_->faiss_index = FaissIndex::Create(faiss_config);
    if (!state_->faiss_index) {
        state_->last_error = "Failed to create FAISS index";
        return false;
    }

    // Train and add vectors
    if (!state_->faiss_index->TrainAndAdd(state_->embeddings, num_reads)) {
        state_->last_error = "Failed to build FAISS index: " +
                             state_->faiss_index->LastError();
        return false;
    }

    // Search for k-NN
    const size_t k = std::min(config_.faiss_k, num_reads - 1);
    auto knn_result = state_->faiss_index->SearchBatch(
        state_->embeddings, num_reads, k);

    if (!knn_result) {
        state_->last_error = "k-NN search failed";
        return false;
    }

    // Build similarity graph from k-NN results
    SimilarityGraphConfig graph_config;
    graph_config.convert_distance_to_similarity = true;
    graph_config.make_symmetric = true;
    graph_config.remove_self_loops = true;

    state_->graph = SimilarityGraph::BuildFromKnn(
        knn_result->indices,
        knn_result->distances,
        num_reads,
        k,
        graph_config);

    if (!state_->graph) {
        state_->last_error = "Failed to build similarity graph";
        return false;
    }

    return true;
}

bool AllpairPipeline::RunClustering() {
    if (!state_->graph) {
        state_->last_error = "Similarity graph not built";
        return false;
    }

    LeidenConfig leiden_config;
    leiden_config.resolution = config_.leiden_resolution;
    leiden_config.max_iterations = config_.leiden_max_iterations;
    leiden_config.min_community_size = config_.min_cluster_size;

    LeidenClustering leiden(leiden_config);
    state_->clustering = leiden.Cluster(*state_->graph);

    if (!state_->clustering) {
        state_->last_error = "Leiden clustering failed";
        return false;
    }

    if (config_.verbose) {
        std::fprintf(stderr, "Found %zu clusters, modularity = %.4f\n",
                     state_->clustering->num_communities,
                     state_->clustering->modularity);
    }

    return true;
}

bool AllpairPipeline::RunSelfWaveCollapse() {
    if (!state_->graph || !state_->clustering) {
        state_->last_error = "Graph or clustering not available";
        return false;
    }

    SelfWaveCollapseConfig swc_config;
    swc_config.max_iterations = config_.swc_max_iterations;
    swc_config.convergence_threshold = config_.swc_convergence_threshold;
    swc_config.collapse_threshold = config_.swc_collapse_threshold;
    swc_config.gamma = config_.swc_gamma;

    SelfWaveCollapse swc(swc_config);
    state_->swc_result = swc.Refine(*state_->graph, *state_->clustering);

    if (!state_->swc_result) {
        state_->last_error = "Self-WaveCollapse refinement failed";
        return false;
    }

    if (config_.verbose) {
        std::fprintf(stderr, "SWC: %zu reads collapsed, %zu clusters\n",
                     state_->swc_result->stats.reads_collapsed,
                     state_->swc_result->num_clusters);
    }

    return true;
}

bool AllpairPipeline::SelectRepresentatives() {
    if (!state_->graph || !state_->clustering) {
        state_->last_error = "Graph or clustering not available";
        return false;
    }

    ClusterRepConfig rep_config;
    rep_config.method = ClusterRepConfig::Method::Medoid;
    rep_config.use_confidence_tiebreaker = true;

    ClusterRepSelector selector(rep_config);

    if (state_->swc_result) {
        state_->rep_result = selector.SelectWithConfidence(
            *state_->graph, *state_->clustering, *state_->swc_result);
    } else {
        state_->rep_result = selector.Select(*state_->graph, *state_->clustering);
    }

    if (!state_->rep_result) {
        state_->last_error = "Representative selection failed";
        return false;
    }

    if (config_.verbose) {
        std::fprintf(stderr, "Selected %zu representatives\n",
                     state_->rep_result->representatives.size());
    }

    return true;
}

bool AllpairPipeline::WriteOutput() {
    if (!state_->swc_result) {
        state_->last_error = "No results to write";
        return false;
    }

    output::ClusterWriterConfig writer_config;
    writer_config.format = output::ClusterOutputFormat::TSV;
    writer_config.include_header = true;

    auto writer = output::ClusterWriter::Create(config_.output, writer_config);
    if (!writer) {
        state_->last_error = "Failed to create output writer";
        return false;
    }

    const size_t num_reads = state_->reads.size();
    std::vector<size_t> cluster_sizes;
    if (state_->swc_result) {
        cluster_sizes = state_->swc_result->GetClusterSizes();
    }

    for (size_t i = 0; i < num_reads; ++i) {
        output::ClusterAssignment assignment;
        assignment.read_id = state_->reads[i].id;
        assignment.read_length = state_->reads[i].sequence.size();

        // Find this read's assignment
        bool found = false;
        for (const auto& a : state_->swc_result->assignments) {
            if (a.read_idx == i) {
                assignment.cluster_id = a.cluster_id;
                assignment.confidence = a.confidence;
                assignment.cluster_size = (a.cluster_id < cluster_sizes.size())
                    ? cluster_sizes[a.cluster_id] : 0;
                found = true;
                break;
            }
        }

        if (!found) {
            // Fallback: use Leiden assignment
            if (state_->clustering && i < state_->clustering->labels.size()) {
                assignment.cluster_id = state_->clustering->labels[i];
                assignment.confidence = 0.5f;
            } else {
                assignment.cluster_id = 0;
                assignment.confidence = 0.0f;
            }
            assignment.cluster_size = 0;
        }

        // Check if representative
        assignment.is_representative = false;
        if (state_->rep_result) {
            assignment.is_representative = state_->rep_result->IsRepresentative(
                static_cast<uint32_t>(i));
        }

        if (!writer->Write(assignment)) {
            state_->last_error = "Failed to write assignment: " + writer->LastError();
            return false;
        }
    }

    writer->Close();
    return true;
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
