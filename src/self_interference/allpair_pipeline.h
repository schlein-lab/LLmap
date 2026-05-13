// LLmap — AllpairPipeline: Stage 1 Self-Interference orchestration.
//
// This pipeline implements the `llmap allpair` command: reads FASTQ, embeds
// sequences, builds similarity graph, clusters with Leiden, refines with
// Self-WaveCollapse, selects representatives, and outputs cluster assignments.
//
// Corresponds to LLmap_SPEC.md §2.5 Stage 1 and §6 CLI `llmap allpair`.

#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace llmap::self_interference {

// Forward declarations
struct ClusterRepResult;
struct ClusteringResult;
struct SelfWaveCollapseResult;

// Configuration for the allpair pipeline
struct AllpairConfig {
    // Input/output paths
    std::filesystem::path reads_a;              // Primary input FASTQ
    std::filesystem::path reads_b;              // Optional second input (for cross-comparison)
    std::filesystem::path output;               // Output file (TSV or Parquet)
    std::filesystem::path model_path;           // Foundation model ONNX path

    // Embedding parameters
    size_t embedding_batch_size = 1024;
    size_t embedding_dim = 256;
    bool use_gpu_embedder = false;
    int gpu_device_id = 0;

    // FAISS ANN parameters
    size_t faiss_k = 50;                        // k-NN neighbors
    size_t faiss_nlist = 1024;                  // IVF cells
    size_t faiss_nprobe = 32;                   // Cells to probe
    bool use_gpu_faiss = false;

    // Leiden clustering parameters
    float leiden_resolution = 1.0f;
    size_t leiden_max_iterations = 100;
    size_t min_cluster_size = 2;

    // Self-WaveCollapse parameters
    size_t swc_max_iterations = 20;
    float swc_convergence_threshold = 1e-5f;
    float swc_collapse_threshold = 0.95f;
    float swc_gamma = 0.3f;

    // Representative selection
    bool select_representatives = true;

    // Performance
    size_t num_threads = 0;                     // 0 = auto-detect
    bool verbose = false;

    // Read filtering
    size_t min_read_length = 0;
    size_t max_read_length = 0;
    size_t max_reads = 0;                       // 0 = no limit
};

// Progress callback for long-running operations
using ProgressCallback = std::function<void(
    const std::string& stage,
    size_t current,
    size_t total,
    const std::string& message)>;

// Statistics from running the pipeline
struct AllpairStats {
    // Input stats
    size_t input_reads = 0;
    size_t total_bases = 0;
    float avg_read_length = 0.0f;

    // Embedding stats
    float embedding_time_ms = 0.0f;
    float embeddings_per_second = 0.0f;

    // ANN search stats
    float ann_build_time_ms = 0.0f;
    float ann_search_time_ms = 0.0f;

    // Graph stats
    size_t graph_nodes = 0;
    size_t graph_edges = 0;
    float graph_build_time_ms = 0.0f;

    // Clustering stats
    size_t num_clusters = 0;
    float modularity = 0.0f;
    float leiden_time_ms = 0.0f;

    // Self-WaveCollapse stats
    size_t swc_iterations = 0;
    size_t reads_collapsed = 0;
    float swc_time_ms = 0.0f;

    // Representative selection stats
    size_t num_representatives = 0;
    float rep_selection_time_ms = 0.0f;

    // Overall
    float total_time_ms = 0.0f;
    float reads_per_second = 0.0f;
};

// Result of the allpair pipeline
struct AllpairResult {
    // Per-read cluster assignments
    std::vector<uint32_t> cluster_ids;

    // Per-read confidence scores
    std::vector<float> confidences;

    // Which reads are representatives
    std::vector<bool> is_representative;

    // Read IDs (from FASTQ)
    std::vector<std::string> read_ids;

    // Cluster sizes (indexed by cluster_id)
    std::vector<size_t> cluster_sizes;

    // Statistics
    AllpairStats stats;

    // Check if successful
    bool success = false;
    std::string error_message;

    // Number of reads
    size_t NumReads() const { return cluster_ids.size(); }

    // Number of clusters
    size_t NumClusters() const { return cluster_sizes.size(); }

    // Get representative read indices
    std::vector<uint32_t> GetRepresentatives() const;
};

// The Stage 1 self-interference pipeline
class AllpairPipeline {
public:
    AllpairPipeline();
    explicit AllpairPipeline(const AllpairConfig& config);

    ~AllpairPipeline();

    // Move-only
    AllpairPipeline(AllpairPipeline&&) noexcept;
    AllpairPipeline& operator=(AllpairPipeline&&) noexcept;
    AllpairPipeline(const AllpairPipeline&) = delete;
    AllpairPipeline& operator=(const AllpairPipeline&) = delete;

    // Run the full pipeline
    std::unique_ptr<AllpairResult> Run();

    // Run with progress callback
    std::unique_ptr<AllpairResult> Run(ProgressCallback callback);

    // Run individual stages (for debugging/testing)
    bool LoadReads();
    bool ComputeEmbeddings();
    bool BuildSimilarityGraph();
    bool RunClustering();
    bool RunSelfWaveCollapse();
    bool SelectRepresentatives();
    bool WriteOutput();

    // Get/set configuration
    const AllpairConfig& GetConfig() const { return config_; }
    void SetConfig(const AllpairConfig& config) { config_ = config; }

    // Get last error
    std::string LastError() const;

private:
    AllpairConfig config_;

    struct InternalState;
    std::unique_ptr<InternalState> state_;

    void ReportProgress(const std::string& stage, size_t current, size_t total,
                        const std::string& message);
    ProgressCallback progress_callback_;
};

// Convenience: run allpair pipeline with default config
std::unique_ptr<AllpairResult> RunAllpairPipeline(
    const std::filesystem::path& reads,
    const std::filesystem::path& output,
    const std::filesystem::path& model_path = {});

}  // namespace llmap::self_interference
