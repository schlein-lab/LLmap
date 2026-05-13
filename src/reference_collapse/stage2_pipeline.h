// LLmap — Stage2Pipeline: reference-based WaveCollapse orchestration.
//
// This pipeline implements Stage 2 of the LLmap algorithm:
//   1. Load reference index (BucketPyramid + embeddings)
//   2. Initialize WaveState from Stage 1 representatives
//   3. Run EM iterations at each level (L0 → L1 → L2)
//   4. Check collapse/convergence after each iteration
//   5. Refine (expand buckets) when level converges
//   6. Propagate final positions to cluster members
//
// Input: Stage 1 results (cluster reps, cluster assignments, similarity graph)
// Output: Genomic positions for all reads

#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "core/wave_state.h"

namespace llmap::self_interference {
struct ClusterRepResult;
struct ClusteringResult;
class SimilarityGraph;
}  // namespace llmap::self_interference

namespace llmap {

class ReferenceIndex;
class WaveState;

// Configuration for Stage 2 pipeline
struct Stage2Config {
    // Reference index path
    std::filesystem::path index_path;

    // EM iteration parameters
    std::uint32_t max_iterations_per_level = 50;
    float gamma = 0.3f;                    // EM damping factor
    float tau_collapse = 0.99f;            // Collapse threshold
    float convergence_delta = 1e-5f;       // Max prob change for convergence
    bool apply_smoothing = true;           // K(b,b') spatial coupling

    // Refinement parameters
    float expansion_threshold = 0.01f;     // Expand buckets with P >= this
    std::uint32_t max_candidates = 50;     // Max candidates after expansion
    float refine_trigger_rate = 0.8f;      // Refine when X% collapsed

    // Member propagation parameters
    float propagation_confidence = 0.8f;   // Base confidence for members
    float similarity_weight = 0.5f;        // Similarity→confidence weight

    // Performance
    std::uint32_t num_threads = 0;         // 0 = auto
    bool verbose = false;
};

// Progress callback
using Stage2ProgressCallback = std::function<void(
    const std::string& stage,
    std::size_t current,
    std::size_t total,
    const std::string& message)>;

// Per-level statistics
struct LevelStats {
    std::uint32_t iterations = 0;
    std::uint32_t reads_processed = 0;
    std::uint32_t reads_collapsed = 0;
    std::uint32_t reads_refined = 0;
    float avg_entropy = 0.0f;
    float time_ms = 0.0f;
};

// Statistics from running Stage 2
struct Stage2Stats {
    // Input
    std::size_t num_representatives = 0;
    std::size_t num_clusters = 0;
    std::size_t total_reads = 0;

    // Reference index
    std::size_t ref_targets = 0;
    std::size_t ref_l2_buckets = 0;

    // Per-level stats
    LevelStats l0_stats;
    LevelStats l1_stats;
    LevelStats l2_stats;

    // Member propagation
    std::uint32_t members_propagated = 0;
    float propagation_time_ms = 0.0f;

    // Overall
    std::uint32_t total_iterations = 0;
    float total_time_ms = 0.0f;
    float collapse_rate = 0.0f;
};

// Per-read alignment result
struct ReadAlignment {
    std::uint32_t read_idx;
    std::string target_name;           // Reference sequence name
    std::uint64_t position;            // Genomic position (0-based)
    float confidence;                  // Alignment confidence [0, 1]
    bool is_representative;            // Was this a Stage 1 rep?
    std::uint32_t cluster_id;          // Original Stage 1 cluster
    std::uint32_t l2_bucket_id;        // Final L2 bucket
};

// Result of Stage 2 pipeline
struct Stage2Result {
    std::vector<ReadAlignment> alignments;
    Stage2Stats stats;

    bool success = false;
    std::string error_message;

    // Convenience accessors
    [[nodiscard]] std::size_t NumAlignments() const { return alignments.size(); }

    // Get alignments for a specific read
    [[nodiscard]] std::vector<ReadAlignment> GetAlignmentsForRead(std::uint32_t read_idx) const;

    // Get alignments above confidence threshold
    [[nodiscard]] std::vector<ReadAlignment> GetHighConfidenceAlignments(float min_confidence) const;
};

// The Stage 2 reference WaveCollapse pipeline
class Stage2Pipeline {
public:
    Stage2Pipeline();
    explicit Stage2Pipeline(const Stage2Config& config);
    ~Stage2Pipeline();

    Stage2Pipeline(Stage2Pipeline&&) noexcept;
    Stage2Pipeline& operator=(Stage2Pipeline&&) noexcept;
    Stage2Pipeline(const Stage2Pipeline&) = delete;
    Stage2Pipeline& operator=(const Stage2Pipeline&) = delete;

    // Run the full Stage 2 pipeline
    std::unique_ptr<Stage2Result> Run(
        const self_interference::ClusteringResult& clustering,
        const self_interference::ClusterRepResult& rep_result,
        const self_interference::SimilarityGraph& graph);

    // Run with progress callback
    std::unique_ptr<Stage2Result> Run(
        const self_interference::ClusteringResult& clustering,
        const self_interference::ClusterRepResult& rep_result,
        const self_interference::SimilarityGraph& graph,
        Stage2ProgressCallback callback);

    // Run stages individually (for testing/debugging)
    bool LoadReferenceIndex();
    bool InitializeWaveState(const self_interference::ClusterRepResult& rep_result);
    bool RunEMLevel(WaveLevel level);
    bool RefineLevel(WaveLevel from_level);
    bool PropagateToMembers(
        const self_interference::ClusteringResult& clustering,
        const self_interference::ClusterRepResult& rep_result,
        const self_interference::SimilarityGraph& graph);
    bool ExtractAlignments();

    // Set reference index directly (alternative to loading from path)
    void SetReferenceIndex(std::unique_ptr<ReferenceIndex> index);

    // Get/set configuration
    [[nodiscard]] const Stage2Config& Config() const { return config_; }
    void SetConfig(const Stage2Config& config) { config_ = config; }

    // Get last error
    [[nodiscard]] std::string LastError() const;

    // Access internal state (for testing)
    [[nodiscard]] const WaveState* GetWaveState() const;

private:
    Stage2Config config_;

    struct InternalState;
    std::unique_ptr<InternalState> state_;

    void ReportProgress(const std::string& stage, std::size_t current,
                        std::size_t total, const std::string& message);
    Stage2ProgressCallback progress_callback_;
};

// Convenience: run Stage 2 with default config
std::unique_ptr<Stage2Result> RunStage2Pipeline(
    const std::filesystem::path& index_path,
    const self_interference::ClusteringResult& clustering,
    const self_interference::ClusterRepResult& rep_result,
    const self_interference::SimilarityGraph& graph);

}  // namespace llmap
