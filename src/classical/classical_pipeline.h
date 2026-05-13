// LLmap — ClassicalPipeline: seed-chain-extend alignment orchestration.
//
// This pipeline implements the classical alignment path (fallback):
//   1. Extract minimizers from query
//   2. Look up minimizer hits in reference index
//   3. Chain colinear anchors
//   4. Extend chains with WFA2 gap-affine alignment
//   5. Report aligned positions with CIGAR
//
// Used when the neural WaveCollapse path doesn't collapse (uncertainty remains),
// or as a validation baseline for benchmarking.

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "classical/chain.h"
#include "classical/minimizer_index.h"
#include "classical/wfa2_aligner.h"

namespace llmap::classical {

// Configuration for the classical pipeline
struct ClassicalPipelineConfig {
    // Seeding
    MinimizerConfig minimizer_config;

    // Chaining
    ChainConfig chain_config;

    // Extension alignment
    WFA2Config extension_config;

    // Filtering
    float min_identity = 0.70f;        // Minimum alignment identity
    int32_t min_aligned_bases = 50;    // Minimum aligned bases
    uint32_t max_alignments = 5;       // Max alignments to report per read

    // Performance
    uint32_t max_chains_to_extend = 10;   // Max chains to try extending
    bool report_secondary = true;          // Include secondary alignments
};

// A single alignment from the classical pipeline
struct ClassicalAlignment {
    std::string query_name;               // Query read name
    std::string ref_name;                 // Reference sequence name
    uint32_t ref_id = 0;                  // Reference sequence index

    int32_t ref_start = 0;                // Start position in reference (0-based)
    int32_t ref_end = 0;                  // End position in reference
    int32_t query_start = 0;              // Start position in query
    int32_t query_end = 0;                // End position in query

    bool is_forward = true;               // Alignment strand
    bool is_primary = true;               // Primary or secondary alignment

    std::vector<CigarElement> cigar;      // CIGAR operations
    int32_t score = 0;                    // Alignment score
    float identity = 0.0f;                // Fraction of matches
    uint32_t mapq = 0;                    // Mapping quality

    // Chain metadata
    uint32_t chain_anchors = 0;           // Number of anchors in chain
    int32_t chain_score = 0;              // Chain score before extension

    std::string CigarString() const;
    uint32_t AlignedBases() const;
};

// Result of aligning a single read
struct ReadAlignmentResult {
    std::string query_name;
    std::vector<ClassicalAlignment> alignments;

    // Statistics
    size_t num_hits = 0;              // Minimizer hits found
    size_t num_chains = 0;            // Chains extracted
    size_t chains_extended = 0;       // Chains with successful extension
    float seeding_time_us = 0.0f;
    float chaining_time_us = 0.0f;
    float extension_time_us = 0.0f;

    bool HasAlignment() const { return !alignments.empty(); }
    const ClassicalAlignment* PrimaryAlignment() const;
};

// Statistics for batch alignment
struct ClassicalPipelineStats {
    size_t total_reads = 0;
    size_t reads_aligned = 0;
    size_t reads_unmapped = 0;

    size_t total_hits = 0;
    size_t total_chains = 0;
    size_t total_extensions = 0;

    float avg_identity = 0.0f;
    float total_time_ms = 0.0f;
    float seeding_time_ms = 0.0f;
    float chaining_time_ms = 0.0f;
    float extension_time_ms = 0.0f;
};

// The classical seed-chain-extend pipeline
class ClassicalPipeline {
public:
    ClassicalPipeline();
    explicit ClassicalPipeline(const ClassicalPipelineConfig& config);
    ~ClassicalPipeline();

    // Non-copyable, movable
    ClassicalPipeline(ClassicalPipeline&&) noexcept;
    ClassicalPipeline& operator=(ClassicalPipeline&&) noexcept;
    ClassicalPipeline(const ClassicalPipeline&) = delete;
    ClassicalPipeline& operator=(const ClassicalPipeline&) = delete;

    // Set the reference index (required before alignment)
    void SetIndex(std::unique_ptr<MinimizerIndex> index);
    void SetIndex(const MinimizerIndex* index);  // Non-owning

    // Check if index is set
    bool HasIndex() const;

    // Align a single read
    ReadAlignmentResult AlignRead(
        std::string_view query_name,
        std::string_view query_seq) const;

    // Align multiple reads
    std::vector<ReadAlignmentResult> AlignReads(
        std::span<const std::string> query_names,
        std::span<const std::string> query_seqs) const;

    // Get statistics from last batch
    const ClassicalPipelineStats& Stats() const { return stats_; }

    // Configuration access
    const ClassicalPipelineConfig& Config() const { return config_; }

private:
    ClassicalPipelineConfig config_;

    std::unique_ptr<MinimizerIndex> owned_index_;
    const MinimizerIndex* index_ = nullptr;

    WFA2Aligner aligner_;
    mutable ClassicalPipelineStats stats_;

    // Internal methods
    std::optional<ClassicalAlignment> ExtendChain(
        std::string_view query_seq,
        const Chain& chain,
        const std::vector<Anchor>& anchors) const;

    std::vector<std::string> GetRefSequences() const;
};

// Convenience: build index and align in one call
std::vector<ReadAlignmentResult> AlignWithClassicalPath(
    std::span<const std::string> ref_names,
    std::span<const std::string> ref_seqs,
    std::span<const std::string> query_names,
    std::span<const std::string> query_seqs,
    const ClassicalPipelineConfig& config = {});

}  // namespace llmap::classical
