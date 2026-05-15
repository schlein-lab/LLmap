// LLmap — cmd_align internal header.
// Shared types and declarations for split cmd_align modules.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "classical/classical_pipeline.h"
#include "claude_agent/pipeline_agent.h"
#include "core/alignment_record.h"
#include "core/thread_pool.h"
#include "io/fastq_reader.h"
#include "output/bam_writer.h"
#include "output/parquet_writer.h"
#include "psv/psv_catalog.h"

namespace llmap::cli::align_internal {

// Read type presets for optimized alignment parameters
enum class Preset {
    None,      // No preset, use defaults
    MapHifi,   // PacBio HiFi: high accuracy, k=19, w=19, tight thresholds
    MapOnt,    // Oxford Nanopore: higher error, k=15, w=10, relaxed identity
    MapPb,     // Legacy PacBio CLR: high error, same as map-ont
    Sr         // Short reads (Illumina): k=21, w=11, high accuracy
};

// Convert preset string to enum
Preset ParsePreset(std::string_view name);

// Get human-readable preset name
const char* PresetName(Preset preset);

// Apply preset defaults to AlignArgs (call after parsing, before alignment)
void ApplyPreset(Preset preset, struct AlignArgs& args);

struct AlignArgs {
    std::string reads;
    std::string reference;
    std::string output;
    std::string parquet_output;
    std::string index;          // Pre-built .llmi index file (optional)
    Preset preset = Preset::None;  // Read type preset (-x)
    bool use_bam = false;
    bool use_sam = true;
    int kmer_size = 15;
    int window_size = 10;
    int min_chain = 30;
    float min_identity = 0.80f;
    int threads = 1;
    int max_chains = 10;
    bool verbose = false;
    bool help = false;

    bool enable_llm = false;
    std::string llm_api_key;
    float llm_threshold = 0.50f;
    std::string llm_work_dir;

    // Layer 3 active checkpoint consult mode: "off" | "auto" | "required".
    // Independent from the legacy --llm diagnostic toggle above. The
    // dispatcher fires inside chain-DP / extension at ambiguous chains,
    // unknown regions, paralog disambiguation, SD expansion, novel insertion.
    std::string llm_mode = "off";
    std::string llm_cache_dir;     // optional override; default ~/.cache/llmap

    // PSV-based paralog disambiguation
    std::string psv_catalog;       // Path to PSV catalog (BED/VCF)
    float psv_weight = 0.5f;       // Weight for PSV vs probabilistic assignment
    float psv_min_posterior = 0.9f; // Min posterior for confident call
    bool psv_only = false;          // Use only PSV (skip prob assignment)

    // Classical-only mode (skip probabilistic framework)
    bool classical_only = false;   // Pure seed-chain-extend, no Wave structures

    // Region-aware mapping. When set, load the .regann annotation file and
    // pass it to the chain DP so per-anchor parameters can be overridden
    // per region (paralog, low-complexity, centromere, ...).
    std::string region_annot;
};

void PrintAlignUsage();
bool ParseAlignArgs(int argc, char** argv, AlignArgs& args);

std::string GetApiKey(const AlignArgs& args);

std::optional<claude_agent::PipelineAgent> CreatePipelineAgent(
    const AlignArgs& args, bool verbose);

void RunLlmDiagnostics(
    claude_agent::PipelineAgent& agent,
    const AlignArgs& args,
    const std::vector<classical::ReadAlignmentResult>& results,
    std::size_t n_mapped, std::size_t n_unmapped,
    float avg_identity);

// PSV integration
std::optional<psv::PsvCatalog> LoadPsvCatalog(const AlignArgs& args);

void ApplyPsvAssignments(
    const psv::PsvCatalog& catalog,
    const AlignArgs& args,
    std::vector<AlignmentRecord>& records,
    const std::vector<std::string>& read_sequences,
    bool verbose);

// Reporting (cmd_align_report.cpp)
void PrintAlignmentSummary(
    const AlignArgs& args,
    const classical::ClassicalPipelineStats& stats,
    std::size_t total_reads,
    std::size_t n_mapped,
    std::size_t n_unmapped,
    float align_time_ms,
    float total_time_ms);

bool ShouldRunLlmDiagnostics(const AlignArgs& args, float mapping_rate);
bool IsLlmEnabledButSkipped(const AlignArgs& args, float mapping_rate);

// Batch processing result (cmd_align_run.cpp)
struct BatchAlignResult {
    classical::ClassicalPipelineStats agg_stats;
    std::size_t total_reads = 0;
    std::size_t n_mapped = 0;
    std::size_t n_unmapped = 0;
    double identity_sum_weighted = 0.0;
    bool error = false;
};

// Convert classical alignment to AlignmentRecord
AlignmentRecord ConvertClassicalAlignment(
    const classical::ClassicalAlignment& aln,
    std::uint32_t read_len);

// Run alignment batches, streaming through the pipeline
BatchAlignResult RunAlignmentBatches(
    classical::ClassicalPipeline& pipeline,
    io::FastqReader& read_reader,
    output::BamWriter& bam_writer,
    output::ParquetWriter* parquet_writer,
    core::ThreadPool* thread_pool,
    const std::optional<psv::PsvCatalog>& psv_catalog,
    const AlignArgs& args);

// Finalize aggregated stats (compute avg_identity)
void FinalizeAlignStats(BatchAlignResult& result);

}  // namespace llmap::cli::align_internal
