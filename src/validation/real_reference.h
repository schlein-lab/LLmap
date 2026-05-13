// LLmap — Real reference validation for Phase 5.3.
//
// Validates LLmap against real hg38 chr14 IGH locus data.
// Supports minimap2 baseline comparison and SLURM job submission.

#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/alignment_record.h"
#include "validation/killswitch.h"

namespace llmap::validation {

// Configuration for real reference validation
struct RealReferenceConfig {
    // Reference paths
    std::filesystem::path reference_fasta;  // hg38 chr14 IGH locus FASTA
    std::filesystem::path reference_index;  // Pre-built LLmap index (optional)

    // Read paths (either real or synthetic)
    std::filesystem::path reads_fastq;      // Input reads
    std::filesystem::path ground_truth_bed; // Ground truth positions (optional)

    // Minimap2 baseline (optional)
    std::filesystem::path minimap2_bam;     // Pre-aligned minimap2 results

    // Validation parameters
    KillSwitchThresholds thresholds;
    PositionTolerance position_tolerance;

    // Pipeline parameters (mirror EndToEndConfig)
    uint8_t minimizer_k = 15;
    uint8_t minimizer_w = 10;
    float min_identity = 0.70f;

    // Output paths
    std::filesystem::path output_dir;
    bool save_alignments = true;
    bool save_metrics = true;

    // GPU settings
    bool use_gpu = false;
    int gpu_device = 0;

    // Validate paths exist
    bool Validate() const;
    std::string ValidationErrors() const;
};

// Ground truth from BED file (real coordinates)
struct RealGroundTruth {
    std::string read_id;
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    std::string strand;
    std::string gene;        // Gene name (if available)
    std::string annotation;  // Additional annotation (if available)
};

// Result of real reference validation
struct RealReferenceResult {
    // Input stats
    size_t n_reference_seqs = 0;
    size_t reference_total_bp = 0;
    size_t n_reads = 0;

    // LLmap results
    size_t llmap_mapped = 0;
    size_t llmap_unmapped = 0;
    float llmap_alignment_rate = 0.0f;

    // Minimap2 baseline (if provided)
    bool has_minimap2_baseline = false;
    size_t minimap2_mapped = 0;
    BaselineComparison baseline_comparison;

    // Position accuracy (if ground truth provided)
    bool has_ground_truth = false;
    ValidationStats validation;

    // Kill-switch verdicts
    bool recall_pass = false;        // ≥99.5% of minimap2
    bool lossless_pass = false;      // No silent drops
    bool overall_pass = false;

    std::string verdict_reason;

    // Timing
    float index_load_time_ms = 0.0f;
    float alignment_time_ms = 0.0f;
    float validation_time_ms = 0.0f;
    float total_time_ms = 0.0f;

    std::string Summary() const;
};

// Parse ground truth from BED file
std::vector<RealGroundTruth> ParseGroundTruthBed(
    const std::filesystem::path& bed_path);

// Parse minimap2 BAM to get mapping status per read
std::unordered_map<std::string, bool> ParseMinimap2Bam(
    const std::filesystem::path& bam_path);

// Run real reference validation
RealReferenceResult RunRealReferenceValidation(const RealReferenceConfig& config);

// SLURM job management
struct SlurmJobConfig {
    std::string job_name = "llmap_validation";
    std::string partition = "gpu";
    int n_gpus = 1;
    int n_cpus = 8;
    std::string memory = "32G";
    std::string time_limit = "4:00:00";
    std::filesystem::path work_dir;
    std::filesystem::path output_log;
    std::filesystem::path error_log;
};

struct SlurmJobStatus {
    std::string job_id;
    std::string state;      // PENDING, RUNNING, COMPLETED, FAILED
    int exit_code = -1;
    bool is_complete = false;
    std::string output_file;
};

// Generate SLURM submission script
std::string GenerateSlurmScript(
    const RealReferenceConfig& validation_config,
    const SlurmJobConfig& slurm_config);

// Submit job and return job ID
std::optional<std::string> SubmitSlurmJob(
    const std::filesystem::path& script_path);

// Check job status
SlurmJobStatus CheckSlurmJob(const std::string& job_id);

// Parse results from completed job output
std::optional<RealReferenceResult> ParseSlurmOutput(
    const std::filesystem::path& output_path);

}  // namespace llmap::validation
