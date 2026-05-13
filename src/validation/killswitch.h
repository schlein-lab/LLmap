// LLmap — Kill-switch validation framework.
//
// Phase 5 validation: the project kill-switch tests.
// This module provides infrastructure for verifying:
//   1. Lossless guarantee: count(input) == count(output)
//   2. Position accuracy: mapped positions are within tolerance of ground truth
//   3. Recall/precision metrics: mapped vs unmapped vs minimap2 baseline
//   4. Paralog accuracy: dup vs canonical discrimination

#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "core/alignment_record.h"
#include "synthetic/igh_locus_generator.h"

namespace llmap::validation {

// Position tolerance for considering a mapping correct
struct PositionTolerance {
    int64_t base_tolerance = 100;   // Absolute bp tolerance
    float relative_tolerance = 0.01f;  // Fraction of read length
};

// Ground truth from synthetic data
struct GroundTruth {
    std::string read_id;
    std::string origin;           // "CAN" or "DUP"
    uint64_t true_start = 0;
    uint64_t true_end = 0;
    std::string source_gene;
    size_t n_psvs_covered = 0;

    static std::optional<GroundTruth> ParseFromReadId(std::string_view read_id);
};

// Validation result for a single read
struct ReadValidation {
    std::string read_id;

    // Status
    bool has_ground_truth = false;
    bool has_alignment_record = false;
    bool is_lossless_consistent = false;

    // Position accuracy (if mapped)
    bool position_correct = false;
    int64_t position_error = 0;  // Actual - expected start

    // Origin accuracy (canonical vs duplicate)
    bool origin_known = false;
    bool origin_correct = false;

    // Status category
    AlignmentStatus mapped_status = AlignmentStatus::Unmapped;
};

// Aggregated validation statistics
struct ValidationStats {
    // Lossless guarantee
    size_t input_reads = 0;
    size_t output_records = 0;
    size_t missing_records = 0;
    bool lossless_pass = false;

    // Status breakdown
    size_t mapped = 0;
    size_t tentative = 0;
    size_t unmapped = 0;

    // Position accuracy (of mapped reads with ground truth)
    size_t position_evaluated = 0;
    size_t position_correct = 0;
    double position_accuracy = 0.0;

    // Origin accuracy (synthetic only: canonical vs duplicate)
    size_t origin_evaluated = 0;
    size_t origin_correct = 0;
    double origin_accuracy = 0.0;

    // Consistency
    size_t lossless_consistent = 0;
    size_t lossless_inconsistent = 0;

    // Kill-switch verdict
    bool kill_switch_pass = false;
    std::string kill_reason;

    void Compute();
    std::string Summary() const;
};

// Kill-switch thresholds (from LLmap_SPEC.md Phase 5)
struct KillSwitchThresholds {
    // Required: recall >= 99.5% of uniquely-mappable (minimap2 baseline)
    double min_recall = 0.995;

    // Required: paralog accuracy > baseline + 10pp
    double min_paralog_gain = 0.10;

    // Required: dup-fraction within ±2% of ground truth
    double max_dup_fraction_error = 0.02;

    // Required: no silent read drops
    bool require_lossless = true;

    // Required: convergence rate >= 80%
    double min_convergence_rate = 0.80;
};

// The main validation harness
class KillSwitchValidator {
public:
    KillSwitchValidator();
    explicit KillSwitchValidator(KillSwitchThresholds thresholds);

    // Load ground truth from synthetic dataset
    void LoadGroundTruth(const synthetic::GeneratedDataset& dataset);
    void LoadGroundTruth(const std::vector<GroundTruth>& truths);

    // Add alignment records for validation
    void AddRecords(const std::vector<AlignmentRecord>& records);
    void AddRecord(const AlignmentRecord& record);

    // Run validation and get results
    ValidationStats Validate() const;
    std::vector<ReadValidation> ValidateDetailed() const;

    // Kill-switch verdict
    bool PassKillSwitch() const;
    std::string KillSwitchVerdict() const;

    // Configuration
    void SetPositionTolerance(PositionTolerance tol) { position_tolerance_ = tol; }
    void SetThresholds(KillSwitchThresholds t) { thresholds_ = t; }

    // Accessors
    size_t NumGroundTruth() const { return ground_truth_.size(); }
    size_t NumRecords() const { return records_.size(); }

private:
    KillSwitchThresholds thresholds_;
    PositionTolerance position_tolerance_;

    std::unordered_map<std::string, GroundTruth> ground_truth_;
    std::unordered_map<std::string, AlignmentRecord> records_;

    ReadValidation ValidateRead(
        const std::string& read_id,
        const GroundTruth* truth,
        const AlignmentRecord* record) const;

    bool CheckPositionCorrect(
        const AlignmentRecord& record,
        const GroundTruth& truth) const;
};

// Convenience: validate end-to-end on synthetic data
ValidationStats ValidateSyntheticRun(
    const synthetic::GeneratedDataset& dataset,
    const std::vector<AlignmentRecord>& records,
    const KillSwitchThresholds& thresholds = {});

// Compare against minimap2 baseline
struct BaselineComparison {
    size_t llmap_mapped = 0;
    size_t minimap2_mapped = 0;
    size_t both_mapped = 0;
    size_t llmap_only = 0;
    size_t minimap2_only = 0;
    size_t neither_mapped = 0;

    double llmap_recall = 0.0;       // llmap_mapped / total
    double minimap2_recall = 0.0;    // minimap2_mapped / total
    double recall_ratio = 0.0;       // llmap / minimap2
};

BaselineComparison CompareToBaseline(
    const std::vector<AlignmentRecord>& llmap_records,
    const std::unordered_map<std::string, bool>& minimap2_mapped);

}  // namespace llmap::validation

// Forward declare classical types
namespace llmap::classical {
struct ClassicalAlignment;
struct ReadAlignmentResult;
struct ClassicalPipelineConfig;
}  // namespace llmap::classical

namespace llmap::validation {

// Convert ClassicalAlignment to AlignmentRecord for validation
AlignmentRecord ConvertToAlignmentRecord(
    const classical::ClassicalAlignment& aln,
    uint32_t read_len);

// Convert batch of ReadAlignmentResult to AlignmentRecords
std::vector<AlignmentRecord> ConvertResults(
    const std::vector<classical::ReadAlignmentResult>& results,
    const std::vector<uint32_t>& read_lengths);

// End-to-end synthetic validation configuration
struct EndToEndConfig {
    // Synthetic data parameters
    uint64_t seed = 42;
    uint64_t locus_length = 20000;
    uint32_t n_psvs = 10;
    uint32_t read_length = 5000;
    float dup_fraction = 0.0f;
    uint32_t n_reads = 50;

    // Classical pipeline tuning
    uint8_t minimizer_k = 15;
    uint8_t minimizer_w = 10;
    float min_identity = 0.70f;

    // Validation thresholds
    double min_position_accuracy = 0.90;  // 90% of mapped reads correctly positioned
    double position_tolerance_bp = 200;   // bp tolerance for position check
    bool require_lossless = true;

    // Create preset configurations
    static EndToEndConfig Minimal();   // Fast test: 10 reads, small locus
    static EndToEndConfig Standard();  // Standard: 50 reads, medium locus
    static EndToEndConfig Stress();    // Stress test: 200 reads, large locus
};

// Result of end-to-end validation
struct EndToEndResult {
    // Synthetic data stats
    size_t n_reads_generated = 0;
    size_t n_canonical = 0;
    size_t n_duplicate = 0;
    float actual_dup_fraction = 0.0f;

    // Alignment stats
    size_t n_aligned = 0;
    size_t n_unmapped = 0;
    float alignment_rate = 0.0f;
    float avg_identity = 0.0f;

    // Validation stats (from KillSwitchValidator)
    ValidationStats validation;

    // Position accuracy details
    double position_accuracy = 0.0;
    int64_t mean_position_error = 0;
    int64_t max_position_error = 0;

    // Overall verdict
    bool passed = false;
    std::string verdict_reason;

    // Timing
    float generation_time_ms = 0.0f;
    float alignment_time_ms = 0.0f;
    float validation_time_ms = 0.0f;
    float total_time_ms = 0.0f;

    std::string Summary() const;
};

// Run complete end-to-end synthetic validation
// 1. Generate synthetic data with known ground truth
// 2. Build reference index from synthetic locus
// 3. Align reads with classical pipeline
// 4. Validate results through kill-switch framework
EndToEndResult RunEndToEndValidation(const EndToEndConfig& config);

// Convenience wrappers
EndToEndResult RunMinimalValidation(uint64_t seed = 42);
EndToEndResult RunStandardValidation(uint64_t seed = 42);

}  // namespace llmap::validation
