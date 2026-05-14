// LLmap — Single-cell paralog quantification QC metrics.
//
// Provides quality control metrics and statistical summaries for single-cell
// paralog assignment results. Includes per-cell, per-paralog, and global
// metrics with entropy, confidence distributions, and expression patterns.

#pragma once

#include "singlecell/per_cell_paralog.h"
#include "psv/psv_types.h"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace llmap::singlecell {

// Per-cell QC metrics
struct CellQcMetrics {
    std::string cell_barcode;
    std::uint32_t total_reads{0};
    std::uint32_t assigned_reads{0};
    std::uint32_t confident_reads{0};
    std::uint32_t paralogs_detected{0};
    float assignment_rate{0.0f};
    float mean_confidence{0.0f};
    float median_confidence{0.0f};
    float assignment_entropy{0.0f};
    float dominance_score{0.0f};
};

// Per-paralog QC metrics
struct ParalogQcMetrics {
    std::string paralog_id;
    std::uint32_t total_cells{0};
    std::uint32_t total_reads{0};
    float mean_expression{0.0f};
    float median_expression{0.0f};
    float expression_variance{0.0f};
    float detection_rate{0.0f};
    float coefficient_of_variation{0.0f};
};

// Confidence distribution bins
struct ConfidenceDistribution {
    std::vector<std::uint32_t> histogram;
    float bin_width{0.1f};
    float mean{0.0f};
    float median{0.0f};
    float std_dev{0.0f};
    float skewness{0.0f};
    std::uint32_t below_threshold{0};
    std::uint32_t above_threshold{0};
};

// Global QC summary
struct GlobalQcSummary {
    // Cell statistics
    std::uint64_t total_cells{0};
    std::uint64_t cells_with_assignments{0};
    std::uint64_t cells_passing_qc{0};
    float mean_paralogs_per_cell{0.0f};
    float median_paralogs_per_cell{0.0f};

    // Read statistics
    std::uint64_t total_reads{0};
    std::uint64_t reads_with_barcode{0};
    std::uint64_t reads_assigned{0};
    std::uint64_t reads_confident{0};
    float global_assignment_rate{0.0f};

    // Paralog statistics
    std::uint64_t total_paralogs{0};
    std::uint64_t paralogs_detected{0};
    float mean_detection_rate{0.0f};

    // Confidence metrics
    ConfidenceDistribution confidence_dist;
    float mean_entropy{0.0f};
    float median_entropy{0.0f};

    // PSV statistics (if available)
    std::uint64_t total_psv_observations{0};
    std::uint64_t informative_psv_observations{0};
    float mean_psvs_per_read{0.0f};
};

// QC thresholds configuration
struct QcThresholds {
    float min_assignment_rate{0.1f};
    float min_confidence{0.5f};
    std::uint32_t min_reads_per_cell{10};
    std::uint32_t min_paralogs_per_cell{1};
    float max_entropy{1.5f};
    float min_detection_rate{0.01f};
};

// QC report with all metrics
struct QcReport {
    GlobalQcSummary global;
    std::vector<CellQcMetrics> cells;
    std::vector<ParalogQcMetrics> paralogs;
    QcThresholds thresholds;
    std::string timestamp;
    std::string sample_id;
};

// Compute per-cell QC metrics from matrix
[[nodiscard]] std::vector<CellQcMetrics> ComputeCellQcMetrics(
    const CellParalogMatrix& matrix,
    const QcThresholds& thresholds = {});

// Compute per-paralog QC metrics from matrix
[[nodiscard]] std::vector<ParalogQcMetrics> ComputeParalogQcMetrics(
    const CellParalogMatrix& matrix);

// Compute global QC summary
[[nodiscard]] GlobalQcSummary ComputeGlobalQcSummary(
    const CellParalogMatrix& matrix,
    const std::vector<CellQcMetrics>& cell_metrics,
    const std::vector<ParalogQcMetrics>& paralog_metrics);

// Compute confidence distribution from raw confidence values
[[nodiscard]] ConfidenceDistribution ComputeConfidenceDistribution(
    const std::vector<float>& confidences,
    float threshold = 0.5f,
    float bin_width = 0.1f);

// Compute entropy of a probability distribution
[[nodiscard]] float ComputeEntropy(const std::vector<float>& probs);

// Compute dominance score (max_prob / second_max_prob)
[[nodiscard]] float ComputeDominance(const std::vector<float>& probs);

// Generate full QC report
[[nodiscard]] QcReport GenerateQcReport(
    const CellParalogMatrix& matrix,
    const QcThresholds& thresholds = {},
    std::string_view sample_id = "");

// Export QC report to JSON
bool ExportQcReportJson(
    const QcReport& report,
    const std::filesystem::path& path);

// Export QC report to TSV (multiple files: cells.tsv, paralogs.tsv, summary.tsv)
bool ExportQcReportTsv(
    const QcReport& report,
    const std::filesystem::path& output_dir);

// Export cell QC metrics to TSV
bool ExportCellQcTsv(
    const std::vector<CellQcMetrics>& metrics,
    const std::filesystem::path& path);

// Export paralog QC metrics to TSV
bool ExportParalogQcTsv(
    const std::vector<ParalogQcMetrics>& metrics,
    const std::filesystem::path& path);

// Export global summary to TSV
bool ExportSummaryTsv(
    const GlobalQcSummary& summary,
    const std::filesystem::path& path);

// Get cells passing QC filters
[[nodiscard]] std::vector<std::string> GetCellsPassingQc(
    const std::vector<CellQcMetrics>& metrics,
    const QcThresholds& thresholds);

// Filter matrix to only include QC-passing cells
[[nodiscard]] CellParalogMatrix FilterMatrixByQc(
    const CellParalogMatrix& matrix,
    const std::vector<std::string>& passing_cells);

}  // namespace llmap::singlecell
