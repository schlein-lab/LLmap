// LLmap — Single-cell paralog QC report export.

#include "singlecell/sc_paralog_qc.h"

#include <fstream>
#include <iomanip>
#include <sstream>

namespace llmap::singlecell {

namespace {

std::string HistogramToJson(const std::vector<std::uint32_t>& histogram) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < histogram.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << histogram[i];
    }
    oss << "]";
    return oss.str();
}

std::string EscapeJson(std::string_view s) {
    std::ostringstream oss;
    for (char c : s) {
        switch (c) {
            case '"': oss << "\\\""; break;
            case '\\': oss << "\\\\"; break;
            case '\n': oss << "\\n"; break;
            case '\r': oss << "\\r"; break;
            case '\t': oss << "\\t"; break;
            default: oss << c;
        }
    }
    return oss.str();
}

}  // namespace

bool ExportQcReportJson(
    const QcReport& report,
    const std::filesystem::path& path) {

    std::ofstream out(path);
    if (!out) return false;

    out << std::fixed << std::setprecision(6);
    out << "{\n";
    out << "  \"sample_id\": \"" << EscapeJson(report.sample_id) << "\",\n";
    out << "  \"timestamp\": \"" << report.timestamp << "\",\n";

    out << "  \"thresholds\": {\n";
    out << "    \"min_assignment_rate\": " << report.thresholds.min_assignment_rate << ",\n";
    out << "    \"min_confidence\": " << report.thresholds.min_confidence << ",\n";
    out << "    \"min_reads_per_cell\": " << report.thresholds.min_reads_per_cell << ",\n";
    out << "    \"min_paralogs_per_cell\": " << report.thresholds.min_paralogs_per_cell << ",\n";
    out << "    \"max_entropy\": " << report.thresholds.max_entropy << ",\n";
    out << "    \"min_detection_rate\": " << report.thresholds.min_detection_rate << "\n";
    out << "  },\n";

    const auto& g = report.global;
    out << "  \"global\": {\n";
    out << "    \"total_cells\": " << g.total_cells << ",\n";
    out << "    \"cells_with_assignments\": " << g.cells_with_assignments << ",\n";
    out << "    \"cells_passing_qc\": " << g.cells_passing_qc << ",\n";
    out << "    \"mean_paralogs_per_cell\": " << g.mean_paralogs_per_cell << ",\n";
    out << "    \"median_paralogs_per_cell\": " << g.median_paralogs_per_cell << ",\n";
    out << "    \"total_reads\": " << g.total_reads << ",\n";
    out << "    \"reads_with_barcode\": " << g.reads_with_barcode << ",\n";
    out << "    \"reads_assigned\": " << g.reads_assigned << ",\n";
    out << "    \"reads_confident\": " << g.reads_confident << ",\n";
    out << "    \"global_assignment_rate\": " << g.global_assignment_rate << ",\n";
    out << "    \"total_paralogs\": " << g.total_paralogs << ",\n";
    out << "    \"paralogs_detected\": " << g.paralogs_detected << ",\n";
    out << "    \"mean_detection_rate\": " << g.mean_detection_rate << ",\n";
    out << "    \"mean_entropy\": " << g.mean_entropy << ",\n";
    out << "    \"median_entropy\": " << g.median_entropy << ",\n";
    out << "    \"confidence_distribution\": {\n";
    out << "      \"histogram\": " << HistogramToJson(g.confidence_dist.histogram) << ",\n";
    out << "      \"bin_width\": " << g.confidence_dist.bin_width << ",\n";
    out << "      \"mean\": " << g.confidence_dist.mean << ",\n";
    out << "      \"median\": " << g.confidence_dist.median << ",\n";
    out << "      \"std_dev\": " << g.confidence_dist.std_dev << ",\n";
    out << "      \"skewness\": " << g.confidence_dist.skewness << ",\n";
    out << "      \"below_threshold\": " << g.confidence_dist.below_threshold << ",\n";
    out << "      \"above_threshold\": " << g.confidence_dist.above_threshold << "\n";
    out << "    }\n";
    out << "  },\n";

    out << "  \"cells\": [\n";
    for (size_t i = 0; i < report.cells.size(); ++i) {
        const auto& c = report.cells[i];
        out << "    {\n";
        out << "      \"cell_barcode\": \"" << EscapeJson(c.cell_barcode) << "\",\n";
        out << "      \"total_reads\": " << c.total_reads << ",\n";
        out << "      \"assigned_reads\": " << c.assigned_reads << ",\n";
        out << "      \"confident_reads\": " << c.confident_reads << ",\n";
        out << "      \"paralogs_detected\": " << c.paralogs_detected << ",\n";
        out << "      \"assignment_rate\": " << c.assignment_rate << ",\n";
        out << "      \"mean_confidence\": " << c.mean_confidence << ",\n";
        out << "      \"median_confidence\": " << c.median_confidence << ",\n";
        out << "      \"assignment_entropy\": " << c.assignment_entropy << ",\n";
        out << "      \"dominance_score\": " << c.dominance_score << "\n";
        out << "    }" << (i + 1 < report.cells.size() ? "," : "") << "\n";
    }
    out << "  ],\n";

    out << "  \"paralogs\": [\n";
    for (size_t i = 0; i < report.paralogs.size(); ++i) {
        const auto& p = report.paralogs[i];
        out << "    {\n";
        out << "      \"paralog_id\": \"" << EscapeJson(p.paralog_id) << "\",\n";
        out << "      \"total_cells\": " << p.total_cells << ",\n";
        out << "      \"total_reads\": " << p.total_reads << ",\n";
        out << "      \"mean_expression\": " << p.mean_expression << ",\n";
        out << "      \"median_expression\": " << p.median_expression << ",\n";
        out << "      \"expression_variance\": " << p.expression_variance << ",\n";
        out << "      \"detection_rate\": " << p.detection_rate << ",\n";
        out << "      \"coefficient_of_variation\": " << p.coefficient_of_variation << "\n";
        out << "    }" << (i + 1 < report.paralogs.size() ? "," : "") << "\n";
    }
    out << "  ]\n";

    out << "}\n";

    return out.good();
}

bool ExportCellQcTsv(
    const std::vector<CellQcMetrics>& metrics,
    const std::filesystem::path& path) {

    std::ofstream out(path);
    if (!out) return false;

    out << std::fixed << std::setprecision(6);
    out << "cell_barcode\ttotal_reads\tassigned_reads\tconfident_reads\t"
        << "paralogs_detected\tassignment_rate\tmean_confidence\t"
        << "median_confidence\tassignment_entropy\tdominance_score\n";

    for (const auto& m : metrics) {
        out << m.cell_barcode << '\t'
            << m.total_reads << '\t'
            << m.assigned_reads << '\t'
            << m.confident_reads << '\t'
            << m.paralogs_detected << '\t'
            << m.assignment_rate << '\t'
            << m.mean_confidence << '\t'
            << m.median_confidence << '\t'
            << m.assignment_entropy << '\t'
            << m.dominance_score << '\n';
    }

    return out.good();
}

bool ExportParalogQcTsv(
    const std::vector<ParalogQcMetrics>& metrics,
    const std::filesystem::path& path) {

    std::ofstream out(path);
    if (!out) return false;

    out << std::fixed << std::setprecision(6);
    out << "paralog_id\ttotal_cells\ttotal_reads\tmean_expression\t"
        << "median_expression\texpression_variance\tdetection_rate\t"
        << "coefficient_of_variation\n";

    for (const auto& m : metrics) {
        out << m.paralog_id << '\t'
            << m.total_cells << '\t'
            << m.total_reads << '\t'
            << m.mean_expression << '\t'
            << m.median_expression << '\t'
            << m.expression_variance << '\t'
            << m.detection_rate << '\t'
            << m.coefficient_of_variation << '\n';
    }

    return out.good();
}

bool ExportSummaryTsv(
    const GlobalQcSummary& summary,
    const std::filesystem::path& path) {

    std::ofstream out(path);
    if (!out) return false;

    out << std::fixed << std::setprecision(6);
    out << "metric\tvalue\n";
    out << "total_cells\t" << summary.total_cells << '\n';
    out << "cells_with_assignments\t" << summary.cells_with_assignments << '\n';
    out << "cells_passing_qc\t" << summary.cells_passing_qc << '\n';
    out << "mean_paralogs_per_cell\t" << summary.mean_paralogs_per_cell << '\n';
    out << "median_paralogs_per_cell\t" << summary.median_paralogs_per_cell << '\n';
    out << "total_reads\t" << summary.total_reads << '\n';
    out << "reads_with_barcode\t" << summary.reads_with_barcode << '\n';
    out << "reads_assigned\t" << summary.reads_assigned << '\n';
    out << "reads_confident\t" << summary.reads_confident << '\n';
    out << "global_assignment_rate\t" << summary.global_assignment_rate << '\n';
    out << "total_paralogs\t" << summary.total_paralogs << '\n';
    out << "paralogs_detected\t" << summary.paralogs_detected << '\n';
    out << "mean_detection_rate\t" << summary.mean_detection_rate << '\n';
    out << "mean_entropy\t" << summary.mean_entropy << '\n';
    out << "median_entropy\t" << summary.median_entropy << '\n';
    out << "confidence_mean\t" << summary.confidence_dist.mean << '\n';
    out << "confidence_median\t" << summary.confidence_dist.median << '\n';
    out << "confidence_std_dev\t" << summary.confidence_dist.std_dev << '\n';
    out << "confidence_skewness\t" << summary.confidence_dist.skewness << '\n';
    out << "confidence_below_threshold\t" << summary.confidence_dist.below_threshold << '\n';
    out << "confidence_above_threshold\t" << summary.confidence_dist.above_threshold << '\n';

    return out.good();
}

bool ExportQcReportTsv(
    const QcReport& report,
    const std::filesystem::path& output_dir) {

    std::filesystem::create_directories(output_dir);

    bool ok = true;
    ok &= ExportCellQcTsv(report.cells, output_dir / "cells_qc.tsv");
    ok &= ExportParalogQcTsv(report.paralogs, output_dir / "paralogs_qc.tsv");
    ok &= ExportSummaryTsv(report.global, output_dir / "summary_qc.tsv");

    return ok;
}

}  // namespace llmap::singlecell
