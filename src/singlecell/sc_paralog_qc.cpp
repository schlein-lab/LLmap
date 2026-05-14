// LLmap — Single-cell paralog QC metric computation.

#include "singlecell/sc_paralog_qc.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace llmap::singlecell {

namespace {

float ComputeMedian(std::vector<float> values) {
    if (values.empty()) return 0.0f;
    const auto mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    if (values.size() % 2 == 0) {
        auto mid_val = values[mid];
        std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
        return (values[mid - 1] + mid_val) / 2.0f;
    }
    return values[mid];
}

float ComputeStdDev(const std::vector<float>& values, float mean) {
    if (values.size() < 2) return 0.0f;
    float sum_sq = 0.0f;
    for (auto v : values) {
        float diff = v - mean;
        sum_sq += diff * diff;
    }
    return std::sqrt(sum_sq / static_cast<float>(values.size() - 1));
}

float ComputeSkewness(const std::vector<float>& values, float mean, float std_dev) {
    if (values.size() < 3 || std_dev < 1e-9f) return 0.0f;
    float sum_cubed = 0.0f;
    for (auto v : values) {
        float z = (v - mean) / std_dev;
        sum_cubed += z * z * z;
    }
    auto n = static_cast<float>(values.size());
    return sum_cubed / n;
}

}  // namespace

float ComputeEntropy(const std::vector<float>& probs) {
    if (probs.empty()) return 0.0f;
    float entropy = 0.0f;
    for (auto p : probs) {
        if (p > 1e-9f) {
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

float ComputeDominance(const std::vector<float>& probs) {
    if (probs.size() < 2) return probs.empty() ? 0.0f : 1.0f;
    float max1 = 0.0f, max2 = 0.0f;
    for (auto p : probs) {
        if (p > max1) {
            max2 = max1;
            max1 = p;
        } else if (p > max2) {
            max2 = p;
        }
    }
    return max2 > 1e-9f ? max1 / max2 : (max1 > 0.0f ? 100.0f : 0.0f);
}

ConfidenceDistribution ComputeConfidenceDistribution(
    const std::vector<float>& confidences,
    float threshold,
    float bin_width) {

    ConfidenceDistribution dist;
    dist.bin_width = bin_width;

    if (confidences.empty()) return dist;

    const auto num_bins = static_cast<size_t>(std::ceil(1.0f / bin_width));
    dist.histogram.resize(num_bins, 0);

    float sum = 0.0f;
    for (auto c : confidences) {
        sum += c;
        auto bin = static_cast<size_t>(c / bin_width);
        if (bin >= num_bins) bin = num_bins - 1;
        dist.histogram[bin]++;

        if (c < threshold) {
            dist.below_threshold++;
        } else {
            dist.above_threshold++;
        }
    }

    dist.mean = sum / static_cast<float>(confidences.size());
    dist.median = ComputeMedian(std::vector<float>(confidences));
    dist.std_dev = ComputeStdDev(confidences, dist.mean);
    dist.skewness = ComputeSkewness(confidences, dist.mean, dist.std_dev);

    return dist;
}

std::vector<CellQcMetrics> ComputeCellQcMetrics(
    const CellParalogMatrix& matrix,
    const QcThresholds& thresholds) {

    std::vector<CellQcMetrics> result;
    const auto entries = matrix.GetEntries();

    std::unordered_map<std::string, std::vector<CellParalogEntry>> by_cell;
    for (const auto& e : entries) {
        by_cell[e.cell_barcode].push_back(e);
    }

    result.reserve(by_cell.size());
    for (const auto& [cell, cell_entries] : by_cell) {
        CellQcMetrics m;
        m.cell_barcode = cell;
        m.paralogs_detected = static_cast<std::uint32_t>(cell_entries.size());

        std::vector<float> confidences;
        std::vector<float> probs;
        confidences.reserve(cell_entries.size());
        probs.reserve(cell_entries.size());

        for (const auto& e : cell_entries) {
            m.total_reads += e.read_count;
            if (e.probability > 0.0f) {
                m.assigned_reads += e.read_count;
                if (e.mean_confidence >= thresholds.min_confidence) {
                    m.confident_reads += e.read_count;
                }
            }
            confidences.push_back(e.mean_confidence);
            probs.push_back(e.probability);
        }

        m.assignment_rate = m.total_reads > 0
            ? static_cast<float>(m.assigned_reads) / static_cast<float>(m.total_reads)
            : 0.0f;

        if (!confidences.empty()) {
            m.mean_confidence = std::accumulate(
                confidences.begin(), confidences.end(), 0.0f) /
                static_cast<float>(confidences.size());
            m.median_confidence = ComputeMedian(confidences);
        }

        m.assignment_entropy = ComputeEntropy(probs);
        m.dominance_score = ComputeDominance(probs);

        result.push_back(std::move(m));
    }

    return result;
}

std::vector<ParalogQcMetrics> ComputeParalogQcMetrics(
    const CellParalogMatrix& matrix) {

    std::vector<ParalogQcMetrics> result;
    const auto entries = matrix.GetEntries();
    const auto cells = matrix.GetCells();
    const auto total_cells = cells.size();

    std::unordered_map<std::string, std::vector<CellParalogEntry>> by_paralog;
    for (const auto& e : entries) {
        by_paralog[e.paralog_id].push_back(e);
    }

    result.reserve(by_paralog.size());
    for (const auto& [paralog, paralog_entries] : by_paralog) {
        ParalogQcMetrics m;
        m.paralog_id = paralog;
        m.total_cells = static_cast<std::uint32_t>(paralog_entries.size());

        std::vector<float> expressions;
        expressions.reserve(paralog_entries.size());

        for (const auto& e : paralog_entries) {
            m.total_reads += e.read_count;
            expressions.push_back(e.probability);
        }

        if (!expressions.empty()) {
            m.mean_expression = std::accumulate(
                expressions.begin(), expressions.end(), 0.0f) /
                static_cast<float>(expressions.size());
            m.median_expression = ComputeMedian(expressions);
            m.expression_variance = ComputeStdDev(expressions, m.mean_expression);
            m.expression_variance *= m.expression_variance;
            m.coefficient_of_variation = m.mean_expression > 1e-9f
                ? std::sqrt(m.expression_variance) / m.mean_expression
                : 0.0f;
        }

        m.detection_rate = total_cells > 0
            ? static_cast<float>(m.total_cells) / static_cast<float>(total_cells)
            : 0.0f;

        result.push_back(std::move(m));
    }

    return result;
}

GlobalQcSummary ComputeGlobalQcSummary(
    const CellParalogMatrix& matrix,
    const std::vector<CellQcMetrics>& cell_metrics,
    const std::vector<ParalogQcMetrics>& paralog_metrics) {

    GlobalQcSummary summary;
    const auto stats = matrix.GetStats();

    summary.total_cells = cell_metrics.size();
    summary.total_paralogs = paralog_metrics.size();
    summary.total_reads = stats.total_reads;
    summary.reads_with_barcode = stats.reads_with_barcode;

    std::vector<float> paralogs_per_cell;
    std::vector<float> all_confidences;
    std::vector<float> all_entropies;
    paralogs_per_cell.reserve(cell_metrics.size());
    all_entropies.reserve(cell_metrics.size());

    for (const auto& cm : cell_metrics) {
        paralogs_per_cell.push_back(static_cast<float>(cm.paralogs_detected));
        all_entropies.push_back(cm.assignment_entropy);
        summary.reads_assigned += cm.assigned_reads;
        summary.reads_confident += cm.confident_reads;
        if (cm.assigned_reads > 0) {
            summary.cells_with_assignments++;
        }
        all_confidences.push_back(cm.mean_confidence);
    }

    if (!paralogs_per_cell.empty()) {
        summary.mean_paralogs_per_cell = std::accumulate(
            paralogs_per_cell.begin(), paralogs_per_cell.end(), 0.0f) /
            static_cast<float>(paralogs_per_cell.size());
        summary.median_paralogs_per_cell = ComputeMedian(paralogs_per_cell);
    }

    if (!all_entropies.empty()) {
        summary.mean_entropy = std::accumulate(
            all_entropies.begin(), all_entropies.end(), 0.0f) /
            static_cast<float>(all_entropies.size());
        summary.median_entropy = ComputeMedian(all_entropies);
    }

    summary.global_assignment_rate = summary.total_reads > 0
        ? static_cast<float>(summary.reads_assigned) /
          static_cast<float>(summary.total_reads)
        : 0.0f;

    std::unordered_set<std::string> detected;
    for (const auto& pm : paralog_metrics) {
        if (pm.total_cells > 0) {
            detected.insert(pm.paralog_id);
            summary.mean_detection_rate += pm.detection_rate;
        }
    }
    summary.paralogs_detected = detected.size();
    if (!paralog_metrics.empty()) {
        summary.mean_detection_rate /= static_cast<float>(paralog_metrics.size());
    }

    summary.confidence_dist = ComputeConfidenceDistribution(all_confidences);

    return summary;
}

std::vector<std::string> GetCellsPassingQc(
    const std::vector<CellQcMetrics>& metrics,
    const QcThresholds& thresholds) {

    std::vector<std::string> passing;
    passing.reserve(metrics.size());

    for (const auto& m : metrics) {
        if (m.total_reads >= thresholds.min_reads_per_cell &&
            m.paralogs_detected >= thresholds.min_paralogs_per_cell &&
            m.assignment_rate >= thresholds.min_assignment_rate &&
            m.mean_confidence >= thresholds.min_confidence &&
            m.assignment_entropy <= thresholds.max_entropy) {
            passing.push_back(m.cell_barcode);
        }
    }

    return passing;
}

CellParalogMatrix FilterMatrixByQc(
    const CellParalogMatrix& matrix,
    const std::vector<std::string>& passing_cells) {

    std::unordered_set<std::string> passing_set(
        passing_cells.begin(), passing_cells.end());

    CellParalogMatrix filtered;
    for (const auto& entry : matrix.GetEntries()) {
        if (passing_set.count(entry.cell_barcode) > 0) {
            AlignmentRecord rec;
            rec.cell_barcode = entry.cell_barcode;
            ParalogCall pc;
            pc.inter_paralog.emplace_back(entry.paralog_id, entry.probability);
            rec.paralog_assignment = pc;
            rec.confidence_scores.push_back(entry.mean_confidence);
            for (std::uint32_t i = 0; i < entry.read_count; ++i) {
                filtered.AddRecord(rec);
            }
        }
    }

    return filtered;
}

QcReport GenerateQcReport(
    const CellParalogMatrix& matrix,
    const QcThresholds& thresholds,
    std::string_view sample_id) {

    QcReport report;
    report.thresholds = thresholds;
    report.sample_id = std::string(sample_id);

    auto now = std::time(nullptr);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    report.timestamp = buf;

    report.cells = ComputeCellQcMetrics(matrix, thresholds);
    report.paralogs = ComputeParalogQcMetrics(matrix);
    report.global = ComputeGlobalQcSummary(matrix, report.cells, report.paralogs);

    auto passing = GetCellsPassingQc(report.cells, thresholds);
    report.global.cells_passing_qc = passing.size();

    return report;
}

}  // namespace llmap::singlecell
