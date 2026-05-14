// LLmap — Per-cell paralog matrix: core implementation.

#include "per_cell_paralog.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace llmap::singlecell {

void ParalogAccumulator::Add(float probability, float confidence) {
    prob_sum += probability;
    prob_max = std::max(prob_max, probability);
    confidence_sum += confidence;
    ++count;
}

float ParalogAccumulator::GetProbability(CellParalogConfig::AggregationMethod method) const {
    if (count == 0) return 0.0f;

    switch (method) {
        case CellParalogConfig::AggregationMethod::Mean:
            return prob_sum / static_cast<float>(count);
        case CellParalogConfig::AggregationMethod::Max:
            return prob_max;
        case CellParalogConfig::AggregationMethod::Sum:
            return prob_sum;
        case CellParalogConfig::AggregationMethod::Weighted:
            if (confidence_sum > 0.0f) {
                return prob_sum / confidence_sum;
            }
            return prob_sum / static_cast<float>(count);
    }
    return 0.0f;
}

float ParalogAccumulator::GetMeanConfidence() const {
    if (count == 0) return 0.0f;
    return confidence_sum / static_cast<float>(count);
}

void CellParalogMatrix::AddRecord(const AlignmentRecord& record) {
    if (finalized_) {
        return;
    }

    if (!record.cell_barcode.has_value() || record.cell_barcode->empty()) {
        return;
    }

    if (!record.paralog_assignment.has_value()) {
        return;
    }

    const auto& cell = *record.cell_barcode;
    const auto& paralog_call = *record.paralog_assignment;

    float base_confidence = 1.0f;
    if (!record.confidence_scores.empty()) {
        base_confidence = record.confidence_scores.front();
    }

    for (const auto& [paralog_id, prob] : paralog_call.inter_paralog) {
        data_[cell][paralog_id].Add(prob, base_confidence * prob);
    }
}

void CellParalogMatrix::AddRecords(const std::vector<AlignmentRecord>& records) {
    for (const auto& record : records) {
        AddRecord(record);
    }
}

void CellParalogMatrix::Finalize(const CellParalogConfig& config) {
    if (finalized_) return;

    config_ = config;
    entries_.clear();
    cells_.clear();
    paralogs_.clear();

    std::unordered_set<std::string> unique_paralogs;

    for (auto& [cell, paralog_map] : data_) {
        std::uint32_t cell_read_count = 0;
        for (const auto& [paralog_id, acc] : paralog_map) {
            cell_read_count += acc.count;
        }

        if (cell_read_count < config.min_reads_per_cell) {
            continue;
        }

        float row_sum = 0.0f;
        std::vector<CellParalogEntry> cell_entries;

        for (const auto& [paralog_id, acc] : paralog_map) {
            if (acc.count < config.min_reads_per_paralog) {
                continue;
            }

            float prob = acc.GetProbability(config.aggregation);
            if (prob < config.min_probability) {
                continue;
            }

            CellParalogEntry entry;
            entry.cell_barcode = cell;
            entry.paralog_id = paralog_id;
            entry.probability = prob;
            entry.read_count = acc.count;
            entry.mean_confidence = acc.GetMeanConfidence();

            row_sum += prob;
            cell_entries.push_back(std::move(entry));
            unique_paralogs.insert(paralog_id);
        }

        if (config.normalize_rows && row_sum > 0.0f) {
            for (auto& entry : cell_entries) {
                entry.probability /= row_sum;
            }
        }

        for (auto& entry : cell_entries) {
            entries_.push_back(std::move(entry));
        }

        cells_.push_back(cell);
    }

    std::sort(cells_.begin(), cells_.end());

    for (const auto& p : unique_paralogs) {
        paralogs_.push_back(p);
    }
    std::sort(paralogs_.begin(), paralogs_.end());

    stats_.total_reads = 0;
    stats_.reads_with_barcode = 0;
    stats_.reads_with_paralog = 0;
    stats_.reads_with_both = 0;
    for (const auto& [cell, paralog_map] : data_) {
        for (const auto& [paralog_id, acc] : paralog_map) {
            stats_.reads_with_both += acc.count;
        }
    }
    stats_.unique_cells = cells_.size();
    stats_.unique_paralogs = paralogs_.size();
    stats_.matrix_entries = entries_.size();

    std::uint64_t dense_size = stats_.unique_cells * stats_.unique_paralogs;
    if (dense_size > 0) {
        stats_.sparsity = 1.0f - static_cast<float>(stats_.matrix_entries) /
                                 static_cast<float>(dense_size);
    }

    finalized_ = true;
}

std::vector<CellParalogEntry> CellParalogMatrix::GetEntries() const {
    return entries_;
}

std::vector<CellParalogEntry> CellParalogMatrix::GetEntriesForCell(
    std::string_view cell_barcode) const {
    std::vector<CellParalogEntry> result;
    for (const auto& entry : entries_) {
        if (entry.cell_barcode == cell_barcode) {
            result.push_back(entry);
        }
    }
    return result;
}

std::vector<std::string> CellParalogMatrix::GetCells() const {
    return cells_;
}

std::vector<std::string> CellParalogMatrix::GetParalogs() const {
    return paralogs_;
}

CellParalogStats CellParalogMatrix::GetStats() const noexcept {
    return stats_;
}

std::tuple<
    std::vector<std::vector<float>>,
    std::vector<std::string>,
    std::vector<std::string>
> CellParalogMatrix::GetDenseMatrix() const {
    std::unordered_map<std::string, size_t> cell_idx;
    std::unordered_map<std::string, size_t> paralog_idx;

    for (size_t i = 0; i < cells_.size(); ++i) {
        cell_idx[cells_[i]] = i;
    }
    for (size_t i = 0; i < paralogs_.size(); ++i) {
        paralog_idx[paralogs_[i]] = i;
    }

    std::vector<std::vector<float>> matrix(
        cells_.size(),
        std::vector<float>(paralogs_.size(), 0.0f)
    );

    for (const auto& entry : entries_) {
        auto ci = cell_idx.find(entry.cell_barcode);
        auto pi = paralog_idx.find(entry.paralog_id);
        if (ci != cell_idx.end() && pi != paralog_idx.end()) {
            matrix[ci->second][pi->second] = entry.probability;
        }
    }

    return {matrix, cells_, paralogs_};
}

void CellParalogMatrix::Clear() noexcept {
    data_.clear();
    entries_.clear();
    cells_.clear();
    paralogs_.clear();
    stats_ = {};
    finalized_ = false;
}

CellParalogMatrixBuilder::CellParalogMatrixBuilder(CellParalogConfig config)
    : config_(std::move(config)) {}

CellParalogMatrixBuilder& CellParalogMatrixBuilder::AddRecord(const AlignmentRecord& record) {
    matrix_.AddRecord(record);
    return *this;
}

CellParalogMatrixBuilder& CellParalogMatrixBuilder::AddRecords(
    const std::vector<AlignmentRecord>& records) {
    matrix_.AddRecords(records);
    return *this;
}

CellParalogMatrix CellParalogMatrixBuilder::Build() {
    matrix_.Finalize(config_);
    return std::move(matrix_);
}

}  // namespace llmap::singlecell
