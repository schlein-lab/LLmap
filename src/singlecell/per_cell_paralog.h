// LLmap — Per-cell paralog probability matrix.
//
// Aggregates paralog assignment probabilities per cell barcode to produce
// a cell × paralog matrix suitable for single-cell analysis workflows.
// Supports export to CSV/TSV and AnnData .h5ad format.

#pragma once

#include "core/alignment_record.h"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::singlecell {

// A single entry in the sparse cell × paralog matrix
struct CellParalogEntry {
    std::string cell_barcode;
    std::string paralog_id;
    float probability{0.0f};
    std::uint32_t read_count{0};
    float mean_confidence{0.0f};
};

// Statistics for matrix construction
struct CellParalogStats {
    std::uint64_t total_reads{0};
    std::uint64_t reads_with_barcode{0};
    std::uint64_t reads_with_paralog{0};
    std::uint64_t reads_with_both{0};
    std::uint64_t unique_cells{0};
    std::uint64_t unique_paralogs{0};
    std::uint64_t matrix_entries{0};
    float sparsity{0.0f};
};

// Configuration for matrix construction
struct CellParalogConfig {
    // Minimum probability to include in matrix
    float min_probability{0.01f};

    // Minimum reads per cell to include
    std::uint32_t min_reads_per_cell{1};

    // Minimum reads supporting a paralog assignment
    std::uint32_t min_reads_per_paralog{1};

    // Aggregation method for multiple reads
    enum class AggregationMethod : std::uint8_t {
        Mean,       // Average probability across reads
        Max,        // Maximum probability
        Sum,        // Sum of probabilities (can exceed 1.0)
        Weighted,   // Weighted by confidence scores
    };
    AggregationMethod aggregation{AggregationMethod::Mean};

    // Whether to normalize rows to sum to 1.0
    bool normalize_rows{false};
};

// Per-cell paralog accumulator (internal, per-cell-paralog)
struct ParalogAccumulator {
    float prob_sum{0.0f};
    float prob_max{0.0f};
    float confidence_sum{0.0f};
    std::uint32_t count{0};

    void Add(float probability, float confidence);
    [[nodiscard]] float GetProbability(CellParalogConfig::AggregationMethod method) const;
    [[nodiscard]] float GetMeanConfidence() const;
};

// Cell × Paralog probability matrix
// Stored as sparse: unordered_map<cell_barcode, map<paralog_id, accumulator>>
class CellParalogMatrix {
public:
    // Add a single alignment record to the matrix
    void AddRecord(const AlignmentRecord& record);

    // Add multiple records
    void AddRecords(const std::vector<AlignmentRecord>& records);

    // Finalize and apply filtering
    void Finalize(const CellParalogConfig& config);

    // Get all entries as a flat vector
    [[nodiscard]] std::vector<CellParalogEntry> GetEntries() const;

    // Get entries for a specific cell
    [[nodiscard]] std::vector<CellParalogEntry> GetEntriesForCell(
        std::string_view cell_barcode) const;

    // Get all unique cell barcodes
    [[nodiscard]] std::vector<std::string> GetCells() const;

    // Get all unique paralog IDs
    [[nodiscard]] std::vector<std::string> GetParalogs() const;

    // Get matrix statistics
    [[nodiscard]] CellParalogStats GetStats() const noexcept;

    // Get dense matrix (cells × paralogs) with ordered indices
    // Returns: (matrix, cell_names, paralog_names)
    [[nodiscard]] std::tuple<
        std::vector<std::vector<float>>,
        std::vector<std::string>,
        std::vector<std::string>
    > GetDenseMatrix() const;

    // Clear all data
    void Clear() noexcept;

    // Check if finalized
    [[nodiscard]] bool IsFinalized() const noexcept { return finalized_; }

private:
    // cell_barcode -> (paralog_id -> accumulator)
    std::unordered_map<std::string, std::unordered_map<std::string, ParalogAccumulator>> data_;

    // Cached after Finalize()
    std::vector<CellParalogEntry> entries_;
    std::vector<std::string> cells_;
    std::vector<std::string> paralogs_;
    CellParalogStats stats_;
    CellParalogConfig config_;
    bool finalized_{false};
};

// Export matrix to CSV file
// Format: cell_barcode,paralog_id,probability,read_count,mean_confidence
bool ExportToCSV(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path);

// Export matrix to TSV file
bool ExportToTSV(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path);

// Export to dense matrix CSV (cells as rows, paralogs as columns)
// First row: header with paralog IDs
// First column: cell barcodes
bool ExportToDenseCSV(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path);

// Export to AnnData .h5ad format (if HDF5 available, else returns false)
// X: sparse CSR matrix of probabilities
// obs: cell metadata (barcode, read_count)
// var: paralog metadata (paralog_id)
bool ExportToH5AD(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path);

// Import matrix from CSV file
std::optional<CellParalogMatrix> ImportFromCSV(
    const std::filesystem::path& path);

// Builder pattern for creating matrix from records
class CellParalogMatrixBuilder {
public:
    explicit CellParalogMatrixBuilder(CellParalogConfig config = {});

    CellParalogMatrixBuilder& AddRecord(const AlignmentRecord& record);
    CellParalogMatrixBuilder& AddRecords(const std::vector<AlignmentRecord>& records);

    [[nodiscard]] CellParalogMatrix Build();

private:
    CellParalogMatrix matrix_;
    CellParalogConfig config_;
};

}  // namespace llmap::singlecell
