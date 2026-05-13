// LLmap — Parquet probabilistic output writer.
//
// Outputs full lossless probabilistic alignments in Apache Parquet format.
// Schema: (read_id, bucket_id, probability, confidence, level, collapsed_at_iter)
//
// Uses Apache Arrow/Parquet C++ if available (LLMAP_HAS_ARROW), otherwise
// writes a simple CSV fallback format that can be loaded into Arrow.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "core/alignment_record.h"

namespace llmap::output {

// Parquet output format
enum class ParquetOutputFormat {
    Parquet,   // Native Parquet (requires Arrow)
    CSV,       // CSV fallback (always available)
};

// A single probability entry for one read × bucket pair
struct ProbabilityEntry {
    std::string read_id;
    std::string bucket_id;       // Target ID or bucket name
    float probability;           // P(bucket|read)
    float confidence;            // Optional confidence score
    std::uint8_t level;          // Bucket level (0=L0, 1=L1, etc.)
    std::uint32_t iteration;     // Iteration at which this probability was recorded
    bool is_collapsed;           // Whether this entry represents collapsed state
};

// Statistics from writing
struct ParquetWriterStats {
    std::size_t records_written = 0;
    std::size_t entries_written = 0;    // Total probability entries
    std::size_t mapped_records = 0;
    std::size_t tentative_records = 0;
    std::size_t unmapped_records = 0;
    float write_time_ms = 0.0f;
    std::size_t bytes_written = 0;
};

// Configuration for the Parquet writer
struct ParquetWriterConfig {
    ParquetOutputFormat format = ParquetOutputFormat::Parquet;

    // Compression
    bool compress = true;               // Use compression (Snappy default)
    std::string compression = "snappy"; // snappy, gzip, zstd, none

    // Output options
    bool include_alternatives = true;   // Include alternative alignments
    bool include_tentative = true;      // Include tentative targets
    bool include_unmapped = true;       // Include unmapped reads
    float min_probability = 0.0f;       // Filter entries below threshold

    // Metadata
    std::string sample_name = "sample";
    std::string program_version = "1.0.0";

    // Batching for memory efficiency
    std::size_t row_group_size = 100000;  // Rows per row group
};

// Forward declaration
class ParquetWriterImpl;

// Parquet writer for lossless probabilistic output
class ParquetWriter {
public:
    // Factory: create a writer for the specified output path
    static std::unique_ptr<ParquetWriter> Create(
        const std::filesystem::path& path,
        const ParquetWriterConfig& config = {});

    ~ParquetWriter();

    // Non-copyable, movable
    ParquetWriter(const ParquetWriter&) = delete;
    ParquetWriter& operator=(const ParquetWriter&) = delete;
    ParquetWriter(ParquetWriter&&) noexcept;
    ParquetWriter& operator=(ParquetWriter&&) noexcept;

    // Write a single alignment record (expands to multiple probability entries)
    bool Write(const AlignmentRecord& record);

    // Write explicit probability entries
    bool WriteEntries(std::span<const ProbabilityEntry> entries);

    // Write a batch of alignment records
    bool WriteBatch(std::span<const AlignmentRecord> records);

    // Finalize and close the file
    bool Close();

    // Get accumulated statistics
    ParquetWriterStats GetStats() const;

    // Get last error
    std::string LastError() const;

    // Check if Arrow/Parquet is available
    static bool ArrowAvailable();

private:
    explicit ParquetWriter(const std::filesystem::path& path,
                           const ParquetWriterConfig& config);
    bool Initialize();

    std::filesystem::path path_;
    ParquetWriterConfig config_;
    std::unique_ptr<ParquetWriterImpl> impl_;
};

// Convenience: write all records to a file in one call
bool WriteParquet(
    const std::filesystem::path& path,
    std::span<const AlignmentRecord> records,
    const ParquetWriterConfig& config = {});

// Convert AlignmentRecord to probability entries
std::vector<ProbabilityEntry> RecordToEntries(
    const AlignmentRecord& record,
    float min_probability = 0.0f);

}  // namespace llmap::output
