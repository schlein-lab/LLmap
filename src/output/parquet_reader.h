// LLmap — Parquet probabilistic reader for round-trip validation.
//
// Reads lossless probabilistic alignments from Parquet or CSV fallback format.
// Used for:
// - Round-trip validation (write→read→compare)
// - External analysis tool integration
// - Lossless invariant verification

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "output/parquet_writer.h"

namespace llmap::output {

// Statistics from reading
struct ParquetReaderStats {
    std::size_t entries_read = 0;
    std::size_t unique_reads = 0;
    std::size_t unique_buckets = 0;
    std::size_t collapsed_entries = 0;
    std::size_t unmapped_entries = 0;
    float read_time_ms = 0.0f;
};

// Reader configuration
struct ParquetReaderConfig {
    float min_probability = 0.0f;    // Filter low-probability entries
    bool skip_unmapped = false;      // Skip entries with bucket_id="*"
};

// Forward declaration
class ParquetReaderImpl;

// Parquet/CSV reader for lossless probabilistic data
class ParquetReader {
public:
    // Factory: open a file for reading
    static std::unique_ptr<ParquetReader> Open(
        const std::filesystem::path& path,
        const ParquetReaderConfig& config = {});

    ~ParquetReader();

    // Non-copyable, movable
    ParquetReader(const ParquetReader&) = delete;
    ParquetReader& operator=(const ParquetReader&) = delete;
    ParquetReader(ParquetReader&&) noexcept;
    ParquetReader& operator=(ParquetReader&&) noexcept;

    // Read all entries into a vector
    std::vector<ProbabilityEntry> ReadAll();

    // Read entries in batches (returns empty when done)
    std::vector<ProbabilityEntry> ReadBatch(std::size_t batch_size);

    // Check if more data is available
    bool HasMore() const;

    // Get accumulated statistics
    ParquetReaderStats GetStats() const;

    // Get last error
    std::string LastError() const;

    // Get the detected format
    ParquetOutputFormat DetectedFormat() const;

private:
    explicit ParquetReader(const std::filesystem::path& path,
                           const ParquetReaderConfig& config);
    bool Initialize();

    std::filesystem::path path_;
    ParquetReaderConfig config_;
    std::unique_ptr<ParquetReaderImpl> impl_;
};

// Convenience: read all entries from a file
std::vector<ProbabilityEntry> ReadParquet(
    const std::filesystem::path& path,
    const ParquetReaderConfig& config = {});

// Parse a single CSV line into a ProbabilityEntry
std::optional<ProbabilityEntry> ParseCSVLine(std::string_view line);

// Group entries by read_id
std::vector<std::vector<ProbabilityEntry>> GroupByReadId(
    const std::vector<ProbabilityEntry>& entries);

// Round-trip validation: compare original entries with re-read entries
struct RoundTripResult {
    bool success = false;
    std::size_t entries_original = 0;
    std::size_t entries_reread = 0;
    std::size_t mismatches = 0;
    std::string error;
};

RoundTripResult ValidateRoundTrip(
    const std::vector<ProbabilityEntry>& original,
    const std::vector<ProbabilityEntry>& reread,
    float probability_tolerance = 1e-5f);

}  // namespace llmap::output
