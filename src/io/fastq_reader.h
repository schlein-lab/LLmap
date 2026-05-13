// LLmap — FASTQ reader for reading sequence files.
//
// Simple streaming FASTQ reader for the allpair command. Does not depend on
// external libraries; pure C++ implementation sufficient for V1.0.

#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::io {

// A single FASTQ record
struct FastqRecord {
    std::string id;         // Read name (without @ prefix)
    std::string sequence;   // Nucleotide sequence
    std::string quality;    // Phred+33 quality string

    // Convenience accessors
    size_t Length() const { return sequence.size(); }
    bool IsValid() const { return !id.empty() && !sequence.empty(); }
};

// Statistics from reading a FASTQ file
struct FastqReaderStats {
    size_t total_records = 0;
    size_t total_bases = 0;
    size_t min_length = 0;
    size_t max_length = 0;
    float avg_length = 0.0f;
    size_t invalid_records = 0;
    float read_time_ms = 0.0f;
};

// Configuration for the FASTQ reader
struct FastqReaderConfig {
    size_t buffer_size = 1 << 20;       // 1 MB read buffer
    size_t max_records = 0;             // 0 = read all
    size_t min_length = 0;              // Filter records shorter than this
    size_t max_length = 0;              // 0 = no limit
    bool validate_quality = false;       // Check quality string length matches sequence
    bool uppercase_sequence = true;      // Convert sequence to uppercase
};

// Forward declaration for PIMPL
class FastqReaderImpl;

// Streaming FASTQ reader
//
// Supports:
// - Plain FASTQ files (.fastq, .fq)
// - Gzip-compressed FASTQ files (.fastq.gz, .fq.gz) via zlib
class FastqReader {
public:
    // Factory: open a FASTQ file
    // Returns nullptr on failure (file not found, permission denied, etc.)
    static std::unique_ptr<FastqReader> Open(
        const std::filesystem::path& path,
        const FastqReaderConfig& config = {});

    ~FastqReader();

    // Non-copyable, movable
    FastqReader(const FastqReader&) = delete;
    FastqReader& operator=(const FastqReader&) = delete;
    FastqReader(FastqReader&&) noexcept;
    FastqReader& operator=(FastqReader&&) noexcept;

    // Read the next record
    // Returns nullopt at end of file
    std::optional<FastqRecord> Next();

    // Read a batch of records
    // Returns empty vector at end of file
    std::vector<FastqRecord> NextBatch(size_t batch_size);

    // Read all remaining records
    std::vector<FastqRecord> ReadAll();

    // Check if more records are available
    bool HasMore() const;

    // Check if the file is gzip-compressed
    bool IsCompressed() const;

    // Get accumulated statistics
    FastqReaderStats GetStats() const;

    // Get last error message
    std::string LastError() const;

private:
    explicit FastqReader(const std::filesystem::path& path,
                         const FastqReaderConfig& config);
    bool Initialize();

    std::filesystem::path path_;
    FastqReaderConfig config_;
    std::unique_ptr<FastqReaderImpl> impl_;
};

// Convenience: read all records from a FASTQ file
std::vector<FastqRecord> ReadFastq(
    const std::filesystem::path& path,
    size_t max_records = 0);

// Convenience: count records in a FASTQ file without loading sequences
size_t CountFastqRecords(const std::filesystem::path& path);

// Utility: detect if a file is gzip-compressed (by magic bytes)
bool IsGzipFile(const std::filesystem::path& path);

// Utility: validate a FASTQ record
bool ValidateFastqRecord(const FastqRecord& record);

}  // namespace llmap::io
