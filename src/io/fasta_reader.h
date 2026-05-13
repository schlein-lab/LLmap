// LLmap — FASTA reader for reading reference files.
//
// Simple streaming FASTA reader for loading reference sequences.
// Pure C++ implementation sufficient for V1.0.

#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace llmap::io {

// A single FASTA record
struct FastaRecord {
    std::string name;       // Sequence name (without > prefix)
    std::string sequence;   // Nucleotide sequence

    size_t Length() const { return sequence.size(); }
    bool IsValid() const { return !name.empty() && !sequence.empty(); }
};

// Configuration for the FASTA reader
struct FastaReaderConfig {
    bool uppercase_sequence = true;   // Convert sequence to uppercase
    bool skip_N_only = false;         // Skip sequences that are only Ns
    size_t max_records = 0;           // 0 = read all
};

// Forward declaration for PIMPL
class FastaReaderImpl;

// Streaming FASTA reader
class FastaReader {
public:
    explicit FastaReader(const std::filesystem::path& path,
                         const FastaReaderConfig& config = {});
    ~FastaReader();

    // Non-copyable, movable
    FastaReader(const FastaReader&) = delete;
    FastaReader& operator=(const FastaReader&) = delete;
    FastaReader(FastaReader&&) noexcept;
    FastaReader& operator=(FastaReader&&) noexcept;

    // Read the next record
    FastaRecord Next();

    // Check if more records are available
    bool HasMore() const;

    // Get total records read
    size_t RecordsRead() const;

    // Get last error message
    std::string LastError() const;

private:
    std::filesystem::path path_;
    FastaReaderConfig config_;
    std::unique_ptr<FastaReaderImpl> impl_;
};

// Convenience: read all records from a FASTA file
std::vector<FastaRecord> ReadFasta(const std::filesystem::path& path);

}  // namespace llmap::io
