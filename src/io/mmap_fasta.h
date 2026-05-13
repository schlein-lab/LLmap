// LLmap — Memory-mapped FASTA reader for large references.
//
// Uses mmap to avoid loading entire reference genomes into heap memory.
// Provides lazy loading: sequences are only accessed when needed.
// Returns string_view for zero-copy access to sequence data.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::io {

// Forward declaration
class MmapFastaImpl;

// Lazy sequence reference — does not own the data
struct MmapSequence {
    std::string_view name;      // Sequence name (without > prefix)
    size_t offset;              // Offset in mmap'd file where sequence starts
    size_t length;              // Total sequence length

    bool IsValid() const { return !name.empty() && length > 0; }
};

// Statistics about the memory-mapped file
struct MmapStats {
    size_t file_size = 0;         // Total file size in bytes
    size_t num_sequences = 0;     // Number of sequences in file
    size_t total_bases = 0;       // Total nucleotide count
    size_t resident_pages = 0;    // Pages currently in RAM (estimate)
    double resident_fraction = 0; // Fraction of file in RAM
};

// Configuration for the mmap FASTA reader
struct MmapFastaConfig {
    bool build_index_on_open = true;  // Scan file to build sequence index
    bool prefault_pages = false;      // Prefault pages into RAM on open
    bool lock_pages = false;          // mlock pages (requires privileges)
    size_t read_ahead_bytes = 0;      // madvise read-ahead hint (0 = none)
};

// Memory-mapped FASTA reader
class MmapFastaReader {
public:
    explicit MmapFastaReader(const std::filesystem::path& path,
                             const MmapFastaConfig& config = {});
    ~MmapFastaReader();

    // Non-copyable, movable
    MmapFastaReader(const MmapFastaReader&) = delete;
    MmapFastaReader& operator=(const MmapFastaReader&) = delete;
    MmapFastaReader(MmapFastaReader&&) noexcept;
    MmapFastaReader& operator=(MmapFastaReader&&) noexcept;

    // Check if file was opened successfully
    bool IsValid() const;

    // Get last error message
    std::string LastError() const;

    // Number of sequences in file
    size_t NumSequences() const;

    // Get sequence metadata by index
    std::optional<MmapSequence> GetSequence(size_t index) const;

    // Get sequence metadata by name
    std::optional<MmapSequence> GetSequence(std::string_view name) const;

    // Get sequence data as string_view (zero-copy, may contain newlines)
    // The view is valid as long as the reader is alive
    std::string_view GetSequenceRaw(size_t index) const;
    std::string_view GetSequenceRaw(std::string_view name) const;

    // Get contiguous sequence data (copies, removes newlines)
    // Use this when you need the actual nucleotide sequence
    std::string GetSequenceData(size_t index) const;
    std::string GetSequenceData(std::string_view name) const;

    // Get a subsequence without loading the entire sequence
    std::string GetSubsequence(size_t index, size_t start, size_t len) const;
    std::string GetSubsequence(std::string_view name, size_t start, size_t len) const;

    // Get all sequence names
    std::vector<std::string_view> SequenceNames() const;

    // Memory statistics
    MmapStats GetStats() const;

    // Advise kernel about access pattern
    void AdviseSequential();    // Will read file sequentially
    void AdviseRandom();        // Random access pattern
    void AdviseWillNeed(size_t index);  // Prefetch a sequence
    void AdviseDontNeed(size_t index);  // Release a sequence from RAM

    // File path
    const std::filesystem::path& Path() const { return path_; }

private:
    std::filesystem::path path_;
    MmapFastaConfig config_;
    std::unique_ptr<MmapFastaImpl> impl_;
};

// Convenience: check if a file is a valid FASTA (by extension and header)
bool IsFastaFile(const std::filesystem::path& path);

}  // namespace llmap::io
