#pragma once

#include "core/arena.h"

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::classical {

// Configuration for minimizer extraction and indexing
struct MinimizerConfig {
    // k-mer size (default matches minimap2 PacBio preset: k=19)
    uint8_t k = 19;

    // Window size for minimizer selection (default matches minimap2: w=19)
    uint8_t w = 19;

    // Use canonical k-mers (forward or reverse-complement, whichever is smaller)
    bool canonical = true;

    // Hash function seed for reproducibility
    uint64_t seed = 0x517cc1b727220a95ULL;

    // Index construction
    size_t max_occ = 500;  // Skip minimizers occurring > max_occ times

    // Memory settings
    size_t bucket_bits = 14;  // Hash table size = 2^bucket_bits
};

// A single minimizer occurrence
struct Minimizer {
    uint64_t hash;       // Hash of the k-mer
    uint32_t pos;        // Position in sequence (0-based)
    bool is_reverse;     // True if reverse-complement was used (canonical mode)
};

// A minimizer hit in the index
struct MinimizerHit {
    uint32_t ref_id;     // Reference sequence ID
    uint32_t ref_pos;    // Position in reference
    uint32_t query_pos;  // Position in query
    bool same_strand;    // True if query and ref are same strand
};

// Statistics from minimizer operations
struct MinimizerStats {
    size_t total_kmers = 0;
    size_t unique_minimizers = 0;
    size_t suppressed_high_occ = 0;
    size_t total_hits = 0;
    float avg_minimizer_spacing = 0.0f;
};

// Reference sequence in the index
struct IndexedSequence {
    std::string name;
    uint32_t length;
};

// MinimizerIndex: k-mer based minimizer index for sequence alignment
//
// This implements the minimizer scheme from minimap2:
// - Extract (k, w)-minimizers using open-syncmer/minimizer hybrid
// - Hash k-mers using invertible hash function
// - Index minimizers by hash for O(1) lookup
// - Support canonical k-mers (strand-agnostic matching)
//
// The index is used for computing L(r|b) in the WaveCollapse update:
// sequences matching many minimizers get higher likelihood.
class MinimizerIndex {
public:
    MinimizerIndex();
    explicit MinimizerIndex(const MinimizerConfig& config);
    ~MinimizerIndex();

    // Move-only
    MinimizerIndex(MinimizerIndex&&) noexcept;
    MinimizerIndex& operator=(MinimizerIndex&&) noexcept;
    MinimizerIndex(const MinimizerIndex&) = delete;
    MinimizerIndex& operator=(const MinimizerIndex&) = delete;

    // Builder pattern for index construction
    class Builder {
    public:
        explicit Builder(const MinimizerConfig& config = {});
        ~Builder();

        Builder(Builder&&) noexcept;
        Builder& operator=(Builder&&) noexcept;

        // Add a reference sequence to the index
        Builder& AddSequence(std::string_view name, std::string_view sequence);

        // Add multiple sequences
        Builder& AddSequences(
            std::span<const std::string> names,
            std::span<const std::string> sequences);

        // Build the final index
        std::unique_ptr<MinimizerIndex> Build();

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };

    // Query the index with a read sequence
    //
    // Returns hits sorted by (ref_id, ref_pos) for chaining
    std::vector<MinimizerHit> Query(std::string_view sequence) const;

    // Query with a maximum hit limit
    std::vector<MinimizerHit> Query(
        std::string_view sequence,
        size_t max_hits) const;

    // Get occurrence count for a minimizer hash
    size_t GetOccurrenceCount(uint64_t hash) const;

    // Get indexed sequence info
    const std::vector<IndexedSequence>& GetSequences() const;

    // Get statistics
    const MinimizerStats& GetStats() const;

    // Get configuration
    const MinimizerConfig& GetConfig() const { return config_; }

    // Serialization
    bool Save(const std::string& path) const;
    static std::unique_ptr<MinimizerIndex> Load(const std::string& path);

    // Check if index is empty
    bool Empty() const;

    // Get total number of minimizers in index
    size_t Size() const;

private:
    MinimizerConfig config_;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Extract minimizers from a sequence (standalone function)
//
// Parameters:
//   sequence: DNA sequence (A/C/G/T/N)
//   config: minimizer extraction parameters
//
// Returns: vector of minimizers sorted by position
std::vector<Minimizer> ExtractMinimizers(
    std::string_view sequence,
    const MinimizerConfig& config = {});

// Extract minimizers into a pre-allocated scratch buffer (zero-allocation hot path)
//
// Parameters:
//   sequence: DNA sequence (A/C/G/T/N)
//   config: minimizer extraction parameters
//   out: scratch buffer to fill with results (cleared, then filled)
//
// Use this version when calling repeatedly to avoid heap allocations.
void ExtractMinimizersInto(
    std::string_view sequence,
    const MinimizerConfig& config,
    core::ScratchBuffer<Minimizer>& out);

// Hash a k-mer using the invertible hash function
uint64_t HashKmer(uint64_t kmer, uint64_t seed);

// Get canonical k-mer (min of forward and reverse complement)
uint64_t CanonicalKmer(uint64_t kmer, uint8_t k);

// Reverse complement a k-mer
uint64_t ReverseComplementKmer(uint64_t kmer, uint8_t k);

}  // namespace llmap::classical
