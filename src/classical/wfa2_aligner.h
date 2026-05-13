#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::classical {

// Gap-affine alignment scoring parameters
struct WFA2Config {
    // Match/mismatch
    int32_t match_score = 0;       // WFA2 uses penalties, not scores (0 for match)
    int32_t mismatch_penalty = 4;  // Penalty for mismatch

    // Affine gap penalties: gap_cost = gap_open + (gap_length - 1) * gap_extend
    int32_t gap_open = 6;          // Penalty for opening a gap
    int32_t gap_extend = 2;        // Penalty per extended base

    // Alignment scope
    bool end_to_end = false;       // true = global, false = extension (semi-global)

    // Memory/computation limits
    int32_t max_alignment_score = 10000;  // Abort if score exceeds this
    size_t max_pattern_length = 100000;   // Max query length
    size_t max_text_length = 100000;      // Max reference length
};

// CIGAR operation types (BAM-style)
enum class CigarOp : uint8_t {
    Match = 0,       // M: alignment match or mismatch
    Insertion = 1,   // I: insertion to reference
    Deletion = 2,    // D: deletion from reference
    SoftClip = 4,    // S: soft clipping (bases present but not aligned)
    Equal = 7,       // =: sequence match
    Diff = 8,        // X: sequence mismatch
};

// A single CIGAR operation
struct CigarElement {
    CigarOp op;
    uint32_t length;

    char OpChar() const;
    std::string ToString() const;
};

// Result of a WFA2 alignment
struct WFA2Result {
    std::vector<CigarElement> cigar;   // CIGAR operations
    int32_t score = 0;                  // Alignment score (lower is better in WFA2)
    int32_t query_start = 0;            // Start position in query
    int32_t query_end = 0;              // End position in query
    int32_t ref_start = 0;              // Start position in reference
    int32_t ref_end = 0;                // End position in reference
    size_t num_matches = 0;             // Number of matching bases
    size_t num_mismatches = 0;          // Number of mismatching bases
    size_t num_insertions = 0;          // Total insertion bases
    size_t num_deletions = 0;           // Total deletion bases
    float identity = 0.0f;              // Fraction of matches
    float alignment_time_us = 0.0f;     // Time for this alignment

    std::string CigarString() const;
    size_t AlignedLength() const;
};

// Forward declaration for PIMPL
class WFA2AlignerImpl;

// WFA2 aligner wrapper
// Wraps WFA2-lib for gap-affine sequence alignment
// Falls back to a simple DP implementation when WFA2-lib is not available
class WFA2Aligner {
public:
    explicit WFA2Aligner(const WFA2Config& config = {});
    ~WFA2Aligner();

    // Non-copyable, movable
    WFA2Aligner(const WFA2Aligner&) = delete;
    WFA2Aligner& operator=(const WFA2Aligner&) = delete;
    WFA2Aligner(WFA2Aligner&&) noexcept;
    WFA2Aligner& operator=(WFA2Aligner&&) noexcept;

    // Align query (pattern) against reference (text)
    // Returns nullopt if alignment fails or exceeds limits
    std::optional<WFA2Result> Align(
        std::string_view query,
        std::string_view reference) const;

    // Extension alignment: extend from anchor positions
    // Aligns query[query_start:] against reference[ref_start:] going forward,
    // or query[:query_end] against reference[:ref_end] going backward
    std::optional<WFA2Result> ExtendRight(
        std::string_view query,
        std::string_view reference,
        int32_t query_start,
        int32_t ref_start) const;

    std::optional<WFA2Result> ExtendLeft(
        std::string_view query,
        std::string_view reference,
        int32_t query_end,
        int32_t ref_end) const;

    // Batch alignment (for efficiency)
    std::vector<std::optional<WFA2Result>> AlignBatch(
        std::span<const std::string_view> queries,
        std::span<const std::string_view> references) const;

    // Configuration
    const WFA2Config& Config() const { return config_; }
    bool IsNativeWFA2() const;  // True if using WFA2-lib, false if fallback

private:
    WFA2Config config_;
    std::unique_ptr<WFA2AlignerImpl> impl_;
};

// Utility: check if WFA2-lib is available at compile time
bool IsWFA2Available();

// Utility: compute alignment identity from CIGAR
float ComputeIdentity(const std::vector<CigarElement>& cigar);

// Utility: parse CIGAR string to elements
std::optional<std::vector<CigarElement>> ParseCigar(std::string_view cigar_str);

}  // namespace llmap::classical
