// LLmap — Reference Index for Stage 2 WaveCollapse.
//
// The reference index bundles everything needed for reference-based alignment:
//   - BucketPyramid hierarchy (L0→L1→L2)
//   - Pre-computed embeddings for each bucket level
//   - Reference metadata (FASTA checksums, target names, sizes)
//   - Biology prior hints from Claude Session A
//
// Built by `llmap index`, loaded by `llmap align`.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "core/bucket_pyramid.h"

namespace llmap {

// Reference sequence metadata (per target/chromosome)
struct ReferenceTarget {
    std::string name;
    std::uint64_t length{0};
    std::string md5;  // MD5 checksum for BAM @SQ line compatibility
};

// Index build configuration
struct ReferenceIndexConfig {
    std::uint64_t l1_granularity = 5'000'000;    // 5 MB default
    std::uint64_t l2_granularity = 50'000;       // 50 kb default
    std::size_t embedding_dim = 256;              // Matches foundation model
    bool include_embeddings = true;               // False for CPU-only fallback
    std::string reference_version;                // e.g., "GRCh38.p14"
};

// Statistics from index building
struct ReferenceIndexStats {
    std::size_t num_targets = 0;
    std::uint64_t total_length = 0;
    std::size_t l0_buckets = 0;
    std::size_t l1_buckets = 0;
    std::size_t l2_buckets = 0;
    float build_time_seconds = 0.0f;
    float embedding_time_seconds = 0.0f;
    std::uint64_t index_size_bytes = 0;
};

// Forward declaration for PIMPL
class ReferenceIndexImpl;

// Main reference index class
// Contains bucket pyramid + embeddings + metadata
// Thread-safe for read operations after construction
class ReferenceIndex {
public:
    // Builder pattern for index construction
    class Builder {
    public:
        explicit Builder(const ReferenceIndexConfig& config = {});
        ~Builder();

        Builder(const Builder&) = delete;
        Builder& operator=(const Builder&) = delete;
        Builder(Builder&&) noexcept;
        Builder& operator=(Builder&&) noexcept;

        // Add reference targets from FASTA
        Builder& AddTarget(const ReferenceTarget& target);
        Builder& AddTargets(std::span<const ReferenceTarget> targets);

        // Set bucket embeddings (pre-computed externally via BucketEmbedder)
        // Embeddings must be [num_buckets x embedding_dim] in row-major order
        Builder& SetL0Embeddings(std::span<const float> embeddings);
        Builder& SetL1Embeddings(std::span<const float> embeddings);
        Builder& SetL2Embeddings(std::span<const float> embeddings);

        // Set biology prior hints (from Claude Session A JSON)
        Builder& SetBiologyPrior(const std::filesystem::path& json_path);
        Builder& SetBiologyHint(std::uint32_t l2_bucket_id, const BiologyHint& hint);

        // Finalize and build the index
        // Returns nullptr on error (check LastError())
        std::unique_ptr<ReferenceIndex> Build();

        // Build stats (available after Build())
        const ReferenceIndexStats& Stats() const { return stats_; }

        std::string LastError() const { return last_error_; }

    private:
        ReferenceIndexConfig config_;
        std::vector<ReferenceTarget> targets_;
        std::vector<float> l0_embeddings_;
        std::vector<float> l1_embeddings_;
        std::vector<float> l2_embeddings_;
        std::unordered_map<std::uint32_t, BiologyHint> biology_hints_;
        ReferenceIndexStats stats_;
        std::string last_error_;
    };

    ~ReferenceIndex();

    ReferenceIndex(const ReferenceIndex&) = delete;
    ReferenceIndex& operator=(const ReferenceIndex&) = delete;
    ReferenceIndex(ReferenceIndex&&) noexcept;
    ReferenceIndex& operator=(ReferenceIndex&&) noexcept;

    // Load from serialized index file
    static std::unique_ptr<ReferenceIndex> Load(const std::filesystem::path& path);

    // Save to file
    bool Save(const std::filesystem::path& path) const;

    // Access bucket pyramid
    const BucketPyramid& Pyramid() const;

    // Access reference targets
    std::span<const ReferenceTarget> Targets() const;
    std::optional<ReferenceTarget> FindTarget(std::string_view name) const;
    std::size_t NumTargets() const;

    // Access embeddings (row-major: [num_buckets x embedding_dim])
    std::span<const float> L0Embeddings() const;
    std::span<const float> L1Embeddings() const;
    std::span<const float> L2Embeddings() const;
    std::size_t EmbeddingDim() const;
    bool HasEmbeddings() const;

    // Get embedding for specific bucket (returns slice of embedding_dim floats)
    std::span<const float> GetL0Embedding(std::uint32_t bucket_id) const;
    std::span<const float> GetL1Embedding(std::uint32_t bucket_id) const;
    std::span<const float> GetL2Embedding(std::uint32_t bucket_id) const;

    // Lookup L2 bucket by genomic position
    std::optional<std::uint32_t> FindL2Bucket(
        std::string_view target_name,
        std::uint64_t position) const;

    // Get L2 buckets overlapping a range
    std::vector<std::uint32_t> FindL2BucketsInRange(
        std::string_view target_name,
        std::uint64_t start,
        std::uint64_t end) const;

    // Index metadata
    const ReferenceIndexConfig& Config() const;
    const ReferenceIndexStats& Stats() const;
    std::string ReferenceVersion() const;

    // File format version for compatibility checking
    static constexpr std::uint32_t FormatVersion() { return 1; }

private:
    friend class Builder;
    explicit ReferenceIndex(std::unique_ptr<ReferenceIndexImpl> impl);

    std::unique_ptr<ReferenceIndexImpl> impl_;
};

// Utility: compute total reference length from targets
std::uint64_t ComputeTotalLength(std::span<const ReferenceTarget> targets);

// Utility: estimate number of L1/L2 buckets for a given reference
std::pair<std::size_t, std::size_t> EstimateBucketCounts(
    std::uint64_t total_length,
    const ReferenceIndexConfig& config);

}  // namespace llmap
