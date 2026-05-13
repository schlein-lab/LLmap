// LLmap — hierarchical bucket pyramid (RG-flow).
//
// Buckets form a 3-level hierarchy:
//   L0: ~1k coarse buckets (chromosomes, repeat families, pseudogene clades)
//   L1: ~600 medium buckets (5 MB windows)
//   L2: ~60k fine buckets (50 kb windows)
//   L3: continuous position (not stored, computed at WFA2 extension time)
//
// The pyramid stores parent indices linking child to parent at each level,
// enabling coarse→fine refinement during WaveCollapse iterations.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#include <span>
#include <filesystem>

namespace llmap {

// --- bucket structs per level ---------------------------------------------

struct L0Bucket {
    std::uint32_t id{0};
    std::string name;  // e.g., "chr14", "repeat:ALU", "clade:IGHG"
    std::uint64_t total_span{0};  // sum of child spans
};

struct L1Bucket {
    std::uint32_t id{0};
    std::string target_id;  // e.g., "chr14"
    std::uint64_t start{0};
    std::uint64_t end{0};  // [start, end), 5 MB nominal
};

struct L2Bucket {
    std::uint32_t id{0};
    std::string target_id;
    std::uint64_t start{0};
    std::uint64_t end{0};  // [start, end), 50 kb nominal
};

// --- biology prior hint from Claude Session A -----------------------------

struct BiologyHint {
    float prior_weight{1.0f};
    std::optional<std::string> annotation;
    std::optional<std::uint32_t> paralog_partner_bucket;
    float expected_coverage_multiplier{1.0f};
};

// --- the pyramid ----------------------------------------------------------

class BucketPyramid {
public:
    BucketPyramid() = default;

    // Build pyramid from reference metadata (placeholder for now)
    void add_l0_bucket(L0Bucket bucket);
    void add_l1_bucket(L1Bucket bucket, std::uint32_t parent_l0_id);
    void add_l2_bucket(L2Bucket bucket, std::uint32_t parent_l1_id);

    // Set biology hint for a bucket (keyed by L2 bucket id)
    void set_biology_hint(std::uint32_t l2_bucket_id, BiologyHint hint);

    // Accessors
    [[nodiscard]] std::span<const L0Bucket> l0_buckets() const noexcept;
    [[nodiscard]] std::span<const L1Bucket> l1_buckets() const noexcept;
    [[nodiscard]] std::span<const L2Bucket> l2_buckets() const noexcept;

    [[nodiscard]] std::uint32_t l1_parent(std::uint32_t l1_id) const;
    [[nodiscard]] std::uint32_t l2_parent(std::uint32_t l2_id) const;

    [[nodiscard]] std::optional<BiologyHint> get_biology_hint(std::uint32_t l2_bucket_id) const;

    [[nodiscard]] std::size_t l0_count() const noexcept { return l0_.size(); }
    [[nodiscard]] std::size_t l1_count() const noexcept { return l1_.size(); }
    [[nodiscard]] std::size_t l2_count() const noexcept { return l2_.size(); }

    // Validation: check parent indices are consistent
    [[nodiscard]] bool validate() const noexcept;

    // Serialization (binary format for speed)
    void serialize(const std::filesystem::path& path) const;
    static BucketPyramid deserialize(const std::filesystem::path& path);

    // Round-trip equality check (for testing)
    [[nodiscard]] bool operator==(const BucketPyramid& other) const noexcept;

private:
    std::vector<L0Bucket> l0_;
    std::vector<L1Bucket> l1_;
    std::vector<L2Bucket> l2_;

    std::vector<std::uint32_t> l1_to_l0_;  // l1_to_l0_[l1_id] = parent l0_id
    std::vector<std::uint32_t> l2_to_l1_;  // l2_to_l1_[l2_id] = parent l1_id

    std::unordered_map<std::uint32_t, BiologyHint> biology_hints_;

    // Future: embeddings (stubbed out until ONNX integration)
    // std::vector<float> l0_embeddings_;  // flattened (l0_count * embed_dim)
    // std::vector<float> l1_embeddings_;
    // std::vector<float> l2_embeddings_;

    static constexpr std::uint32_t MAGIC = 0x4C4C4D50;  // "LLMP"
    static constexpr std::uint32_t VERSION = 1;
};

}  // namespace llmap
