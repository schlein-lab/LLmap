// LLmap — sparse wave-state for Read x Bucket probability matrix.
//
// Compressed Sparse Row (CSR) format stores per-read bucket probabilities
// efficiently: each read retains only top-K candidates (K=10 default).
// GPU mirrors are stubbed until CUDA integration in Phase 1.

#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace llmap {

// Hierarchical level where a read currently operates
enum class WaveLevel : std::uint8_t {
    L0 = 0,
    L1 = 1,
    L2 = 2,
    L3 = 3,
};

// Per-read probability entry (bucket_id, probability) for iteration
struct BucketProb {
    std::uint32_t bucket_id{0};
    float probability{0.0f};
};

// Sparse CSR representation of Read x Bucket probability matrix.
//
// Layout:
//   read_offsets_[r] = start index in bucket_indices_ / probabilities_ for read r
//   read_offsets_[n_reads] = total entries (sentinel)
//   bucket_indices_[i] = bucket id for entry i
//   probabilities_[i] = probability for entry i
//
// For read r, its candidates are:
//   bucket_indices_[read_offsets_[r] .. read_offsets_[r+1])
//   probabilities_[read_offsets_[r] .. read_offsets_[r+1])
class WaveState {
public:
    static constexpr std::uint32_t DEFAULT_TOP_K = 10;

    WaveState() = default;

    // Initialize for n reads, all starting at given level with empty candidates
    explicit WaveState(std::uint32_t n_reads, WaveLevel initial_level = WaveLevel::L0);

    // Reserve capacity hint for total bucket entries
    void reserve_entries(std::size_t total_entries);

    // Set candidates for a read (replaces existing, must be sorted by bucket_id)
    void set_read_candidates(std::uint32_t read_idx, std::span<const BucketProb> candidates);

    // Get candidates for a read
    [[nodiscard]] std::span<const std::uint32_t> bucket_indices_for_read(std::uint32_t read_idx) const;
    [[nodiscard]] std::span<const float> probabilities_for_read(std::uint32_t read_idx) const;

    // Get probability of specific bucket for a read (0.0 if not in candidates)
    [[nodiscard]] float get_probability(std::uint32_t read_idx, std::uint32_t bucket_id) const;

    // Update probability for a specific bucket (no-op if bucket not in candidates)
    void update_probability(std::uint32_t read_idx, std::uint32_t bucket_id, float new_prob);

    // Normalize probabilities for a read so they sum to 1.0
    void normalize_read(std::uint32_t read_idx);

    // Per-read state accessors
    [[nodiscard]] WaveLevel get_level(std::uint32_t read_idx) const;
    void set_level(std::uint32_t read_idx, WaveLevel level);

    [[nodiscard]] bool is_collapsed(std::uint32_t read_idx) const;
    void set_collapsed(std::uint32_t read_idx, bool collapsed);

    [[nodiscard]] std::uint32_t collapsed_bucket(std::uint32_t read_idx) const;

    // Check if read has collapsed (max probability >= threshold)
    [[nodiscard]] bool check_collapse(std::uint32_t read_idx, float threshold = 0.99f) const;

    // Collapse read to its max-probability bucket (marks as collapsed, returns bucket_id)
    std::uint32_t collapse_read(std::uint32_t read_idx);

    // Global accessors
    [[nodiscard]] std::uint32_t n_reads() const noexcept { return n_reads_; }
    [[nodiscard]] std::size_t total_entries() const noexcept { return bucket_indices_.size(); }
    [[nodiscard]] std::uint32_t count_collapsed() const noexcept;
    [[nodiscard]] std::uint32_t count_active() const noexcept;

    // Raw CSR access (for GPU copy / debugging)
    [[nodiscard]] std::span<const std::uint32_t> read_offsets() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> bucket_indices() const noexcept;
    [[nodiscard]] std::span<const float> probabilities() const noexcept;

    // Validation: check CSR invariants
    [[nodiscard]] bool validate() const noexcept;

    // --- GPU stub (Phase 1 will add actual CUDA implementation) ---
    //
    // These placeholders exist so the interface is stable; they are no-ops until
    // CUDA is enabled. Device vectors will be thrust::device_vector<T> when real.

    // Sync host → device (no-op until CUDA)
    void sync_to_device();

    // Sync device → host (no-op until CUDA)
    void sync_from_device();

    // Returns true if device memory is allocated and synced
    [[nodiscard]] bool has_device_copy() const noexcept { return device_synced_; }

private:
    std::uint32_t n_reads_{0};

    // CSR arrays
    std::vector<std::uint32_t> read_offsets_;     // size = n_reads_ + 1
    std::vector<std::uint32_t> bucket_indices_;   // size = sum of candidates across reads
    std::vector<float> probabilities_;            // matched to bucket_indices_

    // Per-read metadata
    std::vector<WaveLevel> current_level_;        // size = n_reads_
    std::vector<bool> collapsed_;                 // size = n_reads_
    std::vector<std::uint32_t> collapsed_bucket_; // bucket id when collapsed, else 0

    // GPU state placeholder
    bool device_synced_{false};

    // Helper: find bucket in read's candidates (returns index in bucket_indices_, or end)
    [[nodiscard]] std::size_t find_bucket_in_read(std::uint32_t read_idx, std::uint32_t bucket_id) const;
};

}  // namespace llmap
