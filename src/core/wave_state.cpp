// LLmap — WaveState implementation.

#include "core/wave_state.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace llmap {

WaveState::WaveState(std::uint32_t n_reads, WaveLevel initial_level)
    : n_reads_(n_reads),
      read_offsets_(n_reads + 1, 0),
      current_level_(n_reads, initial_level),
      collapsed_(n_reads, false),
      collapsed_bucket_(n_reads, 0) {}

void WaveState::reserve_entries(std::size_t total_entries) {
    bucket_indices_.reserve(total_entries);
    probabilities_.reserve(total_entries);
}

void WaveState::set_read_candidates(std::uint32_t read_idx, std::span<const BucketProb> candidates) {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }

    // CSR requires rebuilding when modifying a single row (expensive but correct).
    // In production, batch-build is preferred; this is for flexibility during dev.

    std::uint32_t old_start = read_offsets_[read_idx];
    std::uint32_t old_end = read_offsets_[read_idx + 1];
    std::uint32_t old_count = old_end - old_start;
    auto new_count = static_cast<std::uint32_t>(candidates.size());

    std::int64_t delta = static_cast<std::int64_t>(new_count) - static_cast<std::int64_t>(old_count);

    if (delta == 0) {
        // Same size: overwrite in place
        for (std::uint32_t i = 0; i < new_count; ++i) {
            bucket_indices_[old_start + i] = candidates[i].bucket_id;
            probabilities_[old_start + i] = candidates[i].probability;
        }
    } else if (delta > 0) {
        // Growing: insert space
        bucket_indices_.insert(bucket_indices_.begin() + old_end,
                               static_cast<std::size_t>(delta), 0);
        probabilities_.insert(probabilities_.begin() + old_end,
                              static_cast<std::size_t>(delta), 0.0f);

        for (std::uint32_t i = 0; i < new_count; ++i) {
            bucket_indices_[old_start + i] = candidates[i].bucket_id;
            probabilities_[old_start + i] = candidates[i].probability;
        }

        for (std::uint32_t r = read_idx + 1; r <= n_reads_; ++r) {
            read_offsets_[r] += static_cast<std::uint32_t>(delta);
        }
    } else {
        // Shrinking: erase space
        auto abs_delta = static_cast<std::size_t>(-delta);
        bucket_indices_.erase(bucket_indices_.begin() + old_start,
                              bucket_indices_.begin() + old_start + abs_delta);
        probabilities_.erase(probabilities_.begin() + old_start,
                             probabilities_.begin() + old_start + abs_delta);

        for (std::uint32_t i = 0; i < new_count; ++i) {
            bucket_indices_[old_start + i] = candidates[i].bucket_id;
            probabilities_[old_start + i] = candidates[i].probability;
        }

        for (std::uint32_t r = read_idx + 1; r <= n_reads_; ++r) {
            read_offsets_[r] -= static_cast<std::uint32_t>(abs_delta);
        }
    }

    device_synced_ = false;
}

std::span<const std::uint32_t> WaveState::bucket_indices_for_read(std::uint32_t read_idx) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    std::uint32_t start = read_offsets_[read_idx];
    std::uint32_t end = read_offsets_[read_idx + 1];
    return {bucket_indices_.data() + start, end - start};
}

std::span<const float> WaveState::probabilities_for_read(std::uint32_t read_idx) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    std::uint32_t start = read_offsets_[read_idx];
    std::uint32_t end = read_offsets_[read_idx + 1];
    return {probabilities_.data() + start, end - start};
}

std::size_t WaveState::find_bucket_in_read(std::uint32_t read_idx, std::uint32_t bucket_id) const {
    std::uint32_t start = read_offsets_[read_idx];
    std::uint32_t end = read_offsets_[read_idx + 1];

    // Binary search since bucket_indices should be sorted within each read
    auto it = std::lower_bound(bucket_indices_.begin() + start,
                               bucket_indices_.begin() + end,
                               bucket_id);
    if (it != bucket_indices_.begin() + end && *it == bucket_id) {
        return static_cast<std::size_t>(it - bucket_indices_.begin());
    }
    return static_cast<std::size_t>(end);
}

float WaveState::get_probability(std::uint32_t read_idx, std::uint32_t bucket_id) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    std::size_t idx = find_bucket_in_read(read_idx, bucket_id);
    std::uint32_t end = read_offsets_[read_idx + 1];
    if (idx < end) {
        return probabilities_[idx];
    }
    return 0.0f;
}

void WaveState::update_probability(std::uint32_t read_idx, std::uint32_t bucket_id, float new_prob) {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    std::size_t idx = find_bucket_in_read(read_idx, bucket_id);
    std::uint32_t end = read_offsets_[read_idx + 1];
    if (idx < end) {
        probabilities_[idx] = new_prob;
        device_synced_ = false;
    }
}

void WaveState::normalize_read(std::uint32_t read_idx) {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }

    std::uint32_t start = read_offsets_[read_idx];
    std::uint32_t end = read_offsets_[read_idx + 1];
    if (start == end) {
        return;
    }

    float sum = 0.0f;
    for (std::uint32_t i = start; i < end; ++i) {
        sum += probabilities_[i];
    }

    if (sum > 0.0f) {
        for (std::uint32_t i = start; i < end; ++i) {
            probabilities_[i] /= sum;
        }
        device_synced_ = false;
    }
}

WaveLevel WaveState::get_level(std::uint32_t read_idx) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    return current_level_[read_idx];
}

void WaveState::set_level(std::uint32_t read_idx, WaveLevel level) {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    current_level_[read_idx] = level;
}

bool WaveState::is_collapsed(std::uint32_t read_idx) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    return collapsed_[read_idx];
}

void WaveState::set_collapsed(std::uint32_t read_idx, bool collapsed) {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    collapsed_[read_idx] = collapsed;
}

std::uint32_t WaveState::collapsed_bucket(std::uint32_t read_idx) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    return collapsed_bucket_[read_idx];
}

bool WaveState::check_collapse(std::uint32_t read_idx, float threshold) const {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }
    if (collapsed_[read_idx]) {
        return true;
    }

    std::uint32_t start = read_offsets_[read_idx];
    std::uint32_t end = read_offsets_[read_idx + 1];
    if (start == end) {
        return false;
    }

    float max_prob = *std::max_element(probabilities_.begin() + start,
                                        probabilities_.begin() + end);
    return max_prob >= threshold;
}

std::uint32_t WaveState::collapse_read(std::uint32_t read_idx) {
    if (read_idx >= n_reads_) {
        throw std::out_of_range("read_idx out of range");
    }

    std::uint32_t start = read_offsets_[read_idx];
    std::uint32_t end = read_offsets_[read_idx + 1];
    if (start == end) {
        throw std::runtime_error("cannot collapse read with no candidates");
    }

    auto max_it = std::max_element(probabilities_.begin() + start,
                                   probabilities_.begin() + end);
    std::size_t max_idx = static_cast<std::size_t>(max_it - probabilities_.begin());

    collapsed_[read_idx] = true;
    collapsed_bucket_[read_idx] = bucket_indices_[max_idx];

    return collapsed_bucket_[read_idx];
}

std::uint32_t WaveState::count_collapsed() const noexcept {
    return static_cast<std::uint32_t>(
        std::count(collapsed_.begin(), collapsed_.end(), true));
}

std::uint32_t WaveState::count_active() const noexcept {
    return n_reads_ - count_collapsed();
}

std::span<const std::uint32_t> WaveState::read_offsets() const noexcept {
    return read_offsets_;
}

std::span<const std::uint32_t> WaveState::bucket_indices() const noexcept {
    return bucket_indices_;
}

std::span<const float> WaveState::probabilities() const noexcept {
    return probabilities_;
}

bool WaveState::validate() const noexcept {
    // Empty state (default constructed) is valid
    if (n_reads_ == 0) {
        return read_offsets_.empty() &&
               bucket_indices_.empty() &&
               probabilities_.empty() &&
               current_level_.empty() &&
               collapsed_.empty() &&
               collapsed_bucket_.empty();
    }

    if (read_offsets_.size() != n_reads_ + 1) {
        return false;
    }
    if (current_level_.size() != n_reads_) {
        return false;
    }
    if (collapsed_.size() != n_reads_) {
        return false;
    }
    if (collapsed_bucket_.size() != n_reads_) {
        return false;
    }
    if (bucket_indices_.size() != probabilities_.size()) {
        return false;
    }

    // Check offsets are monotonically non-decreasing and end at correct value
    if (read_offsets_[0] != 0) {
        return false;
    }
    for (std::uint32_t r = 0; r < n_reads_; ++r) {
        if (read_offsets_[r] > read_offsets_[r + 1]) {
            return false;
        }
    }
    if (read_offsets_[n_reads_] != bucket_indices_.size()) {
        return false;
    }

    // Check bucket indices are sorted within each read (for binary search)
    for (std::uint32_t r = 0; r < n_reads_; ++r) {
        std::uint32_t start = read_offsets_[r];
        std::uint32_t end = read_offsets_[r + 1];
        for (std::uint32_t i = start + 1; i < end; ++i) {
            if (bucket_indices_[i - 1] >= bucket_indices_[i]) {
                return false;
            }
        }
    }

    return true;
}

void WaveState::sync_to_device() {
    // GPU stub: no-op until CUDA integration
    device_synced_ = true;
}

void WaveState::sync_from_device() {
    // GPU stub: no-op until CUDA integration
}

}  // namespace llmap
