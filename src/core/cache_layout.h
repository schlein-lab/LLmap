// LLmap — Cache-friendly data layouts for hot paths.
//
// Provides:
//   - SoA (Structure of Arrays) containers for cache-efficient access
//   - Prefetch utilities for sequential and strided access patterns
//   - Cache line alignment helpers
//
// Background: AoS (Array of Structures) like std::vector<Anchor> interleaves
// fields. When processing only one field at a time (e.g., sorting by ref_pos),
// this wastes cache bandwidth loading unused fields. SoA stores each field
// contiguously, improving cache utilization for field-specific operations.

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace llmap::core {

// Common cache line size (64 bytes on x86/ARM)
inline constexpr size_t kCacheLineSize = 64;

// Prefetch distance tuning constants
// These are empirically tuned for typical sequence processing patterns
inline constexpr size_t kPrefetchDistance = 8;       // Elements ahead to prefetch
inline constexpr size_t kPrefetchStride = 512;       // Bytes for strided prefetch

// Prefetch hints for sequential access patterns
//
// Call before entering a loop that will access memory sequentially.
// The CPU will fetch cache lines ahead of actual access.
template <typename T>
inline void PrefetchForRead(const T* ptr) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(ptr, 0, 3);  // Read, high locality
#elif defined(_MSC_VER) && defined(_M_X64)
    _mm_prefetch(reinterpret_cast<const char*>(ptr), _MM_HINT_T0);
#endif
}

template <typename T>
inline void PrefetchForWrite(T* ptr) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(ptr, 1, 3);  // Write, high locality
#elif defined(_MSC_VER) && defined(_M_X64)
    _mm_prefetch(reinterpret_cast<const char*>(ptr), _MM_HINT_T0);
#endif
}

// Prefetch for temporal access (data will be reused soon)
template <typename T>
inline void PrefetchTemporal(const T* ptr) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(ptr, 0, 3);
#endif
}

// Prefetch for non-temporal access (data will not be reused)
template <typename T>
inline void PrefetchNonTemporal(const T* ptr) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(ptr, 0, 0);
#endif
}

// Prefetch a range of memory for sequential read access
inline void PrefetchRange(const void* ptr, size_t bytes) {
    const char* p = static_cast<const char*>(ptr);
    const char* end = p + bytes;
    while (p < end) {
        PrefetchForRead(p);
        p += kCacheLineSize;
    }
}

// Allocate cache-line aligned memory
inline void* AllocateAligned(size_t size, size_t alignment = kCacheLineSize) {
    void* ptr = std::aligned_alloc(alignment, (size + alignment - 1) & ~(alignment - 1));
    if (!ptr) throw std::bad_alloc();
    return ptr;
}

// Free aligned memory
inline void FreeAligned(void* ptr) {
    std::free(ptr);
}

// Check if pointer is cache-line aligned
template <typename T>
inline bool IsCacheAligned(const T* ptr) {
    return (reinterpret_cast<uintptr_t>(ptr) & (kCacheLineSize - 1)) == 0;
}

// SoA container for minimizer data
//
// Instead of vector<Minimizer> where each element has {hash, pos, is_reverse},
// this stores three contiguous arrays: hashes[], positions[], strands[].
//
// Benefits:
// - Sorting by hash only touches the hash array + index array
// - Sequential hash comparisons have better cache locality
// - Position lookups after sorting still access contiguous memory
class MinimizerSoA {
public:
    MinimizerSoA() = default;

    explicit MinimizerSoA(size_t capacity) { Reserve(capacity); }

    ~MinimizerSoA() { Clear(); }

    MinimizerSoA(MinimizerSoA&& other) noexcept
        : hashes_(other.hashes_),
          positions_(other.positions_),
          strands_(other.strands_),
          size_(other.size_),
          capacity_(other.capacity_) {
        other.hashes_ = nullptr;
        other.positions_ = nullptr;
        other.strands_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    MinimizerSoA& operator=(MinimizerSoA&& other) noexcept {
        if (this != &other) {
            Clear();
            hashes_ = other.hashes_;
            positions_ = other.positions_;
            strands_ = other.strands_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.hashes_ = nullptr;
            other.positions_ = nullptr;
            other.strands_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    MinimizerSoA(const MinimizerSoA&) = delete;
    MinimizerSoA& operator=(const MinimizerSoA&) = delete;

    void Reserve(size_t new_capacity) {
        if (new_capacity <= capacity_) return;

        auto* new_hashes = static_cast<uint64_t*>(
            AllocateAligned(new_capacity * sizeof(uint64_t)));
        auto* new_positions = static_cast<uint32_t*>(
            AllocateAligned(new_capacity * sizeof(uint32_t)));
        auto* new_strands = static_cast<uint8_t*>(
            AllocateAligned(new_capacity * sizeof(uint8_t)));

        if (size_ > 0) {
            std::memcpy(new_hashes, hashes_, size_ * sizeof(uint64_t));
            std::memcpy(new_positions, positions_, size_ * sizeof(uint32_t));
            std::memcpy(new_strands, strands_, size_ * sizeof(uint8_t));
        }

        FreeAligned(hashes_);
        FreeAligned(positions_);
        FreeAligned(strands_);

        hashes_ = new_hashes;
        positions_ = new_positions;
        strands_ = new_strands;
        capacity_ = new_capacity;
    }

    void Push(uint64_t hash, uint32_t pos, bool is_reverse) {
        if (size_ == capacity_) {
            Reserve(capacity_ == 0 ? 64 : capacity_ * 2);
        }
        hashes_[size_] = hash;
        positions_[size_] = pos;
        strands_[size_] = is_reverse ? 1 : 0;
        ++size_;
    }

    void Clear() {
        FreeAligned(hashes_);
        FreeAligned(positions_);
        FreeAligned(strands_);
        hashes_ = nullptr;
        positions_ = nullptr;
        strands_ = nullptr;
        size_ = 0;
        capacity_ = 0;
    }

    void Reset() { size_ = 0; }

    size_t Size() const { return size_; }
    size_t Capacity() const { return capacity_; }
    bool Empty() const { return size_ == 0; }

    uint64_t Hash(size_t i) const { return hashes_[i]; }
    uint32_t Position(size_t i) const { return positions_[i]; }
    bool IsReverse(size_t i) const { return strands_[i] != 0; }

    std::span<const uint64_t> Hashes() const { return {hashes_, size_}; }
    std::span<const uint32_t> Positions() const { return {positions_, size_}; }
    std::span<const uint8_t> Strands() const { return {strands_, size_}; }

    uint64_t* HashData() { return hashes_; }
    uint32_t* PositionData() { return positions_; }
    uint8_t* StrandData() { return strands_; }

    void PrefetchNext(size_t current_idx) const {
        if (current_idx + kPrefetchDistance < size_) {
            PrefetchForRead(&hashes_[current_idx + kPrefetchDistance]);
        }
    }

private:
    uint64_t* hashes_ = nullptr;
    uint32_t* positions_ = nullptr;
    uint8_t* strands_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
};

// SoA container for anchor data (used in chaining)
//
// Stores ref_id, ref_pos, query_pos, strand in separate arrays.
// Sorting by (ref_id, ref_pos) for chaining only touches those arrays.
class AnchorSoA {
public:
    AnchorSoA() = default;

    explicit AnchorSoA(size_t capacity) { Reserve(capacity); }

    ~AnchorSoA() { Clear(); }

    AnchorSoA(AnchorSoA&& other) noexcept
        : ref_ids_(other.ref_ids_),
          ref_positions_(other.ref_positions_),
          query_positions_(other.query_positions_),
          strands_(other.strands_),
          size_(other.size_),
          capacity_(other.capacity_) {
        other.ref_ids_ = nullptr;
        other.ref_positions_ = nullptr;
        other.query_positions_ = nullptr;
        other.strands_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    AnchorSoA& operator=(AnchorSoA&& other) noexcept {
        if (this != &other) {
            Clear();
            ref_ids_ = other.ref_ids_;
            ref_positions_ = other.ref_positions_;
            query_positions_ = other.query_positions_;
            strands_ = other.strands_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.ref_ids_ = nullptr;
            other.ref_positions_ = nullptr;
            other.query_positions_ = nullptr;
            other.strands_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    AnchorSoA(const AnchorSoA&) = delete;
    AnchorSoA& operator=(const AnchorSoA&) = delete;

    void Reserve(size_t new_capacity) {
        if (new_capacity <= capacity_) return;

        auto* new_ref_ids = static_cast<uint32_t*>(
            AllocateAligned(new_capacity * sizeof(uint32_t)));
        auto* new_ref_positions = static_cast<uint32_t*>(
            AllocateAligned(new_capacity * sizeof(uint32_t)));
        auto* new_query_positions = static_cast<uint32_t*>(
            AllocateAligned(new_capacity * sizeof(uint32_t)));
        auto* new_strands = static_cast<uint8_t*>(
            AllocateAligned(new_capacity * sizeof(uint8_t)));

        if (size_ > 0) {
            std::memcpy(new_ref_ids, ref_ids_, size_ * sizeof(uint32_t));
            std::memcpy(new_ref_positions, ref_positions_, size_ * sizeof(uint32_t));
            std::memcpy(new_query_positions, query_positions_, size_ * sizeof(uint32_t));
            std::memcpy(new_strands, strands_, size_ * sizeof(uint8_t));
        }

        FreeAligned(ref_ids_);
        FreeAligned(ref_positions_);
        FreeAligned(query_positions_);
        FreeAligned(strands_);

        ref_ids_ = new_ref_ids;
        ref_positions_ = new_ref_positions;
        query_positions_ = new_query_positions;
        strands_ = new_strands;
        capacity_ = new_capacity;
    }

    void Push(uint32_t ref_id, uint32_t ref_pos, uint32_t query_pos, bool same_strand) {
        if (size_ == capacity_) {
            Reserve(capacity_ == 0 ? 64 : capacity_ * 2);
        }
        ref_ids_[size_] = ref_id;
        ref_positions_[size_] = ref_pos;
        query_positions_[size_] = query_pos;
        strands_[size_] = same_strand ? 1 : 0;
        ++size_;
    }

    void Clear() {
        FreeAligned(ref_ids_);
        FreeAligned(ref_positions_);
        FreeAligned(query_positions_);
        FreeAligned(strands_);
        ref_ids_ = nullptr;
        ref_positions_ = nullptr;
        query_positions_ = nullptr;
        strands_ = nullptr;
        size_ = 0;
        capacity_ = 0;
    }

    void Reset() { size_ = 0; }

    size_t Size() const { return size_; }
    size_t Capacity() const { return capacity_; }
    bool Empty() const { return size_ == 0; }

    uint32_t RefId(size_t i) const { return ref_ids_[i]; }
    uint32_t RefPosition(size_t i) const { return ref_positions_[i]; }
    uint32_t QueryPosition(size_t i) const { return query_positions_[i]; }
    bool IsSameStrand(size_t i) const { return strands_[i] != 0; }

    std::span<const uint32_t> RefIds() const { return {ref_ids_, size_}; }
    std::span<const uint32_t> RefPositions() const { return {ref_positions_, size_}; }
    std::span<const uint32_t> QueryPositions() const { return {query_positions_, size_}; }
    std::span<const uint8_t> Strands() const { return {strands_, size_}; }

    uint32_t* RefIdData() { return ref_ids_; }
    uint32_t* RefPositionData() { return ref_positions_; }
    uint32_t* QueryPositionData() { return query_positions_; }
    uint8_t* StrandData() { return strands_; }

    // Sort by (ref_id, ref_pos) using index array for cache efficiency
    // Returns permutation array that can be applied to reorder
    std::vector<size_t> SortPermutation() const {
        std::vector<size_t> perm(size_);
        for (size_t i = 0; i < size_; ++i) perm[i] = i;

        std::sort(perm.begin(), perm.end(), [this](size_t a, size_t b) {
            if (ref_ids_[a] != ref_ids_[b]) return ref_ids_[a] < ref_ids_[b];
            return ref_positions_[a] < ref_positions_[b];
        });

        return perm;
    }

    void PrefetchNext(size_t current_idx) const {
        if (current_idx + kPrefetchDistance < size_) {
            size_t next = current_idx + kPrefetchDistance;
            PrefetchForRead(&ref_ids_[next]);
            PrefetchForRead(&ref_positions_[next]);
        }
    }

private:
    uint32_t* ref_ids_ = nullptr;
    uint32_t* ref_positions_ = nullptr;
    uint32_t* query_positions_ = nullptr;
    uint8_t* strands_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
};

// SoA container for DP state (used in chaining)
//
// Stores score and predecessor arrays separately for cache efficiency
// during the DP computation where scores are accessed linearly.
class DPStateSoA {
public:
    DPStateSoA() = default;

    explicit DPStateSoA(size_t capacity) { Resize(capacity); }

    ~DPStateSoA() {
        FreeAligned(scores_);
        FreeAligned(predecessors_);
    }

    DPStateSoA(DPStateSoA&& other) noexcept
        : scores_(other.scores_),
          predecessors_(other.predecessors_),
          size_(other.size_),
          capacity_(other.capacity_) {
        other.scores_ = nullptr;
        other.predecessors_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    DPStateSoA& operator=(DPStateSoA&& other) noexcept {
        if (this != &other) {
            FreeAligned(scores_);
            FreeAligned(predecessors_);
            scores_ = other.scores_;
            predecessors_ = other.predecessors_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.scores_ = nullptr;
            other.predecessors_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    DPStateSoA(const DPStateSoA&) = delete;
    DPStateSoA& operator=(const DPStateSoA&) = delete;

    void Resize(size_t new_size) {
        if (new_size > capacity_) {
            FreeAligned(scores_);
            FreeAligned(predecessors_);

            capacity_ = new_size;
            scores_ = static_cast<int32_t*>(
                AllocateAligned(capacity_ * sizeof(int32_t)));
            predecessors_ = static_cast<int32_t*>(
                AllocateAligned(capacity_ * sizeof(int32_t)));
        }
        size_ = new_size;
    }

    void Initialize() {
        std::memset(scores_, 0, size_ * sizeof(int32_t));
        std::fill_n(predecessors_, size_, -1);
    }

    size_t Size() const { return size_; }

    int32_t Score(size_t i) const { return scores_[i]; }
    int32_t Predecessor(size_t i) const { return predecessors_[i]; }

    void SetScore(size_t i, int32_t score) { scores_[i] = score; }
    void SetPredecessor(size_t i, int32_t pred) { predecessors_[i] = pred; }

    int32_t* ScoreData() { return scores_; }
    int32_t* PredecessorData() { return predecessors_; }

    void PrefetchScores(size_t current_idx) const {
        if (current_idx + kPrefetchDistance < size_) {
            PrefetchForRead(&scores_[current_idx + kPrefetchDistance]);
        }
    }

private:
    int32_t* scores_ = nullptr;
    int32_t* predecessors_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
};

// Cache-oblivious tile processor for 2D DP
//
// Processes a matrix in cache-friendly tiles rather than row-by-row.
// Useful for WFA2 and other gap-affine DP algorithms.
template <size_t TileSize = 64>
class TiledProcessor {
public:
    // Process tiles of a matrix (rows x cols) calling func(row, col) for each cell
    template <typename F>
    static void Process(size_t rows, size_t cols, F&& func) {
        for (size_t tile_row = 0; tile_row < rows; tile_row += TileSize) {
            for (size_t tile_col = 0; tile_col < cols; tile_col += TileSize) {
                size_t row_end = std::min(tile_row + TileSize, rows);
                size_t col_end = std::min(tile_col + TileSize, cols);

                for (size_t r = tile_row; r < row_end; ++r) {
                    for (size_t c = tile_col; c < col_end; ++c) {
                        func(r, c);
                    }
                }
            }
        }
    }

    // Process with diagonal traversal (for DP with diagonal dependencies)
    template <typename F>
    static void ProcessDiagonal(size_t rows, size_t cols, F&& func) {
        size_t num_diags = rows + cols - 1;
        for (size_t d = 0; d < num_diags; ++d) {
            size_t r_start = (d < cols) ? 0 : d - cols + 1;
            size_t r_end = std::min(d + 1, rows);

            for (size_t r = r_start; r < r_end; ++r) {
                size_t c = d - r;
                func(r, c);
            }
        }
    }
};

// Statistics for cache layout operations
struct CacheLayoutStats {
    size_t aos_to_soa_conversions = 0;
    size_t prefetch_hints = 0;
    size_t aligned_allocations = 0;
    size_t total_bytes_allocated = 0;
};

}  // namespace llmap::core
