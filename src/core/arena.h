// LLmap — Arena allocator: fast bump allocation for hot paths.
//
// Two main components:
//   1. Arena — bump allocator for short-lived allocations within a single
//      operation (e.g., processing one read). Reset() clears all allocations.
//   2. ScratchBuffer<T> — reusable vector-like container that avoids repeated
//      allocations by maintaining capacity across resets.
//
// Thread safety: Arena is NOT thread-safe. Each thread should have its own.
// Use thread_local instances for per-thread scratch space.

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <span>
#include <type_traits>
#include <vector>

namespace llmap::core {

// Simple bump allocator for short-lived, per-operation allocations.
// All allocations are invalidated on Reset(). Does not call destructors.
class Arena {
public:
    static constexpr size_t kDefaultBlockSize = 64 * 1024;  // 64 KB

    explicit Arena(size_t block_size = kDefaultBlockSize)
        : block_size_(block_size),
          current_block_(nullptr),
          current_pos_(0),
          current_end_(0) {}

    ~Arena() { Clear(); }

    // Non-copyable, non-movable
    Arena(const Arena&) = delete;
    Arena& operator=(const Arena&) = delete;
    Arena(Arena&&) = delete;
    Arena& operator=(Arena&&) = delete;

    // Allocate n bytes with given alignment (default 8)
    void* Allocate(size_t n, size_t align = 8) {
        assert(align > 0 && (align & (align - 1)) == 0);  // power of 2

        size_t aligned_pos = (current_pos_ + align - 1) & ~(align - 1);

        if (aligned_pos + n <= current_end_) {
            void* result = current_block_ + aligned_pos;
            current_pos_ = aligned_pos + n;
            total_allocated_ += n;
            return result;
        }

        return AllocateSlow(n, align);
    }

    // Typed allocation
    template <typename T>
    T* Allocate(size_t count = 1) {
        return static_cast<T*>(Allocate(sizeof(T) * count, alignof(T)));
    }

    // Construct in arena
    template <typename T, typename... Args>
    T* Construct(Args&&... args) {
        T* ptr = Allocate<T>();
        new (ptr) T(std::forward<Args>(args)...);
        return ptr;
    }

    // Allocate an array and return as span
    template <typename T>
    std::span<T> AllocateSpan(size_t count) {
        T* ptr = Allocate<T>(count);
        return std::span<T>(ptr, count);
    }

    // Reset arena - all allocations become invalid
    // Keeps the first block for reuse, releases extras
    void Reset() {
        if (blocks_.empty()) {
            current_pos_ = 0;
            total_allocated_ = 0;
            return;
        }

        // Keep first block, free the rest
        for (size_t i = 1; i < blocks_.size(); ++i) {
            std::free(blocks_[i]);
        }
        blocks_.resize(1);
        current_block_ = blocks_[0];
        current_pos_ = 0;
        current_end_ = block_size_;
        total_allocated_ = 0;
    }

    // Clear arena - release all memory
    void Clear() {
        for (void* block : blocks_) {
            std::free(block);
        }
        blocks_.clear();
        current_block_ = nullptr;
        current_pos_ = 0;
        current_end_ = 0;
        total_allocated_ = 0;
    }

    // Statistics
    size_t BytesAllocated() const { return total_allocated_; }
    size_t BlocksUsed() const { return blocks_.size(); }
    size_t BlockSize() const { return block_size_; }
    size_t CapacityBytes() const { return blocks_.size() * block_size_; }

private:
    void* AllocateSlow(size_t n, size_t align) {
        // Need a new block
        size_t alloc_size = std::max(block_size_, n + align);

        void* new_block = std::aligned_alloc(align, alloc_size);
        if (!new_block) {
            throw std::bad_alloc();
        }
        blocks_.push_back(static_cast<char*>(new_block));

        current_block_ = static_cast<char*>(new_block);
        current_pos_ = n;
        current_end_ = alloc_size;
        total_allocated_ += n;

        return new_block;
    }

    size_t block_size_;
    std::vector<char*> blocks_;
    char* current_block_;
    size_t current_pos_;
    size_t current_end_;
    size_t total_allocated_ = 0;
};

// Reusable scratch buffer - maintains capacity across resets.
// Use for hot paths where the same vector is repeatedly filled and cleared.
template <typename T>
class ScratchBuffer {
public:
    ScratchBuffer() = default;

    explicit ScratchBuffer(size_t initial_capacity) {
        reserve(initial_capacity);
    }

    // Clear contents but keep capacity
    void clear() { size_ = 0; }

    // Resize to new size (does not initialize new elements)
    void resize(size_t new_size) {
        if (new_size > capacity_) {
            grow(new_size);
        }
        size_ = new_size;
    }

    // Resize and zero-initialize new elements
    void resize_zero(size_t new_size) {
        if (new_size > capacity_) {
            grow(new_size);
        }
        if (new_size > size_) {
            std::memset(data_ + size_, 0, (new_size - size_) * sizeof(T));
        }
        size_ = new_size;
    }

    void reserve(size_t new_capacity) {
        if (new_capacity > capacity_) {
            grow(new_capacity);
        }
    }

    void push_back(const T& value) {
        if (size_ == capacity_) {
            grow(capacity_ == 0 ? 16 : capacity_ * 2);
        }
        data_[size_++] = value;
    }

    void push_back(T&& value) {
        if (size_ == capacity_) {
            grow(capacity_ == 0 ? 16 : capacity_ * 2);
        }
        data_[size_++] = std::move(value);
    }

    template <typename... Args>
    T& emplace_back(Args&&... args) {
        if (size_ == capacity_) {
            grow(capacity_ == 0 ? 16 : capacity_ * 2);
        }
        new (&data_[size_]) T(std::forward<Args>(args)...);
        return data_[size_++];
    }

    void pop_back() {
        assert(size_ > 0);
        --size_;
    }

    T& operator[](size_t i) {
        assert(i < size_);
        return data_[i];
    }

    const T& operator[](size_t i) const {
        assert(i < size_);
        return data_[i];
    }

    T& front() { return data_[0]; }
    const T& front() const { return data_[0]; }
    T& back() { return data_[size_ - 1]; }
    const T& back() const { return data_[size_ - 1]; }

    T* data() { return data_; }
    const T* data() const { return data_; }

    T* begin() { return data_; }
    const T* begin() const { return data_; }
    T* end() { return data_ + size_; }
    const T* end() const { return data_ + size_; }

    size_t size() const { return size_; }
    size_t capacity() const { return capacity_; }
    bool empty() const { return size_ == 0; }

    std::span<T> span() { return std::span<T>(data_, size_); }
    std::span<const T> span() const { return std::span<const T>(data_, size_); }

    // Move data out to std::vector (transfers ownership)
    std::vector<T> release() {
        std::vector<T> result;
        result.reserve(size_);
        for (size_t i = 0; i < size_; ++i) {
            result.push_back(std::move(data_[i]));
        }
        size_ = 0;
        return result;
    }

    ~ScratchBuffer() {
        std::free(data_);
    }

    // Move semantics
    ScratchBuffer(ScratchBuffer&& other) noexcept
        : data_(other.data_),
          size_(other.size_),
          capacity_(other.capacity_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    ScratchBuffer& operator=(ScratchBuffer&& other) noexcept {
        if (this != &other) {
            std::free(data_);
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    // No copying
    ScratchBuffer(const ScratchBuffer&) = delete;
    ScratchBuffer& operator=(const ScratchBuffer&) = delete;

private:
    void grow(size_t new_capacity) {
        T* new_data = static_cast<T*>(std::malloc(new_capacity * sizeof(T)));
        if (!new_data) {
            throw std::bad_alloc();
        }
        if (data_ && size_ > 0) {
            if constexpr (std::is_trivially_copyable_v<T>) {
                std::memcpy(new_data, data_, size_ * sizeof(T));
            } else {
                for (size_t i = 0; i < size_; ++i) {
                    new (&new_data[i]) T(std::move(data_[i]));
                    data_[i].~T();
                }
            }
        }
        std::free(data_);
        data_ = new_data;
        capacity_ = new_capacity;
    }

    T* data_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
};

// Thread-local scratch space for hot paths
// Provides reusable buffers to avoid repeated allocations
class ScratchSpace {
public:
    // Get thread-local instance
    static ScratchSpace& Local() {
        thread_local ScratchSpace instance;
        return instance;
    }

    // Arena for temporary allocations (reset after each operation)
    Arena& GetArena() { return arena_; }

    // Get a scratch buffer of type T (index selects which buffer)
    template <typename T, size_t Index = 0>
    ScratchBuffer<T>& GetBuffer() {
        static thread_local ScratchBuffer<T> buffer;
        buffer.clear();
        return buffer;
    }

    // Reset all scratch space
    void Reset() {
        arena_.Reset();
    }

private:
    ScratchSpace() = default;
    Arena arena_;
};

// RAII guard to reset scratch space on scope exit
class ScratchGuard {
public:
    explicit ScratchGuard(ScratchSpace& space) : space_(space) {}
    ~ScratchGuard() { space_.Reset(); }

    ScratchGuard(const ScratchGuard&) = delete;
    ScratchGuard& operator=(const ScratchGuard&) = delete;

private:
    ScratchSpace& space_;
};

}  // namespace llmap::core
