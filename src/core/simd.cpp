#include "core/simd.h"
#include "core/simd_internal.h"

#include <bit>

namespace llmap::core::simd {

const CpuFeatures& CpuFeatures::Detect() {
    static CpuFeatures features = []() {
        CpuFeatures f;
#if defined(__x86_64__) || defined(_M_X64)
    #if defined(__SSE4_2__)
        f.has_sse42 = true;
    #endif
    #if defined(__AVX2__)
        f.has_avx2 = true;
    #endif
#endif
        return f;
    }();
    return features;
}

bool AllValid(const uint64_t* valid_mask, size_t start, size_t len, size_t total_len) {
    if (len == 0) return true;

    size_t end = start + len;
    size_t start_word = start / 64;
    size_t end_word = (end - 1) / 64;

    if (start_word == end_word) {
        uint64_t mask = ((1ULL << len) - 1) << (start % 64);
        return (valid_mask[start_word] & mask) == mask;
    }

    size_t first_bit = start % 64;
    if (first_bit != 0) {
        uint64_t mask = ~((1ULL << first_bit) - 1);
        if ((valid_mask[start_word] & mask) != mask) return false;
        ++start_word;
    }

    for (size_t w = start_word; w < end_word; ++w) {
        if (valid_mask[w] != ~0ULL) return false;
    }

    size_t last_bits = end % 64;
    if (last_bits != 0) {
        uint64_t mask = (1ULL << last_bits) - 1;
        if ((valid_mask[end_word] & mask) != mask) return false;
    }

    return true;
}

size_t CountValid(const uint64_t* valid_mask, size_t start, size_t len, size_t total_len) {
    if (len == 0) return 0;

    size_t count = 0;
    size_t end = start + len;
    size_t start_word = start / 64;
    size_t end_word = (end - 1) / 64;

    for (size_t w = start_word; w <= end_word; ++w) {
        uint64_t word = valid_mask[w];

        if (w == start_word && start % 64 != 0) {
            word &= ~((1ULL << (start % 64)) - 1);
        }
        if (w == end_word && end % 64 != 0) {
            word &= (1ULL << (end % 64)) - 1;
        }

        count += std::popcount(word);
    }

    return count;
}

size_t FindNextInvalid(const uint64_t* valid_mask, size_t start, size_t total_len) {
    size_t word = start / 64;
    size_t bit = start % 64;
    size_t num_words = (total_len + 63) / 64;

    if (word < num_words) {
        uint64_t w = valid_mask[word];
        w |= (1ULL << bit) - 1;
        if (word == num_words - 1 && total_len % 64 != 0) {
            w |= ~((1ULL << (total_len % 64)) - 1);
        }
        if (w != ~0ULL) {
            return word * 64 + std::countr_one(w);
        }
        ++word;
    }

    for (; word < num_words; ++word) {
        uint64_t w = valid_mask[word];
        if (word == num_words - 1 && total_len % 64 != 0) {
            w |= ~((1ULL << (total_len % 64)) - 1);
        }
        if (w != ~0ULL) {
            return word * 64 + std::countr_one(w);
        }
    }

    return total_len;
}

size_t PackKmers(
    const uint8_t* encoded,
    size_t len,
    size_t k,
    uint64_t* fwd_kmers,
    uint64_t* rev_kmers)
{
    if (len < k || k > 32) return 0;

    const size_t num_kmers = len - k + 1;
    const uint64_t mask = (1ULL << (2 * k)) - 1;

    uint64_t fwd = 0;
    uint64_t rev = 0;
    for (size_t i = 0; i < k; ++i) {
        uint8_t base = encoded[i];
        fwd = (fwd << 2) | base;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - base) << (2 * (k - 1)));
    }

    fwd_kmers[0] = fwd;
    rev_kmers[0] = rev;

    for (size_t i = 1; i < num_kmers; ++i) {
        uint8_t base = encoded[i + k - 1];
        fwd = ((fwd << 2) | base) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - base) << (2 * (k - 1)));
        fwd_kmers[i] = fwd;
        rev_kmers[i] = rev;
    }

    return num_kmers;
}

}  // namespace llmap::core::simd
