#include "core/simd.h"
#include "core/simd_internal.h"

#include <cstring>

namespace llmap::core::simd {

// Lookup tables defined here
alignas(64) const uint8_t kBaseEncode[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,  // A=65, C=67, G=71
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // T=84
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,  // a=97, c=99, g=103
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // t=116
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

alignas(64) const uint8_t kBaseValid[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  // A=65, C=67, G=71
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // T=84
    0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  // a=97, c=99, g=103
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // t=116
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

namespace detail {

size_t EncodeBasesScalar(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask)
{
    const size_t len = sequence.size();
    const size_t mask_words = (len + 63) / 64;

    std::memset(valid_mask, 0, mask_words * sizeof(uint64_t));

    for (size_t i = 0; i < len; ++i) {
        uint8_t c = static_cast<uint8_t>(sequence[i]);
        encoded[i] = kBaseEncode[c];
        if (kBaseValid[c]) {
            valid_mask[i / 64] |= (1ULL << (i % 64));
        }
    }

    return len;
}

#if defined(LLMAP_HAS_SSE42)

size_t EncodeBasesSse42(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask)
{
    const size_t len = sequence.size();
    const size_t mask_words = (len + 63) / 64;
    const char* data = sequence.data();

    std::memset(valid_mask, 0, mask_words * sizeof(uint64_t));

    const __m128i encode_lo = _mm_setr_epi8(
        0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m128i valid_lo = _mm_setr_epi8(
        0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m128i mask_0f = _mm_set1_epi8(0x0F);
    const __m128i mask_20 = _mm_set1_epi8(0x20);

    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        __m128i chars = _mm_loadu_si128(reinterpret_cast<const __m128i*>(data + i));
        __m128i lower = _mm_or_si128(chars, mask_20);
        __m128i nibbles = _mm_and_si128(lower, mask_0f);

        __m128i enc = _mm_shuffle_epi8(encode_lo, nibbles);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(encoded + i), enc);

        __m128i val = _mm_shuffle_epi8(valid_lo, nibbles);
        uint32_t valid_bits = static_cast<uint32_t>(_mm_movemask_epi8(
            _mm_cmpgt_epi8(val, _mm_setzero_si128())));

        size_t word = i / 64;
        size_t bit = i % 64;
        if (bit + 16 <= 64) {
            valid_mask[word] |= static_cast<uint64_t>(valid_bits) << bit;
        } else {
            valid_mask[word] |= static_cast<uint64_t>(valid_bits) << bit;
            valid_mask[word + 1] |= static_cast<uint64_t>(valid_bits) >> (64 - bit);
        }
    }

    for (; i < len; ++i) {
        uint8_t c = static_cast<uint8_t>(data[i]);
        encoded[i] = kBaseEncode[c];
        if (kBaseValid[c]) {
            valid_mask[i / 64] |= (1ULL << (i % 64));
        }
    }

    return len;
}

#endif  // LLMAP_HAS_SSE42

#if defined(LLMAP_HAS_AVX2)

size_t EncodeBasesAvx2(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask)
{
    const size_t len = sequence.size();
    const size_t mask_words = (len + 63) / 64;
    const char* data = sequence.data();

    std::memset(valid_mask, 0, mask_words * sizeof(uint64_t));

    const __m256i encode_lo = _mm256_setr_epi8(
        0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m256i valid_lo = _mm256_setr_epi8(
        0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m256i mask_0f = _mm256_set1_epi8(0x0F);
    const __m256i mask_20 = _mm256_set1_epi8(0x20);

    size_t i = 0;
    for (; i + 32 <= len; i += 32) {
        __m256i chars = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + i));
        __m256i lower = _mm256_or_si256(chars, mask_20);
        __m256i nibbles = _mm256_and_si256(lower, mask_0f);

        __m256i enc = _mm256_shuffle_epi8(encode_lo, nibbles);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(encoded + i), enc);

        __m256i val = _mm256_shuffle_epi8(valid_lo, nibbles);
        uint32_t valid_bits = static_cast<uint32_t>(_mm256_movemask_epi8(
            _mm256_cmpgt_epi8(val, _mm256_setzero_si256())));

        size_t word = i / 64;
        size_t bit = i % 64;
        if (bit + 32 <= 64) {
            valid_mask[word] |= static_cast<uint64_t>(valid_bits) << bit;
        } else {
            valid_mask[word] |= static_cast<uint64_t>(valid_bits) << bit;
            valid_mask[word + 1] |= static_cast<uint64_t>(valid_bits) >> (64 - bit);
        }
    }

#if defined(LLMAP_HAS_SSE42)
    const __m128i encode_lo_sse = _mm_setr_epi8(
        0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m128i valid_lo_sse = _mm_setr_epi8(
        0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m128i mask_0f_sse = _mm_set1_epi8(0x0F);
    const __m128i mask_20_sse = _mm_set1_epi8(0x20);

    for (; i + 16 <= len; i += 16) {
        __m128i chars = _mm_loadu_si128(reinterpret_cast<const __m128i*>(data + i));
        __m128i lower = _mm_or_si128(chars, mask_20_sse);
        __m128i nibbles = _mm_and_si128(lower, mask_0f_sse);
        __m128i enc = _mm_shuffle_epi8(encode_lo_sse, nibbles);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(encoded + i), enc);

        __m128i val = _mm_shuffle_epi8(valid_lo_sse, nibbles);
        uint32_t valid_bits = static_cast<uint32_t>(_mm_movemask_epi8(
            _mm_cmpgt_epi8(val, _mm_setzero_si128())));

        size_t word = i / 64;
        size_t bit = i % 64;
        if (bit + 16 <= 64) {
            valid_mask[word] |= static_cast<uint64_t>(valid_bits) << bit;
        } else {
            valid_mask[word] |= static_cast<uint64_t>(valid_bits) << bit;
            valid_mask[word + 1] |= static_cast<uint64_t>(valid_bits) >> (64 - bit);
        }
    }
#endif

    for (; i < len; ++i) {
        uint8_t c = static_cast<uint8_t>(data[i]);
        encoded[i] = kBaseEncode[c];
        if (kBaseValid[c]) {
            valid_mask[i / 64] |= (1ULL << (i % 64));
        }
    }

    return len;
}

#endif  // LLMAP_HAS_AVX2

}  // namespace detail

size_t EncodeBases(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask)
{
#if defined(LLMAP_HAS_AVX2)
    return detail::EncodeBasesAvx2(sequence, encoded, valid_mask);
#elif defined(LLMAP_HAS_SSE42)
    return detail::EncodeBasesSse42(sequence, encoded, valid_mask);
#else
    return detail::EncodeBasesScalar(sequence, encoded, valid_mask);
#endif
}

}  // namespace llmap::core::simd
