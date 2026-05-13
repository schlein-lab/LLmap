#include "core/simd.h"

namespace llmap::core::simd {

namespace detail {

void HashKmersBatchScalar(
    const uint64_t* kmers,
    size_t count,
    uint64_t seed,
    uint64_t* hashes)
{
    for (size_t i = 0; i < count; ++i) {
        uint64_t k = kmers[i];
        k = (~k + (k << 21)) ^ seed;
        k = k ^ (k >> 24);
        k = (k + (k << 3)) + (k << 8);
        k = k ^ (k >> 14);
        k = (k + (k << 2)) + (k << 4);
        k = k ^ (k >> 28);
        k = k + (k << 31);
        hashes[i] = k;
    }
}

#if defined(LLMAP_HAS_AVX2)

void HashKmersBatchAvx2(
    const uint64_t* kmers,
    size_t count,
    uint64_t seed,
    uint64_t* hashes)
{
    const __m256i vseed = _mm256_set1_epi64x(static_cast<int64_t>(seed));

    size_t i = 0;
    for (; i + 4 <= count; i += 4) {
        __m256i k = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(kmers + i));

        __m256i not_k = _mm256_xor_si256(k, _mm256_set1_epi64x(-1LL));
        __m256i k_shl21 = _mm256_slli_epi64(k, 21);
        k = _mm256_xor_si256(_mm256_add_epi64(not_k, k_shl21), vseed);

        k = _mm256_xor_si256(k, _mm256_srli_epi64(k, 24));

        __m256i k3 = _mm256_slli_epi64(k, 3);
        __m256i k8 = _mm256_slli_epi64(k, 8);
        k = _mm256_add_epi64(_mm256_add_epi64(k, k3), k8);

        k = _mm256_xor_si256(k, _mm256_srli_epi64(k, 14));

        __m256i k2 = _mm256_slli_epi64(k, 2);
        __m256i k4 = _mm256_slli_epi64(k, 4);
        k = _mm256_add_epi64(_mm256_add_epi64(k, k2), k4);

        k = _mm256_xor_si256(k, _mm256_srli_epi64(k, 28));

        k = _mm256_add_epi64(k, _mm256_slli_epi64(k, 31));

        _mm256_storeu_si256(reinterpret_cast<__m256i*>(hashes + i), k);
    }

    for (; i < count; ++i) {
        uint64_t k = kmers[i];
        k = (~k + (k << 21)) ^ seed;
        k = k ^ (k >> 24);
        k = (k + (k << 3)) + (k << 8);
        k = k ^ (k >> 14);
        k = (k + (k << 2)) + (k << 4);
        k = k ^ (k >> 28);
        k = k + (k << 31);
        hashes[i] = k;
    }
}

#endif  // LLMAP_HAS_AVX2

}  // namespace detail

void HashKmersBatch(
    const uint64_t* kmers,
    size_t count,
    uint64_t seed,
    uint64_t* hashes)
{
#if defined(LLMAP_HAS_AVX2)
    detail::HashKmersBatchAvx2(kmers, count, seed, hashes);
#else
    detail::HashKmersBatchScalar(kmers, count, seed, hashes);
#endif
}

}  // namespace llmap::core::simd
