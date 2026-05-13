// LLmap — SIMD utilities for hot path optimization.
//
// Provides SIMD-accelerated operations for:
//   - Base encoding (ACGT → 2-bit)
//   - Validity checking (detect N/ambiguous bases)
//   - K-mer packing
//
// Uses runtime feature detection to select best implementation.
// Falls back to scalar code on unsupported platforms.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

// Feature detection - check for SSE4.2 and AVX2 at compile time
#if defined(__x86_64__) || defined(_M_X64)
    #if defined(__SSE4_2__)
        #define LLMAP_HAS_SSE42 1
    #endif
    #if defined(__AVX2__)
        #define LLMAP_HAS_AVX2 1
    #endif
#endif

#if defined(LLMAP_HAS_SSE42)
    #include <immintrin.h>
#endif

namespace llmap::core::simd {

// CPU feature flags (detected at runtime)
struct CpuFeatures {
    bool has_sse42 = false;
    bool has_avx2 = false;
    bool has_avx512 = false;

    static const CpuFeatures& Detect();
};

// Encode DNA bases to 2-bit representation.
// A=0, C=1, G=2, T=3. Invalid bases (N, etc.) are encoded as 0 but flagged.
//
// Parameters:
//   sequence: input DNA string
//   encoded: output array (must be at least sequence.size() bytes)
//   valid_mask: output bitmask, bit i is 1 if base i was valid ACGT
//               (must be at least (sequence.size() + 63) / 64 uint64_t values)
//
// Returns: number of bytes written to encoded
size_t EncodeBases(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask);

// Check if all bases in a range are valid (no N or ambiguous)
// Uses the valid_mask from EncodeBases
bool AllValid(const uint64_t* valid_mask, size_t start, size_t len, size_t total_len);

// Count valid bases in a range
size_t CountValid(const uint64_t* valid_mask, size_t start, size_t len, size_t total_len);

// Find next invalid base starting from position
// Returns total_len if all remaining bases are valid
size_t FindNextInvalid(const uint64_t* valid_mask, size_t start, size_t total_len);

// Pack 2-bit encoded bases into k-mers (both forward and reverse complement)
//
// Parameters:
//   encoded: 2-bit encoded bases from EncodeBases
//   len: number of bases
//   k: k-mer size (must be <= 32)
//   fwd_kmers: output forward k-mers (must be at least len - k + 1 values)
//   rev_kmers: output reverse-complement k-mers (same size)
//
// Returns: number of k-mers written
size_t PackKmers(
    const uint8_t* encoded,
    size_t len,
    size_t k,
    uint64_t* fwd_kmers,
    uint64_t* rev_kmers);

// Batch hash k-mers using the minimap2 invertible hash
//
// Parameters:
//   kmers: input k-mers
//   count: number of k-mers
//   seed: hash seed
//   hashes: output hashes (must be at least count values)
void HashKmersBatch(
    const uint64_t* kmers,
    size_t count,
    uint64_t seed,
    uint64_t* hashes);

namespace detail {

// Scalar implementations (always available)
size_t EncodeBasesScalar(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask);

void HashKmersBatchScalar(
    const uint64_t* kmers,
    size_t count,
    uint64_t seed,
    uint64_t* hashes);

#if defined(LLMAP_HAS_SSE42)
// SSE4.2 implementations
size_t EncodeBasesSse42(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask);
#endif

#if defined(LLMAP_HAS_AVX2)
// AVX2 implementations
size_t EncodeBasesAvx2(
    std::string_view sequence,
    uint8_t* encoded,
    uint64_t* valid_mask);

void HashKmersBatchAvx2(
    const uint64_t* kmers,
    size_t count,
    uint64_t seed,
    uint64_t* hashes);
#endif

}  // namespace detail

}  // namespace llmap::core::simd
