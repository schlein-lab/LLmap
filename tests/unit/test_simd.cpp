#include <gtest/gtest.h>

#include "core/simd.h"

#include <array>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace llmap::core::simd;

class SimdTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Seed for reproducible tests
        rng_.seed(42);
    }

    std::string RandomDna(size_t len) {
        static const char bases[] = "ACGT";
        std::string seq;
        seq.reserve(len);
        for (size_t i = 0; i < len; ++i) {
            seq += bases[rng_() % 4];
        }
        return seq;
    }

    std::string RandomDnaWithN(size_t len, double n_rate = 0.05) {
        static const char bases[] = "ACGTN";
        std::string seq;
        seq.reserve(len);
        for (size_t i = 0; i < len; ++i) {
            if (std::uniform_real_distribution<>(0.0, 1.0)(rng_) < n_rate) {
                seq += 'N';
            } else {
                seq += bases[rng_() % 4];
            }
        }
        return seq;
    }

    std::mt19937 rng_;
};

TEST_F(SimdTest, CpuFeaturesDetect) {
    const auto& features = CpuFeatures::Detect();
    // Just check it doesn't crash - actual values depend on CPU
    (void)features.has_sse42;
    (void)features.has_avx2;
    (void)features.has_avx512;
}

TEST_F(SimdTest, EncodeBasesSimple) {
    std::string seq = "ACGT";
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask(1);

    size_t n = EncodeBases(seq, encoded.data(), valid_mask.data());

    EXPECT_EQ(n, 4);
    EXPECT_EQ(encoded[0], 0);  // A
    EXPECT_EQ(encoded[1], 1);  // C
    EXPECT_EQ(encoded[2], 2);  // G
    EXPECT_EQ(encoded[3], 3);  // T
    EXPECT_EQ(valid_mask[0] & 0xF, 0xF);  // All 4 bases valid
}

TEST_F(SimdTest, EncodeBasesLowercase) {
    std::string seq = "acgt";
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask(1);

    size_t n = EncodeBases(seq, encoded.data(), valid_mask.data());

    EXPECT_EQ(n, 4);
    EXPECT_EQ(encoded[0], 0);  // a
    EXPECT_EQ(encoded[1], 1);  // c
    EXPECT_EQ(encoded[2], 2);  // g
    EXPECT_EQ(encoded[3], 3);  // t
    EXPECT_EQ(valid_mask[0] & 0xF, 0xF);
}

TEST_F(SimdTest, EncodeBasesWithN) {
    std::string seq = "ACNGT";
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask(1);

    size_t n = EncodeBases(seq, encoded.data(), valid_mask.data());

    EXPECT_EQ(n, 5);
    EXPECT_EQ(encoded[0], 0);  // A
    EXPECT_EQ(encoded[1], 1);  // C
    // N encodes to 0 but is invalid
    EXPECT_EQ(encoded[3], 2);  // G
    EXPECT_EQ(encoded[4], 3);  // T

    // Bits 0,1,3,4 should be set (A,C,G,T valid), bit 2 (N) should be 0
    EXPECT_EQ(valid_mask[0] & 0x1F, 0x1B);  // 11011 binary
}

TEST_F(SimdTest, EncodeBasesLongSequence) {
    std::string seq = RandomDna(1000);
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask((seq.size() + 63) / 64);

    size_t n = EncodeBases(seq, encoded.data(), valid_mask.data());

    EXPECT_EQ(n, seq.size());

    // Verify each base
    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t expected;
        switch (seq[i]) {
            case 'A': expected = 0; break;
            case 'C': expected = 1; break;
            case 'G': expected = 2; break;
            case 'T': expected = 3; break;
            default: expected = 0; break;
        }
        EXPECT_EQ(encoded[i], expected) << "Mismatch at position " << i;
    }
}

TEST_F(SimdTest, EncodeBasesLongWithN) {
    std::string seq = RandomDnaWithN(500, 0.1);
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask((seq.size() + 63) / 64);

    size_t n = EncodeBases(seq, encoded.data(), valid_mask.data());

    EXPECT_EQ(n, seq.size());

    // Verify validity mask
    for (size_t i = 0; i < seq.size(); ++i) {
        bool expected_valid = (seq[i] == 'A' || seq[i] == 'C' ||
                               seq[i] == 'G' || seq[i] == 'T');
        bool actual_valid = (valid_mask[i / 64] & (1ULL << (i % 64))) != 0;
        EXPECT_EQ(actual_valid, expected_valid)
            << "Validity mismatch at position " << i << " (char=" << seq[i] << ")";
    }
}

TEST_F(SimdTest, AllValidEmpty) {
    uint64_t mask = 0;
    EXPECT_TRUE(AllValid(&mask, 0, 0, 0));
}

TEST_F(SimdTest, AllValidSimple) {
    uint64_t mask = 0xF;  // First 4 bits set
    EXPECT_TRUE(AllValid(&mask, 0, 4, 4));
    EXPECT_FALSE(AllValid(&mask, 0, 5, 8));  // 5th bit not set
}

TEST_F(SimdTest, AllValidLong) {
    std::vector<uint64_t> mask(2, ~0ULL);  // All bits set
    EXPECT_TRUE(AllValid(mask.data(), 0, 100, 128));
    EXPECT_TRUE(AllValid(mask.data(), 50, 50, 128));

    mask[1] &= ~(1ULL << 20);  // Clear bit 84
    EXPECT_FALSE(AllValid(mask.data(), 80, 10, 128));
    EXPECT_TRUE(AllValid(mask.data(), 80, 4, 128));
}

TEST_F(SimdTest, CountValidSimple) {
    uint64_t mask = 0b11011;  // 4 valid
    EXPECT_EQ(CountValid(&mask, 0, 5, 5), 4);
}

TEST_F(SimdTest, CountValidLong) {
    std::string seq = RandomDnaWithN(200, 0.1);
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask((seq.size() + 63) / 64);
    EncodeBases(seq, encoded.data(), valid_mask.data());

    size_t expected = 0;
    for (char c : seq) {
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') ++expected;
    }

    EXPECT_EQ(CountValid(valid_mask.data(), 0, seq.size(), seq.size()), expected);
}

TEST_F(SimdTest, FindNextInvalidNoInvalid) {
    uint64_t mask = ~0ULL;
    EXPECT_EQ(FindNextInvalid(&mask, 0, 64), 64);
    EXPECT_EQ(FindNextInvalid(&mask, 32, 64), 64);
}

TEST_F(SimdTest, FindNextInvalidSimple) {
    uint64_t mask = 0b111011;  // Bit 2 is 0
    EXPECT_EQ(FindNextInvalid(&mask, 0, 6), 2);
    EXPECT_EQ(FindNextInvalid(&mask, 3, 6), 6);
}

TEST_F(SimdTest, FindNextInvalidLong) {
    std::vector<uint64_t> mask(3, ~0ULL);
    mask[1] &= ~(1ULL << 30);  // Clear bit 94

    EXPECT_EQ(FindNextInvalid(mask.data(), 0, 192), 94);
    EXPECT_EQ(FindNextInvalid(mask.data(), 95, 192), 192);
}

TEST_F(SimdTest, PackKmersSimple) {
    // ACGT -> 0b00011011 = 0x1B
    uint8_t encoded[] = {0, 1, 2, 3};  // A C G T
    uint64_t fwd[1], rev[1];

    size_t n = PackKmers(encoded, 4, 4, fwd, rev);

    EXPECT_EQ(n, 1);
    EXPECT_EQ(fwd[0], 0b00011011);  // ACGT
    // ACGT is a palindrome! complement(ACGT) = TGCA, reversed = ACGT
    EXPECT_EQ(rev[0], 0b00011011);  // Same as forward for palindrome
}

TEST_F(SimdTest, PackKmersNonPalindrome) {
    // AACC -> 0b00000101
    uint8_t encoded[] = {0, 0, 1, 1};  // A A C C
    uint64_t fwd[1], rev[1];

    size_t n = PackKmers(encoded, 4, 4, fwd, rev);

    EXPECT_EQ(n, 1);
    EXPECT_EQ(fwd[0], 0b00000101);  // AACC
    // Reverse complement: AACC -> complement TTGG -> reverse GGTT
    // G=2, G=2, T=3, T=3 -> 0b10101111
    EXPECT_EQ(rev[0], 0b10101111);
}

TEST_F(SimdTest, PackKmersSliding) {
    // ACGTACGT
    uint8_t encoded[] = {0, 1, 2, 3, 0, 1, 2, 3};
    uint64_t fwd[5], rev[5];

    size_t n = PackKmers(encoded, 8, 4, fwd, rev);

    EXPECT_EQ(n, 5);  // 8 - 4 + 1 = 5 k-mers
}

TEST_F(SimdTest, PackKmersTooShort) {
    uint8_t encoded[] = {0, 1, 2};
    uint64_t fwd[1], rev[1];

    size_t n = PackKmers(encoded, 3, 4, fwd, rev);
    EXPECT_EQ(n, 0);
}

TEST_F(SimdTest, HashKmersBatchSimple) {
    uint64_t kmers[] = {0x1B, 0x2C, 0x3D, 0x4E};
    uint64_t hashes[4];

    HashKmersBatch(kmers, 4, 0x517cc1b727220a95ULL, hashes);

    // Verify hashes are different
    EXPECT_NE(hashes[0], hashes[1]);
    EXPECT_NE(hashes[1], hashes[2]);
    EXPECT_NE(hashes[2], hashes[3]);
}

TEST_F(SimdTest, HashKmersBatchMatchesScalar) {
    uint64_t kmers[100];
    for (int i = 0; i < 100; ++i) {
        kmers[i] = rng_();
    }

    uint64_t hashes_simd[100];
    uint64_t hashes_scalar[100];

    HashKmersBatch(kmers, 100, 0x517cc1b727220a95ULL, hashes_simd);
    detail::HashKmersBatchScalar(kmers, 100, 0x517cc1b727220a95ULL, hashes_scalar);

    for (int i = 0; i < 100; ++i) {
        EXPECT_EQ(hashes_simd[i], hashes_scalar[i]) << "Mismatch at index " << i;
    }
}

TEST_F(SimdTest, EncodingMatchesScalar) {
    std::string seq = RandomDnaWithN(500, 0.1);

    std::vector<uint8_t> encoded_simd(seq.size());
    std::vector<uint64_t> valid_simd((seq.size() + 63) / 64);

    std::vector<uint8_t> encoded_scalar(seq.size());
    std::vector<uint64_t> valid_scalar((seq.size() + 63) / 64);

    EncodeBases(seq, encoded_simd.data(), valid_simd.data());
    detail::EncodeBasesScalar(seq, encoded_scalar.data(), valid_scalar.data());

    for (size_t i = 0; i < seq.size(); ++i) {
        EXPECT_EQ(encoded_simd[i], encoded_scalar[i])
            << "Encoding mismatch at position " << i;
    }

    for (size_t i = 0; i < valid_simd.size(); ++i) {
        EXPECT_EQ(valid_simd[i], valid_scalar[i])
            << "Validity mismatch at word " << i;
    }
}

TEST_F(SimdTest, EncodeAllBases) {
    std::string seq = "AaCcGgTt";
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask(1);

    EncodeBases(seq, encoded.data(), valid_mask.data());

    EXPECT_EQ(encoded[0], 0);  // A
    EXPECT_EQ(encoded[1], 0);  // a
    EXPECT_EQ(encoded[2], 1);  // C
    EXPECT_EQ(encoded[3], 1);  // c
    EXPECT_EQ(encoded[4], 2);  // G
    EXPECT_EQ(encoded[5], 2);  // g
    EXPECT_EQ(encoded[6], 3);  // T
    EXPECT_EQ(encoded[7], 3);  // t
    EXPECT_EQ(valid_mask[0] & 0xFF, 0xFF);  // All valid
}

TEST_F(SimdTest, EncodeAmbiguousBases) {
    std::string seq = "ACGTRYSWKMBDHVN";
    std::vector<uint8_t> encoded(seq.size());
    std::vector<uint64_t> valid_mask(1);

    EncodeBases(seq, encoded.data(), valid_mask.data());

    // Only A, C, G, T should be valid (bits 0-3)
    EXPECT_EQ(valid_mask[0] & 0xFFFF, 0x000F);
}

TEST_F(SimdTest, LargeBatchHashConsistency) {
    // Test that batched hashing produces consistent results
    std::vector<uint64_t> kmers(1000);
    for (size_t i = 0; i < kmers.size(); ++i) {
        kmers[i] = rng_();
    }

    std::vector<uint64_t> hashes1(1000);
    std::vector<uint64_t> hashes2(1000);

    HashKmersBatch(kmers.data(), 1000, 0x123456789ABCDEFULL, hashes1.data());
    HashKmersBatch(kmers.data(), 1000, 0x123456789ABCDEFULL, hashes2.data());

    for (size_t i = 0; i < 1000; ++i) {
        EXPECT_EQ(hashes1[i], hashes2[i]);
    }
}

TEST_F(SimdTest, EncodeEdgeLengths) {
    // Test various lengths that hit different SIMD code paths
    for (size_t len : {1, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129}) {
        std::string seq = RandomDna(len);
        std::vector<uint8_t> encoded(seq.size());
        std::vector<uint64_t> valid_mask((seq.size() + 63) / 64);

        size_t n = EncodeBases(seq, encoded.data(), valid_mask.data());
        EXPECT_EQ(n, len) << "Length " << len;

        // Verify all valid (no N in RandomDna)
        size_t valid_count = CountValid(valid_mask.data(), 0, len, len);
        EXPECT_EQ(valid_count, len) << "Length " << len;
    }
}
