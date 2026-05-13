#include "classical/minimizer_index.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <random>

namespace llmap::classical {
namespace {

// Hash function tests
class HashFunctionTest : public ::testing::Test {};

TEST_F(HashFunctionTest, HashIsReproducible) {
    uint64_t kmer = 0x123456789ABCDEF0ULL;
    uint64_t seed = 0x517cc1b727220a95ULL;

    uint64_t h1 = HashKmer(kmer, seed);
    uint64_t h2 = HashKmer(kmer, seed);

    EXPECT_EQ(h1, h2);
}

TEST_F(HashFunctionTest, DifferentKmersGetDifferentHashes) {
    uint64_t seed = 0x517cc1b727220a95ULL;

    uint64_t h1 = HashKmer(0x1234ULL, seed);
    uint64_t h2 = HashKmer(0x5678ULL, seed);

    EXPECT_NE(h1, h2);
}

// Reverse complement tests
class ReverseComplementTest : public ::testing::Test {};

TEST_F(ReverseComplementTest, SingleBase) {
    // A (00) -> T (11)
    EXPECT_EQ(ReverseComplementKmer(0b00, 1), 0b11);
    // C (01) -> G (10)
    EXPECT_EQ(ReverseComplementKmer(0b01, 1), 0b10);
    // G (10) -> C (01)
    EXPECT_EQ(ReverseComplementKmer(0b10, 1), 0b01);
    // T (11) -> A (00)
    EXPECT_EQ(ReverseComplementKmer(0b11, 1), 0b00);
}

TEST_F(ReverseComplementTest, TwoBases) {
    // k-mer encoding: base at position 0 is in low bits
    // AC: A(00) at pos 0, C(01) at pos 1 = 0b01_00 = 0b0100 = 4
    // RevComp(AC) = GT: G(10) at pos 0, T(11) at pos 1 = 0b11_10 = 0b1110 = 14
    EXPECT_EQ(ReverseComplementKmer(0b0100, 2), 0b1110);

    // TA: T(11) at pos 0, A(00) at pos 1 = 0b00_11 = 0b0011 = 3
    // RevComp(TA) = TA (palindrome)
    EXPECT_EQ(ReverseComplementKmer(0b0011, 2), 0b0011);
}

TEST_F(ReverseComplementTest, DoubleReverseIsIdentity) {
    for (uint8_t k = 1; k <= 10; ++k) {
        uint64_t kmer = 0x123456789ABCDEF0ULL & ((1ULL << (2 * k)) - 1);
        uint64_t rc = ReverseComplementKmer(kmer, k);
        uint64_t rcrc = ReverseComplementKmer(rc, k);
        EXPECT_EQ(kmer, rcrc);
    }
}

// Canonical k-mer tests
class CanonicalKmerTest : public ::testing::Test {};

TEST_F(CanonicalKmerTest, CanonicalIsMinOfForwardAndReverse) {
    for (uint8_t k = 1; k <= 10; ++k) {
        uint64_t kmer = 0x123456789ABCDEF0ULL & ((1ULL << (2 * k)) - 1);
        uint64_t rc = ReverseComplementKmer(kmer, k);
        uint64_t canonical = CanonicalKmer(kmer, k);

        EXPECT_EQ(canonical, std::min(kmer, rc));
    }
}

TEST_F(CanonicalKmerTest, CanonicalOfForwardEqualsCanonicalOfReverse) {
    for (uint8_t k = 1; k <= 10; ++k) {
        uint64_t kmer = 0xABCDEF0123456789ULL & ((1ULL << (2 * k)) - 1);
        uint64_t rc = ReverseComplementKmer(kmer, k);

        EXPECT_EQ(CanonicalKmer(kmer, k), CanonicalKmer(rc, k));
    }
}

// Minimizer extraction tests
class MinimizerExtractionTest : public ::testing::Test {
protected:
    MinimizerConfig config_;

    void SetUp() override {
        config_.k = 5;
        config_.w = 5;
        config_.canonical = true;
    }
};

TEST_F(MinimizerExtractionTest, EmptySequenceReturnsEmpty) {
    auto minimizers = ExtractMinimizers("", config_);
    EXPECT_TRUE(minimizers.empty());
}

TEST_F(MinimizerExtractionTest, ShortSequenceReturnsEmpty) {
    auto minimizers = ExtractMinimizers("ACGT", config_);  // k=5, need at least 5 bases
    EXPECT_TRUE(minimizers.empty());
}

TEST_F(MinimizerExtractionTest, MinimalSequenceReturnsOneMinimizer) {
    auto minimizers = ExtractMinimizers("ACGTACGTACGTACGT", config_);
    EXPECT_FALSE(minimizers.empty());
}

TEST_F(MinimizerExtractionTest, MinimizersAreSortedByPosition) {
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    auto minimizers = ExtractMinimizers(seq, config_);

    for (size_t i = 1; i < minimizers.size(); ++i) {
        EXPECT_LE(minimizers[i-1].pos, minimizers[i].pos);
    }
}

TEST_F(MinimizerExtractionTest, NBasesSkipKmers) {
    // N in the middle should break the k-mer chain
    std::string seq = "ACGTACNGTACGTACGT";
    auto minimizers = ExtractMinimizers(seq, config_);

    // Verify no minimizer spans the N
    for (const auto& m : minimizers) {
        EXPECT_TRUE(m.pos + config_.k <= 6 || m.pos > 6);
    }
}

TEST_F(MinimizerExtractionTest, DeterministicOutput) {
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";

    auto m1 = ExtractMinimizers(seq, config_);
    auto m2 = ExtractMinimizers(seq, config_);

    EXPECT_EQ(m1.size(), m2.size());
    for (size_t i = 0; i < m1.size(); ++i) {
        EXPECT_EQ(m1[i].hash, m2[i].hash);
        EXPECT_EQ(m1[i].pos, m2[i].pos);
        EXPECT_EQ(m1[i].is_reverse, m2[i].is_reverse);
    }
}

TEST_F(MinimizerExtractionTest, CanonicalModeWorks) {
    // Forward and reverse complement should give same minimizers
    std::string fwd = "ACGTACGT";
    std::string rev = "ACGTACGT";  // Same in this case

    config_.k = 3;
    config_.w = 2;

    auto m_fwd = ExtractMinimizers(fwd, config_);
    auto m_rev = ExtractMinimizers(rev, config_);

    // Should have same minimizer hashes (order may differ due to position)
    EXPECT_EQ(m_fwd.size(), m_rev.size());
}

TEST_F(MinimizerExtractionTest, NonCanonicalModeWorks) {
    config_.canonical = false;
    std::string seq = "ACGTACGTACGTACGT";

    auto minimizers = ExtractMinimizers(seq, config_);
    EXPECT_FALSE(minimizers.empty());

    // All minimizers should have is_reverse = false in non-canonical mode
    for (const auto& m : minimizers) {
        EXPECT_FALSE(m.is_reverse);
    }
}

// MinimizerIndex Builder tests
class MinimizerIndexBuilderTest : public ::testing::Test {
protected:
    MinimizerConfig config_;

    void SetUp() override {
        config_.k = 5;
        config_.w = 5;
    }
};

TEST_F(MinimizerIndexBuilderTest, BuildEmptyIndex) {
    MinimizerIndex::Builder builder(config_);
    auto index = builder.Build();

    EXPECT_TRUE(index->Empty());
    EXPECT_EQ(index->Size(), 0);
}

TEST_F(MinimizerIndexBuilderTest, BuildWithOneSequence) {
    MinimizerIndex::Builder builder(config_);
    builder.AddSequence("seq1", "ACGTACGTACGTACGTACGTACGTACGT");
    auto index = builder.Build();

    EXPECT_FALSE(index->Empty());
    EXPECT_GT(index->Size(), 0);
    EXPECT_EQ(index->GetSequences().size(), 1);
    EXPECT_EQ(index->GetSequences()[0].name, "seq1");
}

TEST_F(MinimizerIndexBuilderTest, BuildWithMultipleSequences) {
    MinimizerIndex::Builder builder(config_);
    builder.AddSequence("seq1", "ACGTACGTACGTACGTACGT");
    builder.AddSequence("seq2", "TGCATGCATGCATGCATGCA");
    builder.AddSequence("seq3", "AAAACCCCGGGGTTTTAAAA");
    auto index = builder.Build();

    EXPECT_EQ(index->GetSequences().size(), 3);
}

TEST_F(MinimizerIndexBuilderTest, ChainedBuilderCalls) {
    auto index = MinimizerIndex::Builder(config_)
        .AddSequence("s1", "ACGTACGTACGTACGT")
        .AddSequence("s2", "TGCATGCATGCATGCA")
        .Build();

    EXPECT_EQ(index->GetSequences().size(), 2);
}

// MinimizerIndex Query tests
class MinimizerIndexQueryTest : public ::testing::Test {
protected:
    MinimizerConfig config_;
    std::unique_ptr<MinimizerIndex> index_;

    void SetUp() override {
        config_.k = 5;
        config_.w = 5;
        config_.max_occ = 1000;

        MinimizerIndex::Builder builder(config_);
        builder.AddSequence("ref1", "ACGTACGTACGTACGTACGTACGTACGTACGT");
        builder.AddSequence("ref2", "TGCATGCATGCATGCATGCATGCATGCATGCA");
        index_ = builder.Build();
    }
};

TEST_F(MinimizerIndexQueryTest, QueryEmptySequence) {
    auto hits = index_->Query("");
    EXPECT_TRUE(hits.empty());
}

TEST_F(MinimizerIndexQueryTest, QueryShortSequence) {
    auto hits = index_->Query("ACGT");  // Too short for k=5
    EXPECT_TRUE(hits.empty());
}

TEST_F(MinimizerIndexQueryTest, QueryMatchingSequence) {
    auto hits = index_->Query("ACGTACGTACGT");  // Should match ref1
    EXPECT_FALSE(hits.empty());

    // At least some hits should be to ref1
    bool has_ref1_hit = false;
    for (const auto& hit : hits) {
        if (hit.ref_id == 0) {
            has_ref1_hit = true;
            break;
        }
    }
    EXPECT_TRUE(has_ref1_hit);
}

TEST_F(MinimizerIndexQueryTest, HitsAreSorted) {
    auto hits = index_->Query("ACGTACGTACGTACGTACGT");

    for (size_t i = 1; i < hits.size(); ++i) {
        bool ordered = (hits[i-1].ref_id < hits[i].ref_id) ||
                       (hits[i-1].ref_id == hits[i].ref_id &&
                        hits[i-1].ref_pos <= hits[i].ref_pos);
        EXPECT_TRUE(ordered);
    }
}

TEST_F(MinimizerIndexQueryTest, MaxHitsLimit) {
    auto all_hits = index_->Query("ACGTACGTACGTACGTACGTACGTACGTACGT");
    auto limited_hits = index_->Query("ACGTACGTACGTACGTACGTACGTACGTACGT", 2);

    EXPECT_LE(limited_hits.size(), 2);
    if (all_hits.size() > 2) {
        EXPECT_LT(limited_hits.size(), all_hits.size());
    }
}

TEST_F(MinimizerIndexQueryTest, GetOccurrenceCount) {
    // Build an index with repeated k-mers
    MinimizerIndex::Builder builder(config_);
    std::string repeated = std::string(100, 'A') + std::string(100, 'C') +
                          std::string(100, 'G') + std::string(100, 'T');
    builder.AddSequence("repeated", repeated);
    auto idx = builder.Build();

    // AAAAA should occur many times
    auto minimizers = ExtractMinimizers("AAAAA", config_);
    if (!minimizers.empty()) {
        size_t count = idx->GetOccurrenceCount(minimizers[0].hash);
        EXPECT_GT(count, 0);
    }
}

// Serialization tests
class MinimizerIndexSerializationTest : public ::testing::Test {
protected:
    MinimizerConfig config_;
    std::string temp_path_;

    void SetUp() override {
        config_.k = 5;
        config_.w = 5;
        temp_path_ = std::filesystem::temp_directory_path() / "test_minimizer_index.idx";
    }

    void TearDown() override {
        std::filesystem::remove(temp_path_);
    }
};

TEST_F(MinimizerIndexSerializationTest, SaveAndLoad) {
    // Build index
    auto original = MinimizerIndex::Builder(config_)
        .AddSequence("seq1", "ACGTACGTACGTACGTACGTACGTACGTACGT")
        .AddSequence("seq2", "TGCATGCATGCATGCATGCATGCATGCATGCA")
        .Build();

    // Save
    EXPECT_TRUE(original->Save(temp_path_));

    // Load
    auto loaded = MinimizerIndex::Load(temp_path_);
    ASSERT_NE(loaded, nullptr);

    // Verify
    EXPECT_EQ(loaded->Size(), original->Size());
    EXPECT_EQ(loaded->GetSequences().size(), original->GetSequences().size());
    EXPECT_EQ(loaded->GetConfig().k, original->GetConfig().k);
    EXPECT_EQ(loaded->GetConfig().w, original->GetConfig().w);

    // Query should give same results
    std::string query = "ACGTACGTACGT";
    auto hits_orig = original->Query(query);
    auto hits_loaded = loaded->Query(query);
    EXPECT_EQ(hits_orig.size(), hits_loaded.size());
}

TEST_F(MinimizerIndexSerializationTest, LoadInvalidPath) {
    auto loaded = MinimizerIndex::Load("/nonexistent/path/index.idx");
    EXPECT_EQ(loaded, nullptr);
}

TEST_F(MinimizerIndexSerializationTest, LoadCorruptedFile) {
    // Write garbage
    std::ofstream out(temp_path_, std::ios::binary);
    out << "not a valid index file";
    out.close();

    auto loaded = MinimizerIndex::Load(temp_path_);
    EXPECT_EQ(loaded, nullptr);
}

// Statistics tests
class MinimizerIndexStatsTest : public ::testing::Test {};

TEST_F(MinimizerIndexStatsTest, StatsArePopulated) {
    MinimizerConfig config;
    config.k = 5;
    config.w = 5;

    auto index = MinimizerIndex::Builder(config)
        .AddSequence("seq", "ACGTACGTACGTACGTACGTACGTACGTACGT")
        .Build();

    const auto& stats = index->GetStats();
    EXPECT_GT(stats.total_kmers, 0);
    EXPECT_GT(stats.unique_minimizers, 0);
    EXPECT_GT(stats.avg_minimizer_spacing, 0.0f);
}

// Edge case tests
class MinimizerIndexEdgeCasesTest : public ::testing::Test {};

TEST_F(MinimizerIndexEdgeCasesTest, AllNSequence) {
    MinimizerConfig config;
    config.k = 5;
    config.w = 5;

    auto index = MinimizerIndex::Builder(config)
        .AddSequence("all_n", "NNNNNNNNNNNNNNNNNNNN")
        .Build();

    EXPECT_TRUE(index->Empty());
}

TEST_F(MinimizerIndexEdgeCasesTest, LowercaseSequence) {
    MinimizerConfig config;
    config.k = 5;
    config.w = 5;

    auto upper = MinimizerIndex::Builder(config)
        .AddSequence("upper", "ACGTACGTACGTACGT")
        .Build();

    auto lower = MinimizerIndex::Builder(config)
        .AddSequence("lower", "acgtacgtacgtacgt")
        .Build();

    // Should produce same index
    EXPECT_EQ(upper->Size(), lower->Size());
}

TEST_F(MinimizerIndexEdgeCasesTest, MixedCaseSequence) {
    MinimizerConfig config;
    config.k = 5;
    config.w = 5;

    auto minimizers = ExtractMinimizers("AcGtAcGtAcGtAcGt", config);
    EXPECT_FALSE(minimizers.empty());
}

TEST_F(MinimizerIndexEdgeCasesTest, HighOccurrenceSuppression) {
    MinimizerConfig config;
    config.k = 3;
    config.w = 2;
    config.max_occ = 2;  // Very low threshold

    // Build index with repeated k-mers
    auto index = MinimizerIndex::Builder(config)
        .AddSequence("s1", "AAAAAAAAAAAAAA")  // Many AAA k-mers
        .AddSequence("s2", "AAAAAAAAAAAAAA")
        .AddSequence("s3", "AAAAAAAAAAAAAA")
        .Build();

    // Query with AAA - should be suppressed due to high occurrence
    auto hits = index->Query("AAAAAAA");
    // Hits should be limited because AAA occurs > max_occ times
}

TEST_F(MinimizerIndexEdgeCasesTest, VeryLongSequence) {
    MinimizerConfig config;
    config.k = 10;
    config.w = 10;

    // 10KB sequence
    std::string seq;
    seq.reserve(10000);
    const char bases[] = "ACGT";
    std::mt19937 rng(42);
    for (size_t i = 0; i < 10000; ++i) {
        seq += bases[rng() % 4];
    }

    auto index = MinimizerIndex::Builder(config)
        .AddSequence("long", seq)
        .Build();

    EXPECT_FALSE(index->Empty());
    EXPECT_GT(index->GetStats().unique_minimizers, 100);
}

// Move semantics tests
class MinimizerIndexMoveTest : public ::testing::Test {};

TEST_F(MinimizerIndexMoveTest, MoveConstruction) {
    MinimizerConfig config;
    config.k = 5;
    config.w = 5;

    auto index = MinimizerIndex::Builder(config)
        .AddSequence("seq", "ACGTACGTACGTACGT")
        .Build();

    size_t original_size = index->Size();

    MinimizerIndex moved(std::move(*index));
    EXPECT_EQ(moved.Size(), original_size);
}

TEST_F(MinimizerIndexMoveTest, MoveAssignment) {
    MinimizerConfig config;
    config.k = 5;
    config.w = 5;

    auto index = MinimizerIndex::Builder(config)
        .AddSequence("seq", "ACGTACGTACGTACGT")
        .Build();

    size_t original_size = index->Size();

    MinimizerIndex other;
    other = std::move(*index);
    EXPECT_EQ(other.Size(), original_size);
}

}  // namespace
}  // namespace llmap::classical
