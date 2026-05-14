#include <gtest/gtest.h>

#include <random>
#include <set>

#include "classical/classical_pipeline.h"

using namespace llmap::classical;

class ClassicalPipelineTest : public ::testing::Test {
protected:
    ClassicalPipelineConfig config;

    void SetUp() override {
        config = ClassicalPipelineConfig{};
        config.minimizer_config.k = 11;
        config.minimizer_config.w = 5;
        config.chain_config.min_chain_anchors = 2;
        config.min_identity = 0.5f;
        config.min_aligned_bases = 20;
    }

    std::string MakeSequence(size_t len, char base = 'A') {
        return std::string(len, base);
    }

    std::string MakeRandomSequence(size_t len, unsigned seed = 42) {
        static const char bases[] = "ACGT";
        std::string seq;
        seq.reserve(len);
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> dist(0, 3);
        for (size_t i = 0; i < len; ++i) {
            seq.push_back(bases[dist(rng)]);
        }
        return seq;
    }
};

TEST_F(ClassicalPipelineTest, DefaultConstruction) {
    ClassicalPipeline pipeline;
    EXPECT_FALSE(pipeline.HasIndex());
}

TEST_F(ClassicalPipelineTest, ConstructWithConfig) {
    ClassicalPipeline pipeline(config);
    EXPECT_FALSE(pipeline.HasIndex());
    EXPECT_EQ(pipeline.Config().minimizer_config.k, 11);
}

TEST_F(ClassicalPipelineTest, SetOwnedIndex) {
    ClassicalPipeline pipeline(config);

    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", MakeRandomSequence(500))
        .Build();

    pipeline.SetIndex(std::move(index));
    EXPECT_TRUE(pipeline.HasIndex());
}

TEST_F(ClassicalPipelineTest, SetNonOwnedIndex) {
    ClassicalPipeline pipeline(config);

    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", MakeRandomSequence(500))
        .Build();

    pipeline.SetIndex(index.get());
    EXPECT_TRUE(pipeline.HasIndex());
}

TEST_F(ClassicalPipelineTest, AlignEmptyReadReturnsNoAlignment) {
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", MakeRandomSequence(500))
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", "");
    EXPECT_FALSE(result.HasAlignment());
    EXPECT_EQ(result.query_name, "read1");
}

TEST_F(ClassicalPipelineTest, AlignWithoutIndexReturnsNoAlignment) {
    ClassicalPipeline pipeline(config);

    auto result = pipeline.AlignRead("read1", MakeRandomSequence(100));
    EXPECT_FALSE(result.HasAlignment());
}

TEST_F(ClassicalPipelineTest, AlignExactMatchProducesAlignment) {
    std::string ref_seq = MakeRandomSequence(500);
    std::string query_seq = ref_seq.substr(100, 200);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query_seq);

    EXPECT_TRUE(result.HasAlignment());
    EXPECT_GT(result.num_hits, 0);
    EXPECT_GT(result.num_chains, 0);
}

TEST_F(ClassicalPipelineTest, AlignmentHasCorrectMetadata) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query_seq = ref_seq.substr(100, 200);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("test_ref", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("test_read", query_seq);

    ASSERT_TRUE(result.HasAlignment());
    const auto* primary = result.PrimaryAlignment();
    ASSERT_NE(primary, nullptr);

    EXPECT_EQ(primary->query_name, "test_read");
    EXPECT_EQ(primary->ref_name, "test_ref");
    EXPECT_TRUE(primary->is_primary);
    EXPECT_GT(primary->mapq, 0);
}

TEST_F(ClassicalPipelineTest, UnrelatedSequenceProducesNoAlignment) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query_seq = MakeRandomSequence(200, 99);  // Different seed

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query_seq);

    // May or may not have hits depending on random sequences
    // but if alignment exists, it should have low identity
    if (result.HasAlignment()) {
        EXPECT_LT(result.alignments[0].identity, 0.9f);
    }
}

TEST_F(ClassicalPipelineTest, AlignBatchMultipleReads) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names = {"read1", "read2", "read3"};
    std::vector<std::string> query_seqs = {
        ref_seq.substr(100, 150),
        ref_seq.substr(300, 150),
        ref_seq.substr(500, 150)
    };

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto results = pipeline.AlignReads(query_names, query_seqs);

    EXPECT_EQ(results.size(), 3);
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_EQ(results[i].query_name, query_names[i]);
    }
}

TEST_F(ClassicalPipelineTest, StatsUpdatedAfterBatch) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names = {"read1", "read2"};
    std::vector<std::string> query_seqs = {
        ref_seq.substr(100, 150),
        MakeRandomSequence(150, 99)  // Unlikely to align
    };

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    pipeline.AlignReads(query_names, query_seqs);
    const auto& stats = pipeline.Stats();

    EXPECT_EQ(stats.total_reads, 2);
    EXPECT_GT(stats.total_time_ms, 0.0f);
}

TEST_F(ClassicalPipelineTest, CigarStringGeneration) {
    ClassicalAlignment aln;
    aln.cigar = {
        {CigarOp::Match, 50},
        {CigarOp::Insertion, 3},
        {CigarOp::Match, 30},
        {CigarOp::Deletion, 2},
        {CigarOp::Match, 20}
    };

    std::string cigar_str = aln.CigarString();
    EXPECT_FALSE(cigar_str.empty());
    EXPECT_NE(cigar_str.find("50M"), std::string::npos);
    EXPECT_NE(cigar_str.find("3I"), std::string::npos);
}

TEST_F(ClassicalPipelineTest, AlignedBasesCalculation) {
    ClassicalAlignment aln;
    aln.cigar = {
        {CigarOp::Match, 50},
        {CigarOp::Insertion, 3},
        {CigarOp::Match, 30}
    };

    EXPECT_EQ(aln.AlignedBases(), 80);  // 50 + 30
}

TEST_F(ClassicalPipelineTest, MultipleReferencesAlignToCorrectOne) {
    std::string ref1 = MakeRandomSequence(500, 1);
    std::string ref2 = MakeRandomSequence(500, 2);
    std::string query = ref2.substr(100, 150);  // From ref2

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref1)
        .AddSequence("ref2", ref2)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    if (result.HasAlignment()) {
        // Should align to ref2
        EXPECT_EQ(result.alignments[0].ref_name, "ref2");
    }
}

TEST_F(ClassicalPipelineTest, SecondaryAlignmentsReported) {
    std::string ref_seq = MakeRandomSequence(1000, 1);
    // Duplicate a region to create secondary alignment possibility
    std::string duplicated = ref_seq;
    duplicated.replace(500, 200, ref_seq.substr(100, 200));

    std::string query = ref_seq.substr(100, 150);

    config.report_secondary = true;
    config.max_alignments = 5;
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", duplicated)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    // May have multiple alignments due to duplication
    if (result.alignments.size() > 1) {
        EXPECT_TRUE(result.alignments[0].is_primary);
        EXPECT_FALSE(result.alignments[1].is_primary);
    }
}

TEST_F(ClassicalPipelineTest, ConvenienceFunctionWorksEndToEnd) {
    std::vector<std::string> ref_names = {"ref1"};
    std::vector<std::string> ref_seqs = {MakeRandomSequence(500, 1)};
    std::vector<std::string> query_names = {"read1"};
    std::vector<std::string> query_seqs = {ref_seqs[0].substr(100, 150)};

    auto results = AlignWithClassicalPath(
        ref_names, ref_seqs, query_names, query_seqs, config);

    EXPECT_EQ(results.size(), 1);
    EXPECT_TRUE(results[0].HasAlignment());
}

TEST_F(ClassicalPipelineTest, TimingStatsAreReasonable) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    EXPECT_GE(result.seeding_time_us, 0.0f);
    EXPECT_GE(result.chaining_time_us, 0.0f);
    EXPECT_GE(result.extension_time_us, 0.0f);
}

TEST_F(ClassicalPipelineTest, MoveConstructor) {
    ClassicalPipeline pipeline1(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", MakeRandomSequence(500))
        .Build();
    pipeline1.SetIndex(std::move(index));

    ClassicalPipeline pipeline2(std::move(pipeline1));
    EXPECT_TRUE(pipeline2.HasIndex());
}

TEST_F(ClassicalPipelineTest, MoveAssignment) {
    ClassicalPipeline pipeline1(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", MakeRandomSequence(500))
        .Build();
    pipeline1.SetIndex(std::move(index));

    ClassicalPipeline pipeline2;
    pipeline2 = std::move(pipeline1);
    EXPECT_TRUE(pipeline2.HasIndex());
}

TEST_F(ClassicalPipelineTest, ReadAlignmentResultPrimaryAccessor) {
    ReadAlignmentResult result;
    EXPECT_EQ(result.PrimaryAlignment(), nullptr);

    ClassicalAlignment aln1;
    aln1.is_primary = false;
    result.alignments.push_back(aln1);

    EXPECT_EQ(result.PrimaryAlignment(), nullptr);

    ClassicalAlignment aln2;
    aln2.is_primary = true;
    result.alignments.push_back(aln2);

    EXPECT_NE(result.PrimaryAlignment(), nullptr);
    EXPECT_TRUE(result.PrimaryAlignment()->is_primary);
}

TEST_F(ClassicalPipelineTest, ChainMetadataPreservedInAlignment) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 200);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    if (result.HasAlignment()) {
        const auto& aln = result.alignments[0];
        EXPECT_GT(aln.chain_anchors, 0);
        EXPECT_GT(aln.chain_score, 0);
    }
}

TEST_F(ClassicalPipelineTest, MaxAlignmentsRespected) {
    std::string ref_seq = MakeRandomSequence(1000, 1);
    // Create repeated regions
    std::string repeated = ref_seq;
    for (int i = 0; i < 5; ++i) {
        repeated += ref_seq.substr(0, 200);
    }

    std::string query = ref_seq.substr(0, 150);

    config.max_alignments = 2;
    config.report_secondary = true;
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", repeated)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    EXPECT_LE(result.alignments.size(), 2);
}

TEST_F(ClassicalPipelineTest, MinIdentityFilterApplied) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = MakeRandomSequence(150, 99);  // Unrelated

    config.min_identity = 0.99f;  // Very high threshold
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    // Should have no alignments passing identity filter
    // (may have hits but alignments filtered out)
    for (const auto& aln : result.alignments) {
        EXPECT_GE(aln.identity, config.min_identity);
    }
}

// === Chain End Extension Tests (Phase A.3) ===

TEST_F(ClassicalPipelineTest, EndExtensionCoversFullQuery) {
    // Query is an exact substring from ref, but anchors don't cover edges
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(50, 200);  // Exact match at offset 50

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    // Set reference sequences for WFA2 extension
    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    ASSERT_TRUE(result.HasAlignment());
    const auto& aln = result.alignments[0];

    // With end extension, query_start should be 0 and query_end should be query length
    EXPECT_EQ(aln.query_start, 0);
    EXPECT_EQ(aln.query_end, static_cast<int32_t>(query.size()));
}

TEST_F(ClassicalPipelineTest, EndExtensionWithRefSeqsProducesAccurateCigar) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    // Take substring that starts a bit into the sequence
    std::string query = ref_seq.substr(100, 150);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    // Set reference sequences for WFA2 extension
    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    ASSERT_TRUE(result.HasAlignment());
    const auto& aln = result.alignments[0];

    // Alignment should have good identity for exact match
    EXPECT_GT(aln.identity, 0.8f);

    // CIGAR should not be empty
    EXPECT_FALSE(aln.cigar.empty());
}

TEST_F(ClassicalPipelineTest, EndExtensionSetsQueryStartToZero) {
    // Verify that query_start is always 0 after end extension
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 130);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    if (result.HasAlignment()) {
        const auto& aln = result.alignments[0];
        // With end extension, query_start should always be 0
        EXPECT_EQ(aln.query_start, 0);
    }
}

TEST_F(ClassicalPipelineTest, EndExtensionSetsQueryEndToQueryLength) {
    // Verify that query_end equals query length after end extension
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 130);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    if (result.HasAlignment()) {
        const auto& aln = result.alignments[0];
        // With end extension, query_end should equal query length
        EXPECT_EQ(aln.query_end, static_cast<int32_t>(query.size()));
    }
}

TEST_F(ClassicalPipelineTest, EndExtensionWithoutRefSeqsSoftClips) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    // Do NOT set reference sequences - should fall back to soft clipping
    auto result = pipeline.AlignRead("read1", query);

    ASSERT_TRUE(result.HasAlignment());
    const auto& aln = result.alignments[0];

    // Should still cover full query (with soft clips at edges if anchors don't cover)
    EXPECT_EQ(aln.query_start, 0);
    EXPECT_EQ(aln.query_end, static_cast<int32_t>(query.size()));
}
