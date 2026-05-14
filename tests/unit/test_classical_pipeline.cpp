#include <gtest/gtest.h>

#include <random>
#include <set>

#include "classical/classical_pipeline.h"
#include "classical/classical_pipeline_internal.h"
#include "core/thread_pool.h"

using namespace llmap::classical;
using namespace llmap::core;

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

// === Parallel Alignment Tests (Phase B.1) ===

TEST_F(ClassicalPipelineTest, ParallelAlignBatchMultipleReads) {
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

    ThreadPool pool(4);
    auto results = pipeline.AlignReadsParallel(query_names, query_seqs, pool);

    EXPECT_EQ(results.size(), 3);
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_EQ(results[i].query_name, query_names[i]);
    }
}

TEST_F(ClassicalPipelineTest, ParallelAndSequentialProduceSameResults) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    for (int i = 0; i < 20; ++i) {
        query_names.push_back("read" + std::to_string(i));
        query_seqs.push_back(ref_seq.substr(50 * i % 800, 150));
    }

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto seq_results = pipeline.AlignReads(query_names, query_seqs);

    ThreadPool pool(4);
    auto par_results = pipeline.AlignReadsParallel(query_names, query_seqs, pool);

    ASSERT_EQ(seq_results.size(), par_results.size());
    for (size_t i = 0; i < seq_results.size(); ++i) {
        EXPECT_EQ(seq_results[i].query_name, par_results[i].query_name);
        EXPECT_EQ(seq_results[i].HasAlignment(), par_results[i].HasAlignment());
        if (seq_results[i].HasAlignment() && par_results[i].HasAlignment()) {
            EXPECT_EQ(seq_results[i].alignments[0].ref_start,
                      par_results[i].alignments[0].ref_start);
        }
    }
}

TEST_F(ClassicalPipelineTest, ParallelStatsAggregatedCorrectly) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    for (int i = 0; i < 10; ++i) {
        query_names.push_back("read" + std::to_string(i));
        query_seqs.push_back(ref_seq.substr(50 * i % 800, 150));
    }

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    ThreadPool pool(4);
    pipeline.AlignReadsParallel(query_names, query_seqs, pool);
    const auto& stats = pipeline.Stats();

    EXPECT_EQ(stats.total_reads, 10);
    EXPECT_EQ(stats.reads_aligned + stats.reads_unmapped, stats.total_reads);
    EXPECT_GT(stats.total_time_ms, 0.0f);
}

TEST_F(ClassicalPipelineTest, ParallelEmptyBatch) {
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", MakeRandomSequence(500))
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> empty_names;
    std::vector<std::string> empty_seqs;

    ThreadPool pool(2);
    auto results = pipeline.AlignReadsParallel(empty_names, empty_seqs, pool);

    EXPECT_TRUE(results.empty());
}

TEST_F(ClassicalPipelineTest, ParallelSingleRead) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> names = {"read1"};
    std::vector<std::string> seqs = {query};

    ThreadPool pool(4);
    auto results = pipeline.AlignReadsParallel(names, seqs, pool);

    ASSERT_EQ(results.size(), 1);
    EXPECT_TRUE(results[0].HasAlignment());
}

TEST_F(ClassicalPipelineTest, ParallelLargeBatch) {
    std::string ref_seq = MakeRandomSequence(2000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    for (int i = 0; i < 100; ++i) {
        query_names.push_back("read" + std::to_string(i));
        query_seqs.push_back(ref_seq.substr((50 * i) % 1800, 150));
    }

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    ThreadPool pool(4);
    auto results = pipeline.AlignReadsParallel(query_names, query_seqs, pool);

    EXPECT_EQ(results.size(), 100);

    size_t aligned_count = 0;
    for (const auto& r : results) {
        if (r.HasAlignment()) ++aligned_count;
    }
    EXPECT_GT(aligned_count, 50);  // Most should align
}

TEST_F(ClassicalPipelineTest, ParallelIdentityStatsAccurate) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    // All exact matches - should have very high identity
    for (int i = 0; i < 10; ++i) {
        query_names.push_back("read" + std::to_string(i));
        query_seqs.push_back(ref_seq.substr(100 + i * 50, 150));
    }

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    ThreadPool pool(4);
    pipeline.AlignReadsParallel(query_names, query_seqs, pool);
    const auto& stats = pipeline.Stats();

    // With exact matches, average identity should be reasonable
    if (stats.reads_aligned > 0) {
        EXPECT_GT(stats.avg_identity, 0.7f);
    }
}

// === Zero-Allocation Chaining Tests (Phase B.2) ===

TEST_F(ClassicalPipelineTest, ZeroAllocChainingProducesCorrectResults) {
    // Verify that zero-allocation chaining produces same results as before
    std::string ref_seq = MakeRandomSequence(1000, 1);
    std::string query = ref_seq.substr(100, 200);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    EXPECT_TRUE(result.HasAlignment());
    EXPECT_GT(result.num_chains, 0);
    EXPECT_GT(result.chaining_time_us, 0.0f);
}

TEST_F(ClassicalPipelineTest, ZeroAllocChainingConsistentAcrossMultipleReads) {
    // Verify scratch buffer reuse doesn't corrupt subsequent alignments
    std::string ref_seq = MakeRandomSequence(2000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    for (int i = 0; i < 50; ++i) {
        query_names.push_back("read" + std::to_string(i));
        query_seqs.push_back(ref_seq.substr((100 * i) % 1800, 150));
    }

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto results = pipeline.AlignReads(query_names, query_seqs);

    ASSERT_EQ(results.size(), 50);

    // Verify each result has correct query name
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_EQ(results[i].query_name, query_names[i]);
    }

    // Count aligned reads
    size_t aligned = 0;
    for (const auto& r : results) {
        if (r.HasAlignment()) ++aligned;
    }
    EXPECT_GT(aligned, 30);  // Most should align
}

TEST_F(ClassicalPipelineTest, ZeroAllocChainingParallelConsistency) {
    // Verify thread-local scratch buffers work correctly in parallel
    std::string ref_seq = MakeRandomSequence(2000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    for (int i = 0; i < 100; ++i) {
        query_names.push_back("read" + std::to_string(i));
        query_seqs.push_back(ref_seq.substr((50 * i) % 1800, 150));
    }

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    ThreadPool pool(8);  // Use 8 threads to stress test thread-local scratch
    auto results = pipeline.AlignReadsParallel(query_names, query_seqs, pool);

    ASSERT_EQ(results.size(), 100);

    // Verify results are consistent (correct query names, no corruption)
    std::set<std::string> seen_names;
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_EQ(results[i].query_name, query_names[i]);
        seen_names.insert(results[i].query_name);
    }
    EXPECT_EQ(seen_names.size(), 100);  // All unique
}

TEST_F(ClassicalPipelineTest, ZeroAllocChainingRepeatedSingleRead) {
    // Verify scratch buffer correctly resets between calls
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    // Align same read multiple times
    for (int i = 0; i < 10; ++i) {
        auto result = pipeline.AlignRead("read" + std::to_string(i), query);
        EXPECT_TRUE(result.HasAlignment());
        EXPECT_GT(result.num_chains, 0);
    }
}

TEST_F(ClassicalPipelineTest, ZeroAllocChainingVaryingSizes) {
    // Test with varying read sizes to stress scratch buffer resizing
    std::string ref_seq = MakeRandomSequence(3000, 1);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    // Vary read lengths significantly
    std::vector<size_t> lengths = {50, 100, 200, 500, 100, 50, 300, 150};

    for (size_t i = 0; i < lengths.size(); ++i) {
        size_t len = lengths[i];
        size_t offset = (i * 100) % (3000 - len);
        std::string query = ref_seq.substr(offset, len);

        auto result = pipeline.AlignRead("read" + std::to_string(i), query);

        if (len >= 50) {  // Long enough to have hits
            EXPECT_GE(result.num_hits, 0);  // May or may not have hits
        }
    }
}

// === Identity Filter Tests (Phase C.1) ===

TEST_F(ClassicalPipelineTest, DefaultMinIdentityIs080) {
    // Verify the new default of 0.80 for precision improvement
    ClassicalPipelineConfig default_config;
    EXPECT_FLOAT_EQ(default_config.min_identity, 0.80f);
}

TEST_F(ClassicalPipelineTest, IdentityFilterTracksFilteredAlignments) {
    // Create a read that will have low identity when aligned to wrong position
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = MakeRandomSequence(150, 99);  // Unrelated sequence

    config.min_identity = 0.50f;  // Low threshold so we might get hits
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    // The result should track filtered alignments
    // (filtered_by_identity + filtered_by_length + alignments) should be
    // related to chains extended
    EXPECT_GE(result.filtered_by_identity, 0);
    EXPECT_GE(result.filtered_by_length, 0);
}

TEST_F(ClassicalPipelineTest, IdentityFilterStrictThresholdFiltersMore) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);  // Exact match

    // First with permissive threshold
    config.min_identity = 0.50f;
    ClassicalPipeline pipeline_permissive(config);
    auto index1 = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline_permissive.SetIndex(std::move(index1));
    auto result_permissive = pipeline_permissive.AlignRead("read1", query);

    // Then with strict threshold
    config.min_identity = 0.99f;
    ClassicalPipeline pipeline_strict(config);
    auto index2 = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline_strict.SetIndex(std::move(index2));
    auto result_strict = pipeline_strict.AlignRead("read1", query);

    // Strict filter should filter more (or equal) alignments
    EXPECT_GE(result_strict.filtered_by_identity,
              result_permissive.filtered_by_identity);
}

TEST_F(ClassicalPipelineTest, IdentityFilterStatsAggregatedInBatch) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    // Mix of exact matches and random sequences
    for (int i = 0; i < 10; ++i) {
        query_names.push_back("read" + std::to_string(i));
        if (i % 2 == 0) {
            query_seqs.push_back(ref_seq.substr(100 + i * 50, 150));  // Exact
        } else {
            query_seqs.push_back(MakeRandomSequence(150, 100 + i));  // Random
        }
    }

    config.min_identity = 0.80f;  // Use new default
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    pipeline.AlignReads(query_names, query_seqs);
    const auto& stats = pipeline.Stats();

    // Stats should track filtering
    EXPECT_EQ(stats.total_reads, 10);
    EXPECT_GE(stats.alignments_filtered_by_identity, 0);
    EXPECT_GE(stats.alignments_filtered_by_length, 0);
}

TEST_F(ClassicalPipelineTest, IdentityFilterStatsAggregatedInParallel) {
    std::string ref_seq = MakeRandomSequence(1000, 1);

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    for (int i = 0; i < 20; ++i) {
        query_names.push_back("read" + std::to_string(i));
        if (i % 3 == 0) {
            query_seqs.push_back(ref_seq.substr(100 + (i * 30) % 700, 150));
        } else {
            query_seqs.push_back(MakeRandomSequence(150, 200 + i));
        }
    }

    config.min_identity = 0.80f;
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    ThreadPool pool(4);
    pipeline.AlignReadsParallel(query_names, query_seqs, pool);
    const auto& stats = pipeline.Stats();

    EXPECT_EQ(stats.total_reads, 20);
    EXPECT_GE(stats.alignments_filtered_by_identity, 0);
    EXPECT_GE(stats.alignments_filtered_by_length, 0);
}

TEST_F(ClassicalPipelineTest, HighIdentityExactMatchPassesFilter) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);  // Exact match

    config.min_identity = 0.80f;  // New default
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    // Exact match should pass the identity filter
    EXPECT_TRUE(result.HasAlignment());
    if (result.HasAlignment()) {
        EXPECT_GE(result.alignments[0].identity, 0.80f);
    }
}

TEST_F(ClassicalPipelineTest, PerReadFilterStatsCorrect) {
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 150);

    config.min_identity = 0.50f;
    config.min_aligned_bases = 20;
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    // Verify filter stats are non-negative and consistent
    EXPECT_GE(result.filtered_by_identity, 0);
    EXPECT_GE(result.filtered_by_length, 0);
    EXPECT_GE(result.chains_extended, 0);
}

// === MAPQ Calculation Tests (Phase C.3) ===

TEST(MapqCalculationTest, UniqueHighScoreHighIdentityGetsMapq60) {
    // High score, high identity, no secondaries = maximum confidence
    uint32_t mapq = ComputeMapq(
        200,    // primary_score: good score
        0,      // secondary_score: no secondary
        0.95f,  // identity: high
        0,      // num_secondaries: none
        150);   // query_len

    EXPECT_EQ(mapq, 60);
}

TEST(MapqCalculationTest, UniqueModerateScoreGetsHighMapq) {
    // Moderate score, decent identity, no secondaries
    uint32_t mapq = ComputeMapq(
        100,    // primary_score
        0,      // secondary_score
        0.85f,  // identity
        0,      // num_secondaries
        150);   // query_len

    EXPECT_GE(mapq, 40);
    EXPECT_LE(mapq, 60);
}

TEST(MapqCalculationTest, MultiMappingWithScoreGapGetsReducedMapq) {
    // Primary better than secondary by a margin
    uint32_t mapq = ComputeMapq(
        200,    // primary_score
        150,    // secondary_score: 50 point gap
        0.90f,  // identity
        1,      // num_secondaries
        150);   // query_len

    // Should have reduced MAPQ due to secondary
    EXPECT_GT(mapq, 0);
    EXPECT_LT(mapq, 60);
}

TEST(MapqCalculationTest, MultiMappingWithNoScoreGapGetsMapq0) {
    // Primary and secondary have same score = ambiguous
    uint32_t mapq = ComputeMapq(
        200,    // primary_score
        200,    // secondary_score: same score
        0.90f,  // identity
        1,      // num_secondaries
        150);   // query_len

    EXPECT_EQ(mapq, 0);
}

TEST(MapqCalculationTest, LowIdentityGetsMapq0) {
    // Even unique alignment with low identity = uncertain
    uint32_t mapq = ComputeMapq(
        200,    // primary_score
        0,      // secondary_score
        0.40f,  // identity: very low
        0,      // num_secondaries
        150);   // query_len

    EXPECT_EQ(mapq, 0);
}

TEST(MapqCalculationTest, ZeroScoreGetsMapq0) {
    uint32_t mapq = ComputeMapq(
        0,      // primary_score: zero
        0,      // secondary_score
        0.90f,  // identity
        0,      // num_secondaries
        150);   // query_len

    EXPECT_EQ(mapq, 0);
}

TEST(MapqCalculationTest, NegativeScoreGetsMapq0) {
    uint32_t mapq = ComputeMapq(
        -50,    // primary_score: negative
        0,      // secondary_score
        0.90f,  // identity
        0,      // num_secondaries
        150);   // query_len

    EXPECT_EQ(mapq, 0);
}

TEST(MapqCalculationTest, ManySecondariesReduceMapq) {
    // Many secondaries should reduce confidence more
    uint32_t mapq_one = ComputeMapq(
        200, 100, 0.90f, 1, 150);

    uint32_t mapq_many = ComputeMapq(
        200, 100, 0.90f, 10, 150);

    // More secondaries = lower MAPQ
    EXPECT_LT(mapq_many, mapq_one);
}

TEST(MapqCalculationTest, LargeScoreGapGivesHighMapq) {
    // Large gap between primary and secondary = confident
    uint32_t mapq = ComputeMapq(
        300,    // primary_score
        50,     // secondary_score: large gap (250)
        0.95f,  // identity
        1,      // num_secondaries
        150);   // query_len

    // With score_diff=250, p_secondary=exp(-250/30)≈0.00024, mapq≈36
    // Reasonably high for multi-mapping case - better than ambiguous
    EXPECT_GE(mapq, 30);  // Should be meaningfully above zero
    EXPECT_LE(mapq, 60);
}

TEST(MapqCalculationTest, MapqInValidRange) {
    // Test various combinations stay in [0, 60]
    for (int32_t score = 0; score <= 500; score += 50) {
        for (int32_t sec_score = 0; sec_score <= score; sec_score += 25) {
            for (float identity = 0.5f; identity <= 1.0f; identity += 0.1f) {
                for (uint32_t num_sec = 0; num_sec <= 5; ++num_sec) {
                    uint32_t mapq = ComputeMapq(
                        score, sec_score, identity, num_sec, 150);
                    EXPECT_GE(mapq, 0);
                    EXPECT_LE(mapq, 60);
                }
            }
        }
    }
}

TEST_F(ClassicalPipelineTest, MapqUniqueAlignmentHighConfidence) {
    // Test that unique alignment to exact match gets high MAPQ
    std::string ref_seq = MakeRandomSequence(500, 1);
    std::string query = ref_seq.substr(100, 180);  // Exact match

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> refs = {ref_seq};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    ASSERT_TRUE(result.HasAlignment());
    const auto& aln = result.alignments[0];

    // Unique exact match should have high MAPQ
    EXPECT_GE(aln.mapq, 40);
}

TEST_F(ClassicalPipelineTest, MapqSecondaryAlwaysZero) {
    // Secondary alignments should always have MAPQ=0
    std::string ref_seq = MakeRandomSequence(1000, 1);
    // Duplicate region to create secondaries
    std::string repeated = ref_seq;
    repeated.replace(500, 200, ref_seq.substr(100, 200));
    std::string query = ref_seq.substr(100, 150);

    config.report_secondary = true;
    config.max_alignments = 5;
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", repeated)
        .Build();
    pipeline.SetIndex(std::move(index));

    auto result = pipeline.AlignRead("read1", query);

    // Check all secondary alignments have MAPQ=0
    for (const auto& aln : result.alignments) {
        if (!aln.is_primary) {
            EXPECT_EQ(aln.mapq, 0);
        }
    }
}

TEST_F(ClassicalPipelineTest, MapqMultiMappingReducesConfidence) {
    // When there are multiple good alignments, MAPQ should be lower
    std::string ref_seq = MakeRandomSequence(1000, 1);
    // Create exact duplicate to force multi-mapping
    std::string ref_with_dup = ref_seq + ref_seq.substr(100, 200);
    std::string query = ref_seq.substr(100, 150);

    config.report_secondary = true;
    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_with_dup)
        .Build();
    pipeline.SetIndex(std::move(index));

    std::vector<std::string> refs = {ref_with_dup};
    pipeline.SetReferenceSequences(refs);

    auto result = pipeline.AlignRead("read1", query);

    if (result.alignments.size() > 1) {
        // Primary with secondaries should have reduced MAPQ
        EXPECT_LT(result.alignments[0].mapq, 60);
    }
}

TEST_F(ClassicalPipelineTest, MapqIsInValidRange) {
    // Verify MAPQ is always in [0, 60] for various inputs
    std::string ref_seq = MakeRandomSequence(500, 1);

    ClassicalPipeline pipeline(config);
    auto index = MinimizerIndex::Builder(config.minimizer_config)
        .AddSequence("ref1", ref_seq)
        .Build();
    pipeline.SetIndex(std::move(index));

    // Test with various query types
    std::vector<std::string> queries = {
        ref_seq.substr(100, 150),           // Exact match
        MakeRandomSequence(150, 99),        // Unrelated
        ref_seq.substr(0, 50),              // Short
    };

    for (const auto& query : queries) {
        auto result = pipeline.AlignRead("read", query);
        for (const auto& aln : result.alignments) {
            EXPECT_GE(aln.mapq, 0);
            EXPECT_LE(aln.mapq, 60);
        }
    }
}
