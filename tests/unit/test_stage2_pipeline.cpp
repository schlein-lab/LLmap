// LLmap — Unit tests for Stage 2 pipeline orchestrator.

#include <gtest/gtest.h>

#include "reference_collapse/stage2_pipeline.h"
#include "reference_collapse/reference_index.h"
#include "self_interference/cluster_rep.h"
#include "core/wave_state.h"

#include <memory>
#include <vector>

namespace llmap {
namespace {

// ========== Configuration Tests ==========

TEST(Stage2ConfigTest, DefaultValues) {
    Stage2Config config;
    EXPECT_EQ(config.max_iterations_per_level, 50u);
    EXPECT_FLOAT_EQ(config.gamma, 0.3f);
    EXPECT_FLOAT_EQ(config.tau_collapse, 0.99f);
    EXPECT_FLOAT_EQ(config.convergence_delta, 1e-5f);
    EXPECT_TRUE(config.apply_smoothing);
    EXPECT_FLOAT_EQ(config.expansion_threshold, 0.01f);
    EXPECT_EQ(config.max_candidates, 50u);
    EXPECT_FLOAT_EQ(config.refine_trigger_rate, 0.8f);
    EXPECT_FLOAT_EQ(config.propagation_confidence, 0.8f);
    EXPECT_FLOAT_EQ(config.similarity_weight, 0.5f);
    EXPECT_EQ(config.num_threads, 0u);
    EXPECT_FALSE(config.verbose);
}

TEST(Stage2ConfigTest, ModifyAllFields) {
    Stage2Config config;
    config.index_path = "/some/path";
    config.max_iterations_per_level = 100;
    config.gamma = 0.5f;
    config.tau_collapse = 0.95f;
    config.convergence_delta = 1e-6f;
    config.apply_smoothing = false;
    config.expansion_threshold = 0.02f;
    config.max_candidates = 100;
    config.refine_trigger_rate = 0.9f;
    config.propagation_confidence = 0.9f;
    config.similarity_weight = 0.7f;
    config.num_threads = 8;
    config.verbose = true;

    EXPECT_EQ(config.index_path, "/some/path");
    EXPECT_EQ(config.max_iterations_per_level, 100u);
    EXPECT_FLOAT_EQ(config.gamma, 0.5f);
    EXPECT_FLOAT_EQ(config.tau_collapse, 0.95f);
    EXPECT_FLOAT_EQ(config.convergence_delta, 1e-6f);
    EXPECT_FALSE(config.apply_smoothing);
    EXPECT_FLOAT_EQ(config.expansion_threshold, 0.02f);
    EXPECT_EQ(config.max_candidates, 100u);
    EXPECT_FLOAT_EQ(config.refine_trigger_rate, 0.9f);
    EXPECT_FLOAT_EQ(config.propagation_confidence, 0.9f);
    EXPECT_FLOAT_EQ(config.similarity_weight, 0.7f);
    EXPECT_EQ(config.num_threads, 8u);
    EXPECT_TRUE(config.verbose);
}

// ========== Stats Structure Tests ==========

TEST(LevelStatsTest, DefaultValues) {
    LevelStats stats;
    EXPECT_EQ(stats.iterations, 0u);
    EXPECT_EQ(stats.reads_processed, 0u);
    EXPECT_EQ(stats.reads_collapsed, 0u);
    EXPECT_EQ(stats.reads_refined, 0u);
    EXPECT_FLOAT_EQ(stats.avg_entropy, 0.0f);
    EXPECT_FLOAT_EQ(stats.time_ms, 0.0f);
}

TEST(Stage2StatsTest, DefaultValues) {
    Stage2Stats stats;
    EXPECT_EQ(stats.num_representatives, 0u);
    EXPECT_EQ(stats.num_clusters, 0u);
    EXPECT_EQ(stats.total_reads, 0u);
    EXPECT_EQ(stats.ref_targets, 0u);
    EXPECT_EQ(stats.ref_l2_buckets, 0u);
    EXPECT_EQ(stats.members_propagated, 0u);
    EXPECT_FLOAT_EQ(stats.propagation_time_ms, 0.0f);
    EXPECT_EQ(stats.total_iterations, 0u);
    EXPECT_FLOAT_EQ(stats.total_time_ms, 0.0f);
    EXPECT_FLOAT_EQ(stats.collapse_rate, 0.0f);
}

// ========== ReadAlignment Tests ==========

TEST(ReadAlignmentTest, FieldsCanBeSet) {
    ReadAlignment alignment{};
    alignment.read_idx = 42;
    alignment.target_name = "chr1";
    alignment.position = 1000000;
    alignment.confidence = 0.95f;
    alignment.is_representative = true;
    alignment.cluster_id = 5;
    alignment.l2_bucket_id = 100;

    EXPECT_EQ(alignment.read_idx, 42u);
    EXPECT_EQ(alignment.target_name, "chr1");
    EXPECT_EQ(alignment.position, 1000000u);
    EXPECT_FLOAT_EQ(alignment.confidence, 0.95f);
    EXPECT_TRUE(alignment.is_representative);
    EXPECT_EQ(alignment.cluster_id, 5u);
    EXPECT_EQ(alignment.l2_bucket_id, 100u);
}

// ========== Stage2Result Tests ==========

TEST(Stage2ResultTest, DefaultValues) {
    Stage2Result result;
    EXPECT_TRUE(result.alignments.empty());
    EXPECT_FALSE(result.success);
    EXPECT_TRUE(result.error_message.empty());
    EXPECT_EQ(result.NumAlignments(), 0u);
}

TEST(Stage2ResultTest, GetAlignmentsForRead) {
    Stage2Result result;

    ReadAlignment a1{};
    a1.read_idx = 0;
    a1.target_name = "chr1";
    a1.position = 100;

    ReadAlignment a2{};
    a2.read_idx = 1;
    a2.target_name = "chr2";
    a2.position = 200;

    ReadAlignment a3{};
    a3.read_idx = 0;
    a3.target_name = "chr1";
    a3.position = 150;

    result.alignments.push_back(a1);
    result.alignments.push_back(a2);
    result.alignments.push_back(a3);

    auto read0_alignments = result.GetAlignmentsForRead(0);
    EXPECT_EQ(read0_alignments.size(), 2u);

    auto read1_alignments = result.GetAlignmentsForRead(1);
    EXPECT_EQ(read1_alignments.size(), 1u);

    auto read99_alignments = result.GetAlignmentsForRead(99);
    EXPECT_TRUE(read99_alignments.empty());
}

TEST(Stage2ResultTest, GetHighConfidenceAlignments) {
    Stage2Result result;

    ReadAlignment a1{}, a2{}, a3{}, a4{};
    a1.read_idx = 0; a1.confidence = 0.5f;
    a2.read_idx = 1; a2.confidence = 0.8f;
    a3.read_idx = 2; a3.confidence = 0.95f;
    a4.read_idx = 3; a4.confidence = 0.99f;

    result.alignments.push_back(a1);
    result.alignments.push_back(a2);
    result.alignments.push_back(a3);
    result.alignments.push_back(a4);

    auto high_conf = result.GetHighConfidenceAlignments(0.9f);
    EXPECT_EQ(high_conf.size(), 2u);

    auto all_above_zero = result.GetHighConfidenceAlignments(0.0f);
    EXPECT_EQ(all_above_zero.size(), 4u);

    auto none = result.GetHighConfidenceAlignments(1.0f);
    EXPECT_TRUE(none.empty());
}

// ========== Pipeline Construction Tests ==========

TEST(Stage2PipelineTest, DefaultConstruction) {
    Stage2Pipeline pipeline;
    auto config = pipeline.Config();
    EXPECT_TRUE(config.index_path.empty());
}

TEST(Stage2PipelineTest, ConstructWithConfig) {
    Stage2Config config;
    config.gamma = 0.5f;
    config.max_iterations_per_level = 100;

    Stage2Pipeline pipeline(config);
    EXPECT_FLOAT_EQ(pipeline.Config().gamma, 0.5f);
    EXPECT_EQ(pipeline.Config().max_iterations_per_level, 100u);
}

TEST(Stage2PipelineTest, SetConfig) {
    Stage2Pipeline pipeline;

    Stage2Config config;
    config.tau_collapse = 0.95f;
    pipeline.SetConfig(config);

    EXPECT_FLOAT_EQ(pipeline.Config().tau_collapse, 0.95f);
}

TEST(Stage2PipelineTest, MoveConstruction) {
    Stage2Config config;
    config.gamma = 0.6f;

    Stage2Pipeline pipeline1(config);
    Stage2Pipeline pipeline2(std::move(pipeline1));

    EXPECT_FLOAT_EQ(pipeline2.Config().gamma, 0.6f);
}

TEST(Stage2PipelineTest, MoveAssignment) {
    Stage2Config config;
    config.gamma = 0.7f;

    Stage2Pipeline pipeline1(config);
    Stage2Pipeline pipeline2;
    pipeline2 = std::move(pipeline1);

    EXPECT_FLOAT_EQ(pipeline2.Config().gamma, 0.7f);
}

// ========== LoadReferenceIndex Tests ==========

TEST(Stage2PipelineTest, LoadReferenceIndex_NoPath) {
    Stage2Pipeline pipeline;
    EXPECT_FALSE(pipeline.LoadReferenceIndex());
    EXPECT_FALSE(pipeline.LastError().empty());
}

TEST(Stage2PipelineTest, LoadReferenceIndex_InvalidPath) {
    Stage2Config config;
    config.index_path = "/nonexistent/path/to/index.llmap";

    Stage2Pipeline pipeline(config);
    EXPECT_FALSE(pipeline.LoadReferenceIndex());
    EXPECT_FALSE(pipeline.LastError().empty());
}

// ========== SetReferenceIndex Tests ==========

TEST(Stage2PipelineTest, SetReferenceIndex) {
    ReferenceIndexConfig idx_config;
    idx_config.l1_granularity = 1'000'000;
    idx_config.l2_granularity = 100'000;

    ReferenceIndex::Builder builder(idx_config);
    builder.AddTarget({"chr1", 10'000'000, "abc123"});

    auto index = builder.Build();
    ASSERT_NE(index, nullptr);

    Stage2Pipeline pipeline;
    pipeline.SetReferenceIndex(std::move(index));

    EXPECT_TRUE(pipeline.LoadReferenceIndex());
}

// ========== InitializeWaveState Tests ==========

TEST(Stage2PipelineTest, InitializeWaveState_NoIndex) {
    Stage2Pipeline pipeline;

    self_interference::ClusterRepResult rep_result;
    rep_result.representatives.push_back({0, 0, 0.9f, 1, 0.0f, 0.0f});
    rep_result.representatives.push_back({1, 1, 0.8f, 1, 0.0f, 0.0f});
    rep_result.representatives.push_back({2, 2, 0.7f, 1, 0.0f, 0.0f});

    EXPECT_FALSE(pipeline.InitializeWaveState(rep_result));
}

TEST(Stage2PipelineTest, InitializeWaveState_WithIndex) {
    ReferenceIndexConfig idx_config;
    idx_config.l1_granularity = 1'000'000;
    idx_config.l2_granularity = 100'000;

    ReferenceIndex::Builder builder(idx_config);
    builder.AddTarget({"chr1", 5'000'000, "abc123"});

    auto index = builder.Build();
    ASSERT_NE(index, nullptr);

    Stage2Pipeline pipeline;
    pipeline.SetReferenceIndex(std::move(index));
    ASSERT_TRUE(pipeline.LoadReferenceIndex());

    self_interference::ClusterRepResult rep_result;
    rep_result.representatives.push_back({0, 0, 0.9f, 1, 0.0f, 0.0f});
    rep_result.representatives.push_back({1, 1, 0.8f, 1, 0.0f, 0.0f});
    rep_result.representatives.push_back({2, 2, 0.7f, 1, 0.0f, 0.0f});

    EXPECT_TRUE(pipeline.InitializeWaveState(rep_result));

    auto* ws = pipeline.GetWaveState();
    ASSERT_NE(ws, nullptr);
    EXPECT_EQ(ws->n_reads(), 3u);
}

TEST(Stage2PipelineTest, InitializeWaveState_EmptyReps) {
    ReferenceIndexConfig idx_config;
    idx_config.l1_granularity = 1'000'000;
    idx_config.l2_granularity = 100'000;

    ReferenceIndex::Builder builder(idx_config);
    builder.AddTarget({"chr1", 5'000'000, "abc"});

    auto index = builder.Build();
    ASSERT_NE(index, nullptr);

    Stage2Pipeline pipeline;
    pipeline.SetReferenceIndex(std::move(index));
    ASSERT_TRUE(pipeline.LoadReferenceIndex());

    self_interference::ClusterRepResult rep_result;
    // Empty representatives

    EXPECT_TRUE(pipeline.InitializeWaveState(rep_result));

    auto* ws = pipeline.GetWaveState();
    ASSERT_NE(ws, nullptr);
    EXPECT_EQ(ws->n_reads(), 0u);
}

// ========== RunEMLevel Tests ==========

TEST(Stage2PipelineTest, RunEMLevel_NoWaveState) {
    Stage2Pipeline pipeline;
    EXPECT_FALSE(pipeline.RunEMLevel(WaveLevel::L0));
}

TEST(Stage2PipelineTest, RunEMLevel_WithWaveState) {
    ReferenceIndexConfig idx_config;
    idx_config.l1_granularity = 1'000'000;
    idx_config.l2_granularity = 100'000;

    ReferenceIndex::Builder builder(idx_config);
    builder.AddTarget({"chr1", 2'000'000, "abc"});

    auto index = builder.Build();
    ASSERT_NE(index, nullptr);

    Stage2Config config;
    config.max_iterations_per_level = 2;  // Quick test

    Stage2Pipeline pipeline(config);
    pipeline.SetReferenceIndex(std::move(index));
    ASSERT_TRUE(pipeline.LoadReferenceIndex());

    self_interference::ClusterRepResult rep_result;
    rep_result.representatives.push_back({0, 0, 0.9f, 1, 0.0f, 0.0f});

    ASSERT_TRUE(pipeline.InitializeWaveState(rep_result));
    EXPECT_TRUE(pipeline.RunEMLevel(WaveLevel::L0));
}

// ========== RefineLevel Tests ==========

TEST(Stage2PipelineTest, RefineLevel_NoWaveState) {
    Stage2Pipeline pipeline;
    EXPECT_FALSE(pipeline.RefineLevel(WaveLevel::L0));
}

TEST(Stage2PipelineTest, RefineLevel_L2NotAllowed) {
    ReferenceIndexConfig idx_config;
    idx_config.l1_granularity = 1'000'000;
    idx_config.l2_granularity = 100'000;

    ReferenceIndex::Builder builder(idx_config);
    builder.AddTarget({"chr1", 2'000'000, "abc"});

    auto index = builder.Build();
    ASSERT_NE(index, nullptr);

    Stage2Pipeline pipeline;
    pipeline.SetReferenceIndex(std::move(index));
    ASSERT_TRUE(pipeline.LoadReferenceIndex());

    self_interference::ClusterRepResult rep_result;
    rep_result.representatives.push_back({0, 0, 0.9f, 1, 0.0f, 0.0f});

    ASSERT_TRUE(pipeline.InitializeWaveState(rep_result));

    EXPECT_FALSE(pipeline.RefineLevel(WaveLevel::L2));
}

// ========== ExtractAlignments Tests ==========

TEST(Stage2PipelineTest, ExtractAlignments_NoState) {
    Stage2Pipeline pipeline;
    EXPECT_FALSE(pipeline.ExtractAlignments());
}

}  // namespace
}  // namespace llmap
