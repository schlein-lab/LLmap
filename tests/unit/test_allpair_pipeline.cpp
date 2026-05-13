// Unit tests for AllpairPipeline (Stage 1 Self-Interference)

#include <gtest/gtest.h>

#include "self_interference/allpair_pipeline.h"
#include "self_interference/faiss_wrapper.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace llmap::self_interference {
namespace {

// Macro to skip tests when FAISS is not available
#define SKIP_IF_NO_FAISS() \
    if (!IsFaissAvailable()) { \
        GTEST_SKIP() << "FAISS not available - skipping test"; \
    }

class AllpairPipelineTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_allpair_test";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path CreateTestFastq(const std::string& name, size_t num_reads,
                                           size_t read_length = 100) {
        auto path = test_dir_ / name;
        std::ofstream out(path);

        // Create reads with some similarity patterns
        std::string bases = "ACGT";
        for (size_t i = 0; i < num_reads; ++i) {
            out << "@read" << i << "\n";

            // Generate sequence with some pattern
            std::string seq;
            size_t cluster = i % 5;  // 5 implicit clusters
            for (size_t j = 0; j < read_length; ++j) {
                // Reads in same cluster share a prefix pattern
                if (j < 20) {
                    seq += bases[(cluster + j) % 4];
                } else {
                    seq += bases[(i + j) % 4];
                }
            }
            out << seq << "\n";
            out << "+\n";
            out << std::string(read_length, 'I') << "\n";
        }

        return path;
    }

    std::vector<std::string> ReadOutputLines(const std::filesystem::path& path) {
        std::vector<std::string> lines;
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line)) {
            lines.push_back(line);
        }
        return lines;
    }

    std::filesystem::path test_dir_;
};

// ========== Configuration Tests ==========

TEST_F(AllpairPipelineTest, DefaultConfig) {
    AllpairConfig config;

    EXPECT_EQ(config.embedding_batch_size, 1024);
    EXPECT_EQ(config.embedding_dim, 256);
    EXPECT_EQ(config.faiss_k, 50);
    EXPECT_EQ(config.faiss_nlist, 1024);
    EXPECT_FLOAT_EQ(config.leiden_resolution, 1.0f);
    EXPECT_EQ(config.min_cluster_size, 2);
    EXPECT_EQ(config.swc_max_iterations, 20);
    EXPECT_FLOAT_EQ(config.swc_collapse_threshold, 0.95f);
    EXPECT_FLOAT_EQ(config.swc_gamma, 0.3f);
    EXPECT_TRUE(config.select_representatives);
    EXPECT_FALSE(config.verbose);
    EXPECT_FALSE(config.use_gpu_embedder);
    EXPECT_FALSE(config.use_gpu_faiss);
}

TEST_F(AllpairPipelineTest, ConfigModification) {
    AllpairConfig config;
    config.faiss_k = 100;
    config.leiden_resolution = 0.5f;
    config.min_cluster_size = 5;
    config.verbose = true;

    AllpairPipeline pipeline(config);
    const auto& actual = pipeline.GetConfig();

    EXPECT_EQ(actual.faiss_k, 100);
    EXPECT_FLOAT_EQ(actual.leiden_resolution, 0.5f);
    EXPECT_EQ(actual.min_cluster_size, 5);
    EXPECT_TRUE(actual.verbose);
}

// ========== Pipeline Execution ==========

TEST_F(AllpairPipelineTest, RunWithSmallDataset) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("small.fastq", 20, 50);
    auto output = test_dir_ / "small_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 64;  // Smaller for faster test
    config.faiss_k = 5;
    config.min_cluster_size = 1;

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(result->success) << "Error: " << result->error_message;

    // Check basic result properties
    EXPECT_EQ(result->NumReads(), 20);
    EXPECT_GT(result->NumClusters(), 0);
    EXPECT_LE(result->NumClusters(), 20);

    // All reads should have assignments
    EXPECT_EQ(result->cluster_ids.size(), 20);
    EXPECT_EQ(result->confidences.size(), 20);
    EXPECT_EQ(result->read_ids.size(), 20);

    // Check output file was created
    EXPECT_TRUE(std::filesystem::exists(output));

    auto lines = ReadOutputLines(output);
    EXPECT_EQ(lines.size(), 21);  // header + 20 records
}

TEST_F(AllpairPipelineTest, RunWithProgressCallback) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("progress.fastq", 10, 30);
    auto output = test_dir_ / "progress_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 32;
    config.faiss_k = 3;
    config.min_cluster_size = 1;

    std::vector<std::string> stages_seen;
    auto progress = [&stages_seen](const std::string& stage, size_t current,
                                     size_t total, const std::string& /*msg*/) {
        stages_seen.push_back(stage);
    };

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run(progress);

    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(result->success);

    // Should have seen all stages
    EXPECT_GE(stages_seen.size(), 6);

    // Check expected stages
    bool saw_load = false, saw_embed = false, saw_cluster = false;
    for (const auto& s : stages_seen) {
        if (s == "load") saw_load = true;
        if (s == "embed") saw_embed = true;
        if (s == "cluster") saw_cluster = true;
    }
    EXPECT_TRUE(saw_load);
    EXPECT_TRUE(saw_embed);
    EXPECT_TRUE(saw_cluster);
}

TEST_F(AllpairPipelineTest, StatisticsTracking) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("stats.fastq", 50, 100);
    auto output = test_dir_ / "stats_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 64;
    config.faiss_k = 10;
    config.min_cluster_size = 1;

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(result->success);

    // Check statistics
    EXPECT_EQ(result->stats.input_reads, 50);
    EXPECT_GT(result->stats.graph_nodes, 0);
    EXPECT_GT(result->stats.graph_edges, 0);
    EXPECT_GT(result->stats.num_clusters, 0);
    EXPECT_GT(result->stats.total_time_ms, 0);
    EXPECT_GT(result->stats.reads_per_second, 0);
}

TEST_F(AllpairPipelineTest, RepresentativeSelection) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("reps.fastq", 30, 80);
    auto output = test_dir_ / "reps_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 64;
    config.faiss_k = 8;
    config.min_cluster_size = 1;
    config.select_representatives = true;

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(result->success);

    // Check representatives were selected
    auto reps = result->GetRepresentatives();
    EXPECT_GT(reps.size(), 0);
    EXPECT_LE(reps.size(), result->NumClusters());

    // Each representative should be marked
    for (uint32_t rep : reps) {
        EXPECT_TRUE(result->is_representative[rep]);
    }

    // Check stat
    EXPECT_EQ(result->stats.num_representatives, reps.size());
}

// ========== Error Handling ==========

TEST_F(AllpairPipelineTest, MissingInputFile) {
    AllpairConfig config;
    config.reads_a = test_dir_ / "nonexistent.fastq";
    config.output = test_dir_ / "output.tsv";

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    EXPECT_FALSE(result->success);
    EXPECT_FALSE(result->error_message.empty());
}

TEST_F(AllpairPipelineTest, EmptyInputFile) {
    auto fastq = test_dir_ / "empty.fastq";
    { std::ofstream out(fastq); }  // Create empty file

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = test_dir_ / "output.tsv";

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    EXPECT_FALSE(result->success);
}

TEST_F(AllpairPipelineTest, SingleRead) {
    auto fastq = CreateTestFastq("single.fastq", 1, 50);
    auto output = test_dir_ / "single_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 32;
    config.faiss_k = 1;
    config.min_cluster_size = 1;

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    // Single read should fail as we need at least 2 for similarity graph
    EXPECT_FALSE(result->success);
}

TEST_F(AllpairPipelineTest, TwoReads) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("two.fastq", 2, 50);
    auto output = test_dir_ / "two_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 32;
    config.faiss_k = 1;
    config.min_cluster_size = 1;

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(result->success) << result->error_message;
    EXPECT_EQ(result->NumReads(), 2);
}

// ========== Configuration Effects ==========

TEST_F(AllpairPipelineTest, MaxReadsLimit) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("many.fastq", 100, 50);
    auto output = test_dir_ / "limited_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 32;
    config.faiss_k = 5;
    config.min_cluster_size = 1;
    config.max_reads = 25;  // Limit to 25 reads

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(result->success);

    EXPECT_EQ(result->NumReads(), 25);
    EXPECT_EQ(result->stats.input_reads, 25);
}

TEST_F(AllpairPipelineTest, HighResolutionMoreClusters) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("res.fastq", 40, 60);

    // Low resolution run
    AllpairConfig low_config;
    low_config.reads_a = fastq;
    low_config.output = test_dir_ / "low_res.tsv";
    low_config.embedding_dim = 64;
    low_config.faiss_k = 8;
    low_config.min_cluster_size = 1;
    low_config.leiden_resolution = 0.5f;

    AllpairPipeline low_pipeline(low_config);
    auto low_result = low_pipeline.Run();
    ASSERT_TRUE(low_result && low_result->success);

    // High resolution run
    AllpairConfig high_config;
    high_config.reads_a = fastq;
    high_config.output = test_dir_ / "high_res.tsv";
    high_config.embedding_dim = 64;
    high_config.faiss_k = 8;
    high_config.min_cluster_size = 1;
    high_config.leiden_resolution = 2.0f;

    AllpairPipeline high_pipeline(high_config);
    auto high_result = high_pipeline.Run();
    ASSERT_TRUE(high_result && high_result->success);

    // Higher resolution typically produces more clusters
    EXPECT_GE(high_result->NumClusters(), low_result->NumClusters());
}

// ========== Individual Stage Tests ==========

TEST_F(AllpairPipelineTest, LoadReadsStage) {
    auto fastq = CreateTestFastq("load.fastq", 15, 40);

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = test_dir_ / "dummy.tsv";
    config.embedding_dim = 32;

    AllpairPipeline pipeline(config);

    EXPECT_TRUE(pipeline.LoadReads());
}

TEST_F(AllpairPipelineTest, ComputeEmbeddingsStage) {
    auto fastq = CreateTestFastq("embed.fastq", 10, 30);

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = test_dir_ / "dummy.tsv";
    config.embedding_dim = 32;

    AllpairPipeline pipeline(config);

    ASSERT_TRUE(pipeline.LoadReads());
    EXPECT_TRUE(pipeline.ComputeEmbeddings());
}

// ========== Convenience Function ==========

TEST_F(AllpairPipelineTest, RunAllpairPipelineConvenience) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("conv.fastq", 15, 50);
    auto output = test_dir_ / "conv_output.tsv";

    auto result = RunAllpairPipeline(fastq, output);

    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(result->success);
    EXPECT_EQ(result->NumReads(), 15);
}

// ========== Output Verification ==========

TEST_F(AllpairPipelineTest, OutputContainsAllReads) {
    SKIP_IF_NO_FAISS();

    auto fastq = CreateTestFastq("verify.fastq", 25, 60);
    auto output = test_dir_ / "verify_output.tsv";

    AllpairConfig config;
    config.reads_a = fastq;
    config.output = output;
    config.embedding_dim = 48;
    config.faiss_k = 6;
    config.min_cluster_size = 1;

    AllpairPipeline pipeline(config);
    auto result = pipeline.Run();

    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(result->success);

    // Read output and verify
    auto lines = ReadOutputLines(output);
    EXPECT_EQ(lines.size(), 26);  // header + 25 records

    // Check all read IDs are present
    std::set<std::string> output_ids;
    for (size_t i = 1; i < lines.size(); ++i) {
        auto tab_pos = lines[i].find('\t');
        if (tab_pos != std::string::npos) {
            output_ids.insert(lines[i].substr(0, tab_pos));
        }
    }

    for (size_t i = 0; i < 25; ++i) {
        EXPECT_TRUE(output_ids.count("read" + std::to_string(i)) > 0)
            << "Missing read" << i;
    }
}

}  // namespace
}  // namespace llmap::self_interference
