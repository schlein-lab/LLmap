// LLmap — Unit tests for benchmark truth generation (Phase 11.2).

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "synthetic/benchmark_truth.h"

namespace llmap::synthetic {
namespace {

class BenchmarkTruthTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() /
                    ("llmap_bench_test_" + std::to_string(getpid()) + "_" +
                     std::to_string(test_counter_++));
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path test_dir_;
    static int test_counter_;
};

int BenchmarkTruthTest::test_counter_ = 0;

// ============================================================================
// T1 Generation Tests
// ============================================================================

TEST_F(BenchmarkTruthTest, T1ConfigDefaults) {
    T1Config config;
    EXPECT_EQ(config.seed, 42);
    EXPECT_EQ(config.n_reads, 1000000);
    EXPECT_EQ(config.coverage, 30);
    EXPECT_EQ(config.read_length, 10000);
    EXPECT_EQ(config.chromosomes.size(), 2);
    EXPECT_EQ(config.chromosomes[0], "chr20");
    EXPECT_EQ(config.chromosomes[1], "chr14");
}

TEST_F(BenchmarkTruthTest, T1GenerateSmall) {
    T1Config config;
    config.seed = 123;
    config.n_reads = 100;
    config.coverage = 10;
    config.region_length = 10000;

    auto dataset = BenchmarkGenerator::generate_t1(config);

    EXPECT_EQ(dataset.task, BenchmarkTask::T1_WGS);
    EXPECT_GT(dataset.reads.size(), 0);
    EXPECT_LE(dataset.reads.size(), config.n_reads);
    EXPECT_EQ(dataset.references.size(), 2);
    EXPECT_EQ(dataset.total_reads, dataset.reads.size());
    EXPECT_GT(dataset.total_bases, 0);
}

TEST_F(BenchmarkTruthTest, T1ReadsHaveGroundTruth) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 50;
    config.coverage = 5;
    config.region_length = 5000;

    auto dataset = BenchmarkGenerator::generate_t1(config);

    for (const auto& read : dataset.reads) {
        EXPECT_FALSE(read.id.empty());
        EXPECT_FALSE(read.sequence.empty());
        EXPECT_EQ(read.sequence.size(), read.quality.size());
        EXPECT_FALSE(read.chrom.empty());
        EXPECT_GE(read.pos, 0);
    }
}

TEST_F(BenchmarkTruthTest, T1ReferenceSequencesGenerated) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 10;
    config.region_length = 1000;

    auto dataset = BenchmarkGenerator::generate_t1(config);

    ASSERT_EQ(dataset.references.size(), 2);
    for (const auto& [name, seq] : dataset.references) {
        EXPECT_FALSE(name.empty());
        EXPECT_EQ(seq.size(), config.region_length);
        for (char c : seq) {
            EXPECT_TRUE(c == 'A' || c == 'C' || c == 'G' || c == 'T');
        }
    }
}

TEST_F(BenchmarkTruthTest, T1Deterministic) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 50;
    config.region_length = 5000;

    auto dataset1 = BenchmarkGenerator::generate_t1(config);
    auto dataset2 = BenchmarkGenerator::generate_t1(config);

    ASSERT_EQ(dataset1.reads.size(), dataset2.reads.size());
    for (size_t i = 0; i < dataset1.reads.size(); ++i) {
        EXPECT_EQ(dataset1.reads[i].id, dataset2.reads[i].id);
        EXPECT_EQ(dataset1.reads[i].sequence, dataset2.reads[i].sequence);
        EXPECT_EQ(dataset1.reads[i].chrom, dataset2.reads[i].chrom);
        EXPECT_EQ(dataset1.reads[i].pos, dataset2.reads[i].pos);
    }
}

// ============================================================================
// T2 Generation Tests
// ============================================================================

TEST_F(BenchmarkTruthTest, T2ConfigDefaults) {
    T2Config config;
    EXPECT_EQ(config.seed, 42);
    EXPECT_EQ(config.n_reads, 500000);
    EXPECT_EQ(config.coverage, 30);
    EXPECT_EQ(config.read_length, 10000);
    EXPECT_TRUE(config.families.empty());
}

TEST_F(BenchmarkTruthTest, T2GenerateWithIGH) {
    T2Config config;
    config.seed = 42;
    config.n_reads = 100;
    config.families.push_back(paralog_presets::igh_constant());

    auto dataset = BenchmarkGenerator::generate_t2(config);

    EXPECT_EQ(dataset.task, BenchmarkTask::T2_ParalogStress);
    EXPECT_GT(dataset.reads.size(), 0);
    EXPECT_EQ(dataset.references.size(), 5);  // IGHG1, IGHG2, IGHG3, IGHG4, IGHGP
}

TEST_F(BenchmarkTruthTest, T2GenerateWithMultipleFamilies) {
    T2Config config;
    config.seed = 42;
    config.n_reads = 150;
    config.families.push_back(paralog_presets::igh_constant());
    config.families.push_back(paralog_presets::nphp1());

    auto dataset = BenchmarkGenerator::generate_t2(config);

    EXPECT_EQ(dataset.references.size(), 7);  // 5 IGH + 2 NPHP1
    EXPECT_GT(dataset.reads.size(), 0);
}

TEST_F(BenchmarkTruthTest, T2ReadsHaveParalogTruth) {
    T2Config config;
    config.seed = 42;
    config.n_reads = 50;
    config.families.push_back(paralog_presets::igh_constant());

    auto dataset = BenchmarkGenerator::generate_t2(config);

    for (const auto& read : dataset.reads) {
        EXPECT_FALSE(read.id.empty());
        EXPECT_FALSE(read.sequence.empty());
        EXPECT_FALSE(read.chrom.empty());
        EXPECT_FALSE(read.paralog.empty());
        EXPECT_EQ(read.chrom, read.paralog);  // For T2, chrom == paralog name
    }
}

// ============================================================================
// Paralog Preset Tests
// ============================================================================

TEST_F(BenchmarkTruthTest, ParalogPresetIGH) {
    auto family = paralog_presets::igh_constant();
    EXPECT_EQ(family.name, "IGH_constant");
    EXPECT_EQ(family.members.size(), 5);
    EXPECT_EQ(family.members[0], "IGHG1");
    EXPECT_EQ(family.members[4], "IGHGP");
    EXPECT_EQ(family.region_length, 50000);
}

TEST_F(BenchmarkTruthTest, ParalogPresetNPHP1) {
    auto family = paralog_presets::nphp1();
    EXPECT_EQ(family.name, "NPHP1");
    EXPECT_EQ(family.members.size(), 2);
    EXPECT_EQ(family.members[0], "NPHP1");
    EXPECT_EQ(family.members[1], "NPHP1B");
}

TEST_F(BenchmarkTruthTest, ParalogPresetMHC) {
    auto family = paralog_presets::mhc_class1();
    EXPECT_EQ(family.name, "MHC_class1");
    EXPECT_EQ(family.members.size(), 2);
}

// ============================================================================
// File Output Tests
// ============================================================================

TEST_F(BenchmarkTruthTest, WriteFastq) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 10;
    config.region_length = 1000;

    auto dataset = BenchmarkGenerator::generate_t1(config);
    std::string path = test_dir_ / "reads.fastq";

    BenchmarkGenerator::write_fastq(dataset, path);

    EXPECT_TRUE(std::filesystem::exists(path));
    EXPECT_GT(std::filesystem::file_size(path), 0);

    // Verify FASTQ format.
    std::ifstream in(path);
    std::string line;
    int line_num = 0;
    while (std::getline(in, line)) {
        switch (line_num % 4) {
            case 0: EXPECT_EQ(line[0], '@'); break;
            case 2: EXPECT_EQ(line[0], '+'); break;
        }
        ++line_num;
    }
    EXPECT_EQ(line_num, static_cast<int>(dataset.reads.size() * 4));
}

TEST_F(BenchmarkTruthTest, WriteReference) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 10;
    config.region_length = 500;

    auto dataset = BenchmarkGenerator::generate_t1(config);
    std::string path = test_dir_ / "reference.fa";

    BenchmarkGenerator::write_reference(dataset, path);

    EXPECT_TRUE(std::filesystem::exists(path));
    EXPECT_GT(std::filesystem::file_size(path), 0);

    // Verify FASTA format.
    std::ifstream in(path);
    std::string line;
    int seq_count = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            ++seq_count;
        }
    }
    EXPECT_EQ(seq_count, static_cast<int>(dataset.references.size()));
}

TEST_F(BenchmarkTruthTest, WriteTruthTsv) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 20;
    config.region_length = 1000;

    auto dataset = BenchmarkGenerator::generate_t1(config);
    std::string path = test_dir_ / "truth.tsv";

    BenchmarkGenerator::write_truth_tsv(dataset, path);

    EXPECT_TRUE(std::filesystem::exists(path));

    // Verify TSV format: read_id<TAB>chrom<TAB>pos
    std::ifstream in(path);
    std::string line;
    int data_lines = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        ++data_lines;
        int tab_count = std::count(line.begin(), line.end(), '\t');
        EXPECT_EQ(tab_count, 2) << "Line: " << line;
    }
    EXPECT_EQ(data_lines, static_cast<int>(dataset.reads.size()));
}

TEST_F(BenchmarkTruthTest, WriteParalogTruthTsv) {
    T2Config config;
    config.seed = 42;
    config.n_reads = 30;
    config.families.push_back(paralog_presets::nphp1());

    auto dataset = BenchmarkGenerator::generate_t2(config);
    std::string path = test_dir_ / "truth_paralog.tsv";

    BenchmarkGenerator::write_paralog_truth_tsv(dataset, path);

    EXPECT_TRUE(std::filesystem::exists(path));

    // Verify TSV format: read_id<TAB>paralog<TAB>chrom<TAB>pos
    std::ifstream in(path);
    std::string line;
    int data_lines = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        ++data_lines;
        int tab_count = std::count(line.begin(), line.end(), '\t');
        EXPECT_EQ(tab_count, 3) << "Line: " << line;
    }
    EXPECT_EQ(data_lines, static_cast<int>(dataset.reads.size()));
}

TEST_F(BenchmarkTruthTest, TruthTsvComputeCompatible) {
    // Verify truth.tsv format matches benchmarks/metrics/compute.py expectations.
    T1Config config;
    config.seed = 42;
    config.n_reads = 5;
    config.region_length = 1000;

    auto dataset = BenchmarkGenerator::generate_t1(config);
    std::string path = test_dir_ / "truth.tsv";

    BenchmarkGenerator::write_truth_tsv(dataset, path);

    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string read_id, chrom;
        int64_t pos;
        std::getline(iss, read_id, '\t');
        std::getline(iss, chrom, '\t');
        iss >> pos;

        EXPECT_FALSE(read_id.empty());
        EXPECT_FALSE(chrom.empty());
        EXPECT_GE(pos, 0);
        EXPECT_LT(pos, static_cast<int64_t>(config.region_length));
    }
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_F(BenchmarkTruthTest, T1ZeroReads) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 0;
    config.region_length = 1000;

    auto dataset = BenchmarkGenerator::generate_t1(config);
    EXPECT_EQ(dataset.reads.size(), 0);
    EXPECT_EQ(dataset.total_reads, 0);
    EXPECT_EQ(dataset.total_bases, 0);
    EXPECT_EQ(dataset.references.size(), 2);  // References still generated
}

TEST_F(BenchmarkTruthTest, T2EmptyFamilies) {
    T2Config config;
    config.seed = 42;
    config.n_reads = 100;
    // No families added.

    auto dataset = BenchmarkGenerator::generate_t2(config);
    EXPECT_EQ(dataset.reads.size(), 0);
    EXPECT_EQ(dataset.references.size(), 0);
}

TEST_F(BenchmarkTruthTest, T1VeryShortRegion) {
    T1Config config;
    config.seed = 42;
    config.n_reads = 10;
    config.region_length = 200;  // Shorter than default read length
    config.read_length = 100;

    auto dataset = BenchmarkGenerator::generate_t1(config);
    EXPECT_GT(dataset.reads.size(), 0);
    for (const auto& read : dataset.reads) {
        EXPECT_LE(static_cast<int64_t>(read.sequence.size()),
                  static_cast<int64_t>(config.region_length));
    }
}

TEST_F(BenchmarkTruthTest, T1DifferentSeeds) {
    T1Config config1, config2;
    config1.seed = 1;
    config2.seed = 2;
    config1.n_reads = config2.n_reads = 20;
    config1.region_length = config2.region_length = 5000;

    auto dataset1 = BenchmarkGenerator::generate_t1(config1);
    auto dataset2 = BenchmarkGenerator::generate_t1(config2);

    // Different seeds should produce different data.
    bool all_same = true;
    for (size_t i = 0; i < std::min(dataset1.reads.size(), dataset2.reads.size()); ++i) {
        if (dataset1.reads[i].sequence != dataset2.reads[i].sequence) {
            all_same = false;
            break;
        }
    }
    EXPECT_FALSE(all_same);
}

}  // namespace
}  // namespace llmap::synthetic
