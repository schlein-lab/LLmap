// LLmap Phase 5.3 — Real reference validation tests.
//
// These tests verify the real reference validation infrastructure
// including BED parsing, BAM parsing, SLURM script generation, and
// the validation runner.

#include "validation/real_reference.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>

using namespace llmap::validation;
namespace fs = std::filesystem;

// =============================================================================
// Test fixtures
// =============================================================================

class RealReferenceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create temp directory for test files
        temp_dir_ = fs::temp_directory_path() / "llmap_test_real_ref";
        fs::create_directories(temp_dir_);
    }

    void TearDown() override {
        fs::remove_all(temp_dir_);
    }

    fs::path CreateTempFile(const std::string& name,
                            const std::string& content) {
        auto path = temp_dir_ / name;
        std::ofstream out(path);
        out << content;
        return path;
    }

    fs::path temp_dir_;
};

// =============================================================================
// RealReferenceConfig tests
// =============================================================================

TEST_F(RealReferenceTest, ConfigValidationRequiresReference) {
    RealReferenceConfig cfg;
    EXPECT_FALSE(cfg.Validate());
    EXPECT_NE(cfg.ValidationErrors().find("reference_fasta"),
              std::string::npos);
}

TEST_F(RealReferenceTest, ConfigValidationRequiresReads) {
    auto fasta = CreateTempFile("ref.fa", ">chr14\nACGT");
    RealReferenceConfig cfg;
    cfg.reference_fasta = fasta;

    EXPECT_FALSE(cfg.Validate());
    EXPECT_NE(cfg.ValidationErrors().find("reads_fastq"), std::string::npos);
}

TEST_F(RealReferenceTest, ConfigValidationPassesWithRequired) {
    auto fasta = CreateTempFile("ref.fa", ">chr14\nACGT");
    auto fastq = CreateTempFile("reads.fq", "@read1\nACGT\n+\nIIII");

    RealReferenceConfig cfg;
    cfg.reference_fasta = fasta;
    cfg.reads_fastq = fastq;

    EXPECT_TRUE(cfg.Validate());
    EXPECT_TRUE(cfg.ValidationErrors().empty());
}

TEST_F(RealReferenceTest, ConfigValidationChecksOptionalFiles) {
    auto fasta = CreateTempFile("ref.fa", ">chr14\nACGT");
    auto fastq = CreateTempFile("reads.fq", "@read1\nACGT\n+\nIIII");

    RealReferenceConfig cfg;
    cfg.reference_fasta = fasta;
    cfg.reads_fastq = fastq;
    cfg.minimap2_bam = "/nonexistent/baseline.bam";

    EXPECT_FALSE(cfg.Validate());
    EXPECT_NE(cfg.ValidationErrors().find("minimap2_bam"), std::string::npos);
}

TEST_F(RealReferenceTest, ConfigDefaultValues) {
    RealReferenceConfig cfg;

    EXPECT_EQ(cfg.minimizer_k, 15);
    EXPECT_EQ(cfg.minimizer_w, 10);
    EXPECT_FLOAT_EQ(cfg.min_identity, 0.70f);
    EXPECT_FALSE(cfg.use_gpu);
    EXPECT_EQ(cfg.gpu_device, 0);
}

// =============================================================================
// Ground truth BED parsing tests
// =============================================================================

TEST_F(RealReferenceTest, ParseGroundTruthBedBasic) {
    auto bed = CreateTempFile("truth.bed",
        "chr14\t1000\t2000\tread_1\t100\t+\n"
        "chr14\t3000\t4000\tread_2\t200\t-\n");

    auto truths = ParseGroundTruthBed(bed);

    ASSERT_EQ(truths.size(), 2u);

    EXPECT_EQ(truths[0].chrom, "chr14");
    EXPECT_EQ(truths[0].start, 1000u);
    EXPECT_EQ(truths[0].end, 2000u);
    EXPECT_EQ(truths[0].read_id, "read_1");
    EXPECT_EQ(truths[0].strand, "+");

    EXPECT_EQ(truths[1].chrom, "chr14");
    EXPECT_EQ(truths[1].start, 3000u);
    EXPECT_EQ(truths[1].end, 4000u);
    EXPECT_EQ(truths[1].read_id, "read_2");
    EXPECT_EQ(truths[1].strand, "-");
}

TEST_F(RealReferenceTest, ParseGroundTruthBedSkipsComments) {
    auto bed = CreateTempFile("truth.bed",
        "# This is a comment\n"
        "chr14\t1000\t2000\tread_1\t100\t+\n"
        "# Another comment\n"
        "chr14\t3000\t4000\tread_2\t200\t-\n");

    auto truths = ParseGroundTruthBed(bed);

    EXPECT_EQ(truths.size(), 2u);
}

TEST_F(RealReferenceTest, ParseGroundTruthBedSkipsEmptyLines) {
    auto bed = CreateTempFile("truth.bed",
        "chr14\t1000\t2000\tread_1\t100\t+\n"
        "\n"
        "chr14\t3000\t4000\tread_2\t200\t-\n");

    auto truths = ParseGroundTruthBed(bed);

    EXPECT_EQ(truths.size(), 2u);
}

TEST_F(RealReferenceTest, ParseGroundTruthBedSkipsMalformed) {
    auto bed = CreateTempFile("truth.bed",
        "chr14\t1000\n"  // Too few fields
        "chr14\t1000\t2000\tread_valid\t100\t+\n"
        "short\n");  // Too few fields

    auto truths = ParseGroundTruthBed(bed);

    EXPECT_EQ(truths.size(), 1u);
    EXPECT_EQ(truths[0].read_id, "read_valid");
}

TEST_F(RealReferenceTest, ParseGroundTruthBedNonexistent) {
    auto truths = ParseGroundTruthBed("/nonexistent/path.bed");
    EXPECT_TRUE(truths.empty());
}

TEST_F(RealReferenceTest, ParseGroundTruthBedWithGene) {
    auto bed = CreateTempFile("truth.bed",
        "chr14\t1000\t2000\tread_1\t100\t+\tIGHV3-23\n");

    auto truths = ParseGroundTruthBed(bed);

    ASSERT_EQ(truths.size(), 1u);
    EXPECT_EQ(truths[0].gene, "IGHV3-23");
}

// =============================================================================
// SLURM script generation tests
// =============================================================================

TEST_F(RealReferenceTest, GenerateSlurmScriptBasic) {
    RealReferenceConfig val_cfg;
    val_cfg.reference_fasta = "/data/ref.fa";
    val_cfg.reads_fastq = "/data/reads.fq";
    val_cfg.minimizer_k = 15;
    val_cfg.minimizer_w = 10;
    val_cfg.min_identity = 0.70f;

    SlurmJobConfig slurm_cfg;
    slurm_cfg.job_name = "test_job";
    slurm_cfg.partition = "gpu";
    slurm_cfg.n_gpus = 1;
    slurm_cfg.n_cpus = 8;
    slurm_cfg.memory = "32G";
    slurm_cfg.time_limit = "4:00:00";

    auto script = GenerateSlurmScript(val_cfg, slurm_cfg);

    EXPECT_NE(script.find("#!/bin/bash"), std::string::npos);
    EXPECT_NE(script.find("--job-name=test_job"), std::string::npos);
    EXPECT_NE(script.find("--partition=gpu"), std::string::npos);
    EXPECT_NE(script.find("--gres=gpu:1"), std::string::npos);
    EXPECT_NE(script.find("--cpus-per-task=8"), std::string::npos);
    EXPECT_NE(script.find("--mem=32G"), std::string::npos);
    EXPECT_NE(script.find("--time=4:00:00"), std::string::npos);
    EXPECT_NE(script.find("--reference /data/ref.fa"), std::string::npos);
    EXPECT_NE(script.find("--reads /data/reads.fq"), std::string::npos);
}

TEST_F(RealReferenceTest, GenerateSlurmScriptWithBaseline) {
    RealReferenceConfig val_cfg;
    val_cfg.reference_fasta = "/data/ref.fa";
    val_cfg.reads_fastq = "/data/reads.fq";
    val_cfg.minimap2_bam = "/data/baseline.bam";

    SlurmJobConfig slurm_cfg;
    slurm_cfg.job_name = "with_baseline";

    auto script = GenerateSlurmScript(val_cfg, slurm_cfg);

    EXPECT_NE(script.find("--baseline /data/baseline.bam"), std::string::npos);
}

TEST_F(RealReferenceTest, GenerateSlurmScriptWithGpu) {
    RealReferenceConfig val_cfg;
    val_cfg.reference_fasta = "/data/ref.fa";
    val_cfg.reads_fastq = "/data/reads.fq";
    val_cfg.use_gpu = true;
    val_cfg.gpu_device = 2;

    SlurmJobConfig slurm_cfg;
    slurm_cfg.job_name = "gpu_job";

    auto script = GenerateSlurmScript(val_cfg, slurm_cfg);

    EXPECT_NE(script.find("--gpu 2"), std::string::npos);
}

TEST_F(RealReferenceTest, GenerateSlurmScriptWithWorkDir) {
    RealReferenceConfig val_cfg;
    val_cfg.reference_fasta = "/data/ref.fa";
    val_cfg.reads_fastq = "/data/reads.fq";

    SlurmJobConfig slurm_cfg;
    slurm_cfg.job_name = "workdir_job";
    slurm_cfg.work_dir = "/work/llmap";

    auto script = GenerateSlurmScript(val_cfg, slurm_cfg);

    EXPECT_NE(script.find("cd /work/llmap"), std::string::npos);
}

TEST_F(RealReferenceTest, GenerateSlurmScriptModuleLoads) {
    RealReferenceConfig val_cfg;
    val_cfg.reference_fasta = "/data/ref.fa";
    val_cfg.reads_fastq = "/data/reads.fq";

    SlurmJobConfig slurm_cfg;

    auto script = GenerateSlurmScript(val_cfg, slurm_cfg);

    EXPECT_NE(script.find("module load cuda"), std::string::npos);
    EXPECT_NE(script.find("module load gcc"), std::string::npos);
}

// =============================================================================
// SlurmJobConfig tests
// =============================================================================

TEST(SlurmJobConfigTest, DefaultValues) {
    SlurmJobConfig cfg;

    EXPECT_EQ(cfg.job_name, "llmap_validation");
    EXPECT_EQ(cfg.partition, "gpu");
    EXPECT_EQ(cfg.n_gpus, 1);
    EXPECT_EQ(cfg.n_cpus, 8);
    EXPECT_EQ(cfg.memory, "32G");
    EXPECT_EQ(cfg.time_limit, "4:00:00");
}

// =============================================================================
// RealReferenceResult tests
// =============================================================================

TEST(RealReferenceResultTest, SummaryContainsAllSections) {
    RealReferenceResult result;
    result.n_reference_seqs = 1;
    result.reference_total_bp = 1000000;
    result.n_reads = 100;
    result.llmap_mapped = 95;
    result.llmap_unmapped = 5;
    result.llmap_alignment_rate = 0.95f;
    result.has_minimap2_baseline = true;
    result.minimap2_mapped = 93;
    result.baseline_comparison.both_mapped = 90;
    result.baseline_comparison.llmap_only = 5;
    result.baseline_comparison.minimap2_only = 3;
    result.baseline_comparison.recall_ratio = 1.02f;
    result.index_load_time_ms = 100.0f;
    result.alignment_time_ms = 500.0f;
    result.validation_time_ms = 50.0f;
    result.total_time_ms = 650.0f;
    result.recall_pass = true;
    result.lossless_pass = true;
    result.overall_pass = true;

    std::string summary = result.Summary();

    EXPECT_NE(summary.find("Reference:"), std::string::npos);
    EXPECT_NE(summary.find("Reads:"), std::string::npos);
    EXPECT_NE(summary.find("LLmap Results:"), std::string::npos);
    EXPECT_NE(summary.find("Minimap2 Baseline:"), std::string::npos);
    EXPECT_NE(summary.find("Timing:"), std::string::npos);
    EXPECT_NE(summary.find("Kill-Switch Verdicts:"), std::string::npos);
    EXPECT_NE(summary.find("PASS"), std::string::npos);
}

TEST(RealReferenceResultTest, SummaryShowsFailure) {
    RealReferenceResult result;
    result.recall_pass = false;
    result.lossless_pass = true;
    result.overall_pass = false;
    result.verdict_reason = "Recall below threshold";

    std::string summary = result.Summary();

    EXPECT_NE(summary.find("FAIL"), std::string::npos);
    EXPECT_NE(summary.find("Recall below threshold"), std::string::npos);
}

TEST(RealReferenceResultTest, SummaryHidesBaselineIfNotPresent) {
    RealReferenceResult result;
    result.has_minimap2_baseline = false;

    std::string summary = result.Summary();

    EXPECT_EQ(summary.find("Minimap2 Baseline:"), std::string::npos);
}

TEST(RealReferenceResultTest, SummaryHidesGroundTruthIfNotPresent) {
    RealReferenceResult result;
    result.has_ground_truth = false;

    std::string summary = result.Summary();

    EXPECT_EQ(summary.find("Position Accuracy:"), std::string::npos);
}

// =============================================================================
// RunRealReferenceValidation tests
// =============================================================================

TEST_F(RealReferenceTest, ValidationFailsWithInvalidConfig) {
    RealReferenceConfig cfg;
    // Empty config - should fail validation

    auto result = RunRealReferenceValidation(cfg);

    EXPECT_FALSE(result.overall_pass);
    EXPECT_NE(result.verdict_reason.find("Invalid config"), std::string::npos);
}

TEST_F(RealReferenceTest, ValidationRunsWithMinimalData) {
    // Create minimal reference
    auto fasta = CreateTempFile("ref.fa",
        ">synthetic_locus\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");

    // Create minimal reads matching the reference
    auto fastq = CreateTempFile("reads.fq",
        "@read_1\n"
        "ACGTACGTACGTACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIII\n"
        "@read_2\n"
        "GTACGTACGTACGTACGTACGTAC\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIII\n");

    RealReferenceConfig cfg;
    cfg.reference_fasta = fasta;
    cfg.reads_fastq = fastq;
    cfg.minimizer_k = 7;  // Smaller k for short sequences
    cfg.minimizer_w = 5;

    auto result = RunRealReferenceValidation(cfg);

    // Should complete without error
    EXPECT_EQ(result.n_reference_seqs, 1u);
    EXPECT_EQ(result.n_reads, 2u);
    EXPECT_GT(result.total_time_ms, 0.0f);
}

TEST_F(RealReferenceTest, ValidationReportsLosslessCorrectly) {
    auto fasta = CreateTempFile("ref.fa",
        ">locus\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");

    auto fastq = CreateTempFile("reads.fq",
        "@read_1\n"
        "ACGTACGTACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIII\n");

    RealReferenceConfig cfg;
    cfg.reference_fasta = fasta;
    cfg.reads_fastq = fastq;
    cfg.minimizer_k = 7;
    cfg.minimizer_w = 5;

    auto result = RunRealReferenceValidation(cfg);

    // Lossless: 1 read in, should produce 1 record (mapped or unmapped)
    EXPECT_EQ(result.llmap_mapped + result.llmap_unmapped, result.n_reads);
    EXPECT_TRUE(result.lossless_pass);
}

TEST_F(RealReferenceTest, ValidationWithGroundTruth) {
    auto fasta = CreateTempFile("ref.fa",
        ">chr14\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");

    auto fastq = CreateTempFile("reads.fq",
        "@read_1\n"
        "ACGTACGTACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIII\n");

    auto bed = CreateTempFile("truth.bed",
        "chr14\t0\t20\tread_1\t100\t+\n");

    RealReferenceConfig cfg;
    cfg.reference_fasta = fasta;
    cfg.reads_fastq = fastq;
    cfg.ground_truth_bed = bed;
    cfg.minimizer_k = 7;
    cfg.minimizer_w = 5;

    auto result = RunRealReferenceValidation(cfg);

    // Ground truth should be loaded
    EXPECT_TRUE(result.has_ground_truth);
}

// =============================================================================
// ParseSlurmOutput tests
// =============================================================================

TEST_F(RealReferenceTest, ParseSlurmOutputBasic) {
    auto output = CreateTempFile("slurm_output.log",
        "=== Real Reference Validation ===\n"
        "\n"
        "Reference:\n"
        "  Sequences:    5\n"
        "  Total bp:     1000000\n"
        "\n"
        "LLmap Results:\n"
        "  Mapped:       950\n"
        "  Unmapped:     50\n"
        "\n"
        "Overall:                 PASS\n");

    auto result = ParseSlurmOutput(output);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->n_reference_seqs, 5u);
    EXPECT_EQ(result->reference_total_bp, 1000000u);
    EXPECT_EQ(result->llmap_mapped, 950u);
    EXPECT_TRUE(result->overall_pass);
}

TEST_F(RealReferenceTest, ParseSlurmOutputFailed) {
    auto output = CreateTempFile("slurm_output.log",
        "=== Real Reference Validation ===\n"
        "\n"
        "Overall:                 FAIL (Recall below threshold)\n");

    auto result = ParseSlurmOutput(output);

    ASSERT_TRUE(result.has_value());
    EXPECT_FALSE(result->overall_pass);
}

TEST_F(RealReferenceTest, ParseSlurmOutputNonexistent) {
    auto result = ParseSlurmOutput("/nonexistent/output.log");
    EXPECT_FALSE(result.has_value());
}

// =============================================================================
// RealGroundTruth tests
// =============================================================================

TEST(RealGroundTruthTest, DefaultValues) {
    RealGroundTruth truth;

    EXPECT_TRUE(truth.read_id.empty());
    EXPECT_TRUE(truth.chrom.empty());
    EXPECT_EQ(truth.start, 0u);
    EXPECT_EQ(truth.end, 0u);
    EXPECT_TRUE(truth.strand.empty());
    EXPECT_TRUE(truth.gene.empty());
}

// =============================================================================
// SlurmJobStatus tests
// =============================================================================

TEST(SlurmJobStatusTest, DefaultValues) {
    SlurmJobStatus status;

    EXPECT_TRUE(status.job_id.empty());
    EXPECT_TRUE(status.state.empty());
    EXPECT_EQ(status.exit_code, -1);
    EXPECT_FALSE(status.is_complete);
}
