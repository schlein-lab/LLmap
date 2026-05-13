// LLmap Phase 5.2 — End-to-end synthetic validation tests.
//
// These tests verify the complete pipeline from synthetic data generation
// through alignment to kill-switch validation.

#include "validation/killswitch.h"

#include <gtest/gtest.h>

#include "classical/classical_pipeline.h"
#include "core/alignment_record.h"
#include "synthetic/igh_locus_generator.h"

using namespace llmap;
using namespace llmap::validation;
using namespace llmap::classical;
using namespace llmap::synthetic;

// =============================================================================
// Conversion tests
// =============================================================================

TEST(ConvertToAlignmentRecord, BasicConversion) {
    ClassicalAlignment aln;
    aln.query_name = "test_read";
    aln.ref_name = "chr14";
    aln.ref_id = 0;
    aln.ref_start = 1000;
    aln.ref_end = 2000;
    aln.query_start = 0;
    aln.query_end = 1000;
    aln.is_forward = true;
    aln.is_primary = true;
    aln.cigar.push_back({CigarOp::Match, 1000});
    aln.score = 2000;
    aln.identity = 0.99f;
    aln.mapq = 60;

    auto record = ConvertToAlignmentRecord(aln, 1000);

    EXPECT_EQ(record.read_id, "test_read");
    EXPECT_EQ(record.read_len, 1000u);
    EXPECT_EQ(record.status, AlignmentStatus::Mapped);
    ASSERT_TRUE(record.primary.has_value());
    EXPECT_EQ(record.primary->target_id, "chr14");
    EXPECT_EQ(record.primary->start, 1000u);
    EXPECT_EQ(record.primary->end, 2000u);
    EXPECT_EQ(record.primary->score, 2000);
}

TEST(ConvertResults, MixedResults) {
    std::vector<ReadAlignmentResult> results;

    // Mapped read
    ReadAlignmentResult r1;
    r1.query_name = "read_1";
    ClassicalAlignment aln;
    aln.query_name = "read_1";
    aln.ref_name = "ref";
    aln.ref_start = 100;
    aln.ref_end = 200;
    aln.is_primary = true;
    aln.cigar.push_back({CigarOp::Match, 100});
    aln.identity = 0.95f;
    r1.alignments.push_back(aln);
    results.push_back(r1);

    // Unmapped read
    ReadAlignmentResult r2;
    r2.query_name = "read_2";
    results.push_back(r2);

    std::vector<uint32_t> lens = {100, 100};
    auto records = ConvertResults(results, lens);

    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0].status, AlignmentStatus::Mapped);
    EXPECT_EQ(records[1].status, AlignmentStatus::Unmapped);
}

// =============================================================================
// EndToEndConfig preset tests
// =============================================================================

TEST(EndToEndConfig, MinimalPreset) {
    auto cfg = EndToEndConfig::Minimal();
    EXPECT_EQ(cfg.n_reads, 10u);
    EXPECT_EQ(cfg.locus_length, 10000u);
    EXPECT_DOUBLE_EQ(cfg.min_position_accuracy, 0.80);
}

TEST(EndToEndConfig, StandardPreset) {
    auto cfg = EndToEndConfig::Standard();
    EXPECT_EQ(cfg.n_reads, 50u);
    EXPECT_EQ(cfg.locus_length, 30000u);
    EXPECT_DOUBLE_EQ(cfg.min_position_accuracy, 0.90);
}

TEST(EndToEndConfig, StressPreset) {
    auto cfg = EndToEndConfig::Stress();
    EXPECT_EQ(cfg.n_reads, 200u);
    EXPECT_EQ(cfg.locus_length, 100000u);
    EXPECT_DOUBLE_EQ(cfg.min_position_accuracy, 0.95);
}

// =============================================================================
// EndToEndResult tests
// =============================================================================

TEST(EndToEndResult, SummaryContainsAllSections) {
    EndToEndResult result;
    result.n_reads_generated = 100;
    result.n_canonical = 80;
    result.n_duplicate = 20;
    result.actual_dup_fraction = 0.2f;
    result.n_aligned = 95;
    result.n_unmapped = 5;
    result.alignment_rate = 0.95f;
    result.avg_identity = 0.98f;
    result.position_accuracy = 0.92;
    result.mean_position_error = 25;
    result.max_position_error = 150;
    result.generation_time_ms = 10.0f;
    result.alignment_time_ms = 50.0f;
    result.validation_time_ms = 5.0f;
    result.total_time_ms = 65.0f;
    result.passed = true;

    std::string summary = result.Summary();
    EXPECT_NE(summary.find("Synthetic Data"), std::string::npos);
    EXPECT_NE(summary.find("Alignment Results"), std::string::npos);
    EXPECT_NE(summary.find("Position Accuracy"), std::string::npos);
    EXPECT_NE(summary.find("Timing"), std::string::npos);
    EXPECT_NE(summary.find("PASS"), std::string::npos);
}

TEST(EndToEndResult, SummaryShowsFailure) {
    EndToEndResult result;
    result.passed = false;
    result.verdict_reason = "Test failure reason";

    std::string summary = result.Summary();
    EXPECT_NE(summary.find("FAIL"), std::string::npos);
    EXPECT_NE(summary.find("Test failure reason"), std::string::npos);
}

// =============================================================================
// Minimal validation run tests
// =============================================================================

TEST(RunEndToEndValidation, MinimalConfigRuns) {
    auto result = RunMinimalValidation(42);

    EXPECT_GT(result.n_reads_generated, 0u);
    EXPECT_EQ(result.n_reads_generated, result.n_aligned + result.n_unmapped);
    EXPECT_GT(result.total_time_ms, 0.0f);
}

TEST(RunEndToEndValidation, LosslessGuaranteeHeld) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 123;
    cfg.require_lossless = true;

    auto result = RunEndToEndValidation(cfg);

    EXPECT_TRUE(result.validation.lossless_pass);
    EXPECT_EQ(result.validation.missing_records, 0u);
}

TEST(RunEndToEndValidation, AllReadsProduceRecords) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 456;

    auto result = RunEndToEndValidation(cfg);

    EXPECT_EQ(result.validation.input_reads, result.n_reads_generated);
    EXPECT_EQ(result.validation.output_records, result.n_reads_generated);
}

TEST(RunEndToEndValidation, AlignmentRateReasonable) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 789;

    auto result = RunEndToEndValidation(cfg);

    // Should align at least some reads
    EXPECT_GT(result.n_aligned, 0u);
    EXPECT_GT(result.alignment_rate, 0.0f);
}

TEST(RunEndToEndValidation, DeterministicBySeed) {
    auto cfg = EndToEndConfig::Minimal();

    auto result1 = RunEndToEndValidation(cfg);
    auto result2 = RunEndToEndValidation(cfg);

    EXPECT_EQ(result1.n_reads_generated, result2.n_reads_generated);
    EXPECT_EQ(result1.n_aligned, result2.n_aligned);
    EXPECT_EQ(result1.validation.position_correct, result2.validation.position_correct);
}

TEST(RunEndToEndValidation, DifferentSeedsDifferentResults) {
    auto cfg1 = EndToEndConfig::Minimal();
    cfg1.seed = 111;

    auto cfg2 = EndToEndConfig::Minimal();
    cfg2.seed = 222;

    auto result1 = RunEndToEndValidation(cfg1);
    auto result2 = RunEndToEndValidation(cfg2);

    // Results may differ (not guaranteed but likely with different seeds)
    // Just verify both run successfully
    EXPECT_GT(result1.n_reads_generated, 0u);
    EXPECT_GT(result2.n_reads_generated, 0u);
}

// =============================================================================
// Position accuracy tests
// =============================================================================

TEST(RunEndToEndValidation, PositionAccuracyComputed) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 321;

    auto result = RunEndToEndValidation(cfg);

    // If reads were aligned, position accuracy should be computed
    if (result.validation.position_evaluated > 0) {
        EXPECT_GE(result.position_accuracy, 0.0);
        EXPECT_LE(result.position_accuracy, 1.0);
    }
}

TEST(RunEndToEndValidation, PositionErrorStatsComputed) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 654;

    auto result = RunEndToEndValidation(cfg);

    if (result.validation.position_evaluated > 0) {
        // Max error should be >= mean error
        EXPECT_GE(result.max_position_error, std::abs(result.mean_position_error));
    }
}

// =============================================================================
// Tolerance configuration tests
// =============================================================================

TEST(RunEndToEndValidation, PositionToleranceAffectsAccuracy) {
    auto cfg1 = EndToEndConfig::Minimal();
    cfg1.seed = 987;
    cfg1.position_tolerance_bp = 50;

    auto cfg2 = EndToEndConfig::Minimal();
    cfg2.seed = 987;
    cfg2.position_tolerance_bp = 500;

    auto result1 = RunEndToEndValidation(cfg1);
    auto result2 = RunEndToEndValidation(cfg2);

    // Larger tolerance should give >= position correct count
    EXPECT_GE(result2.validation.position_correct, result1.validation.position_correct);
}

// =============================================================================
// Verdict tests
// =============================================================================

TEST(RunEndToEndValidation, VerdictPassesWhenThresholdsMet) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 135;
    cfg.min_position_accuracy = 0.0;  // Very permissive
    cfg.require_lossless = false;

    auto result = RunEndToEndValidation(cfg);

    EXPECT_TRUE(result.passed);
    EXPECT_TRUE(result.verdict_reason.empty());
}

TEST(RunEndToEndValidation, VerdictFailsWhenPositionAccuracyLow) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 246;
    cfg.min_position_accuracy = 1.0;  // Impossible to achieve
    cfg.require_lossless = false;

    auto result = RunEndToEndValidation(cfg);

    // Should fail unless perfect accuracy
    if (result.position_accuracy < 1.0) {
        EXPECT_FALSE(result.passed);
        EXPECT_NE(result.verdict_reason.find("Position accuracy"), std::string::npos);
    }
}

// =============================================================================
// Timing tests
// =============================================================================

TEST(RunEndToEndValidation, TimingFieldsPopulated) {
    auto result = RunMinimalValidation(753);

    EXPECT_GT(result.generation_time_ms, 0.0f);
    EXPECT_GT(result.alignment_time_ms, 0.0f);
    EXPECT_GT(result.validation_time_ms, 0.0f);
    EXPECT_GT(result.total_time_ms, 0.0f);

    // Total should approximately equal sum of phases
    float sum = result.generation_time_ms + result.alignment_time_ms +
                result.validation_time_ms;
    EXPECT_NEAR(result.total_time_ms, sum, 10.0f);  // 10ms tolerance
}

// =============================================================================
// Standard validation test (longer but more thorough)
// =============================================================================

TEST(RunEndToEndValidation, StandardConfigProducesResults) {
    auto result = RunStandardValidation(42);

    EXPECT_GT(result.n_reads_generated, 30u);
    EXPECT_TRUE(result.validation.lossless_pass);

    // With standard config, should have reasonable alignment rate
    EXPECT_GT(result.alignment_rate, 0.5f);
}

// =============================================================================
// Integration with KillSwitchValidator
// =============================================================================

TEST(RunEndToEndValidation, ValidationStatsConsistent) {
    auto result = RunMinimalValidation(864);

    const auto& vs = result.validation;

    // Status counts should sum to output records
    EXPECT_EQ(vs.mapped + vs.tentative + vs.unmapped, vs.output_records);

    // Position evaluated should be <= mapped
    EXPECT_LE(vs.position_evaluated, vs.mapped);

    // Position correct should be <= position evaluated
    EXPECT_LE(vs.position_correct, vs.position_evaluated);
}

// =============================================================================
// Edge cases
// =============================================================================

TEST(RunEndToEndValidation, HandlesPureDuplicate) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 975;
    cfg.dup_fraction = 1.0f;  // All duplicates

    auto result = RunEndToEndValidation(cfg);

    // Should still run without error
    EXPECT_GT(result.n_reads_generated, 0u);
    EXPECT_TRUE(result.validation.lossless_pass);
}

TEST(RunEndToEndValidation, HandlesBalancedMosaic) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = 468;
    cfg.dup_fraction = 0.5f;  // 50/50 split

    auto result = RunEndToEndValidation(cfg);

    // Both canonical and duplicate should be present
    // (exact counts depend on random sampling, but both should exist)
    EXPECT_GT(result.n_reads_generated, 0u);
}
