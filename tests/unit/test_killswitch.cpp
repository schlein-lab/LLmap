// LLmap Phase 5.1 — Kill-switch validation framework tests.
//
// These tests verify the validation infrastructure used for Phase 5 go/no-go.

#include "validation/killswitch.h"

#include <gtest/gtest.h>

#include "core/alignment_record.h"
#include "synthetic/igh_locus_generator.h"

using namespace llmap;
using namespace llmap::validation;

namespace {

AlignmentHit MakeHit(uint64_t start, uint64_t end) {
    AlignmentHit h;
    h.target_id = "synthetic_locus";
    h.start = start;
    h.end = end;
    h.cigar.ops = std::to_string(end - start) + "M";
    h.score = static_cast<int32_t>((end - start) * 2);
    h.nm = 5;
    return h;
}

ParalogCall MakeParalogCall(bool is_canonical) {
    ParalogCall pc;
    if (is_canonical) {
        pc.p_canonical = 0.95f;
        pc.p_dup = 0.05f;
    } else {
        pc.p_canonical = 0.05f;
        pc.p_dup = 0.95f;
    }
    pc.n_discriminating_psvs = 5;
    return pc;
}

}  // namespace

// =============================================================================
// GroundTruth parsing tests
// =============================================================================

TEST(GroundTruth, ParseValidReadId) {
    auto truth = GroundTruth::ParseFromReadId("synth_42_CAN_1000-2000_IGHG1");
    ASSERT_TRUE(truth.has_value());
    EXPECT_EQ(truth->read_id, "synth_42_CAN_1000-2000_IGHG1");
    EXPECT_EQ(truth->origin, "CAN");
    EXPECT_EQ(truth->true_start, 1000u);
    EXPECT_EQ(truth->true_end, 2000u);
    EXPECT_EQ(truth->source_gene, "IGHG1");
}

TEST(GroundTruth, ParseDuplicateOrigin) {
    auto truth = GroundTruth::ParseFromReadId("synth_7_DUP_500-1500_IGHG4");
    ASSERT_TRUE(truth.has_value());
    EXPECT_EQ(truth->origin, "DUP");
    EXPECT_EQ(truth->true_start, 500u);
    EXPECT_EQ(truth->true_end, 1500u);
}

TEST(GroundTruth, ParseInvalidReadId) {
    EXPECT_FALSE(GroundTruth::ParseFromReadId("invalid_read_id").has_value());
    EXPECT_FALSE(GroundTruth::ParseFromReadId("read_123").has_value());
    EXPECT_FALSE(GroundTruth::ParseFromReadId("").has_value());
}

TEST(GroundTruth, ParseMinimalParts) {
    auto truth = GroundTruth::ParseFromReadId("synth_0_CAN_100-200");
    ASSERT_TRUE(truth.has_value());
    EXPECT_EQ(truth->origin, "CAN");
    EXPECT_EQ(truth->true_start, 100u);
    EXPECT_TRUE(truth->source_gene.empty());
}

// =============================================================================
// KillSwitchValidator basic tests
// =============================================================================

TEST(KillSwitchValidator, EmptyValidatorPassesLossless) {
    KillSwitchValidator validator;
    auto stats = validator.Validate();
    EXPECT_EQ(stats.input_reads, 0u);
    EXPECT_EQ(stats.output_records, 0u);
    EXPECT_TRUE(stats.lossless_pass);
}

TEST(KillSwitchValidator, MissingRecordFailsLossless) {
    KillSwitchValidator validator;

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 100, 200, "IGHG1", 2});
    truths.push_back({"read_2", "DUP", 300, 400, "IGHG2", 1});
    validator.LoadGroundTruth(truths);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.input_reads, 2u);
    EXPECT_EQ(stats.output_records, 0u);
    EXPECT_EQ(stats.missing_records, 2u);
    EXPECT_FALSE(stats.lossless_pass);
}

TEST(KillSwitchValidator, MatchingRecordsPassLossless) {
    KillSwitchValidator validator;

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 100, 200, "IGHG1", 2});
    truths.push_back({"read_2", "DUP", 300, 400, "IGHG2", 1});
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;
    records.push_back(make_mapped("read_1", 100, MakeHit(100, 200)));
    records.push_back(make_mapped("read_2", 100, MakeHit(300, 400)));
    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.input_reads, 2u);
    EXPECT_EQ(stats.output_records, 2u);
    EXPECT_EQ(stats.missing_records, 0u);
    EXPECT_TRUE(stats.lossless_pass);
}

TEST(KillSwitchValidator, StatusBreakdownCorrect) {
    KillSwitchValidator validator;

    std::vector<GroundTruth> truths;
    for (int i = 0; i < 10; ++i) {
        truths.push_back({"read_" + std::to_string(i), "CAN", 100, 200, "IGHG1", 0});
    }
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;
    for (int i = 0; i < 5; ++i) {
        records.push_back(make_mapped("read_" + std::to_string(i), 100, MakeHit(100, 200)));
    }
    for (int i = 5; i < 8; ++i) {
        TentativeTarget tt;
        tt.target_id = "chr14";
        tt.final_probability = 0.4f;
        records.push_back(make_tentative(
            "read_" + std::to_string(i), 100, {tt}, RejectionReason::DidNotConverge));
    }
    for (int i = 8; i < 10; ++i) {
        records.push_back(make_unmapped("read_" + std::to_string(i), 100, RejectionReason::NoSeeds));
    }
    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.mapped, 5u);
    EXPECT_EQ(stats.tentative, 3u);
    EXPECT_EQ(stats.unmapped, 2u);
    EXPECT_TRUE(stats.lossless_pass);
}

// =============================================================================
// Position accuracy tests
// =============================================================================

TEST(KillSwitchValidator, CorrectPositionWithinTolerance) {
    KillSwitchValidator validator;

    PositionTolerance tol;
    tol.base_tolerance = 50;
    validator.SetPositionTolerance(tol);

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 1000, 2000, "IGHG1", 0});
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;
    records.push_back(make_mapped("read_1", 1000, MakeHit(1020, 2020)));
    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.position_evaluated, 1u);
    EXPECT_EQ(stats.position_correct, 1u);
    EXPECT_DOUBLE_EQ(stats.position_accuracy, 1.0);
}

TEST(KillSwitchValidator, IncorrectPositionOutsideTolerance) {
    KillSwitchValidator validator;

    PositionTolerance tol;
    tol.base_tolerance = 50;
    validator.SetPositionTolerance(tol);

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 1000, 2000, "IGHG1", 0});
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;
    records.push_back(make_mapped("read_1", 1000, MakeHit(1100, 2100)));
    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.position_evaluated, 1u);
    EXPECT_EQ(stats.position_correct, 0u);
    EXPECT_DOUBLE_EQ(stats.position_accuracy, 0.0);
}

TEST(KillSwitchValidator, RelativeToleranceUsed) {
    KillSwitchValidator validator;

    PositionTolerance tol;
    tol.base_tolerance = 10;
    tol.relative_tolerance = 0.05f;
    validator.SetPositionTolerance(tol);

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 0, 1000, "IGHG1", 0});
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;
    records.push_back(make_mapped("read_1", 1000, MakeHit(40, 1040)));
    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.position_evaluated, 1u);
    EXPECT_EQ(stats.position_correct, 1u);
}

// =============================================================================
// Origin accuracy tests
// =============================================================================

TEST(KillSwitchValidator, OriginAccuracyComputed) {
    KillSwitchValidator validator;

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 100, 200, "IGHG1", 0});
    truths.push_back({"read_2", "DUP", 300, 400, "IGHG2", 0});
    truths.push_back({"read_3", "CAN", 500, 600, "IGHG3", 0});
    truths.push_back({"read_4", "DUP", 700, 800, "IGHG4", 0});
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;

    auto r1 = make_mapped("read_1", 100, MakeHit(100, 200));
    r1.paralog_assignment = MakeParalogCall(true);
    records.push_back(r1);

    auto r2 = make_mapped("read_2", 100, MakeHit(300, 400));
    r2.paralog_assignment = MakeParalogCall(false);
    records.push_back(r2);

    auto r3 = make_mapped("read_3", 100, MakeHit(500, 600));
    r3.paralog_assignment = MakeParalogCall(false);
    records.push_back(r3);

    auto r4 = make_mapped("read_4", 100, MakeHit(700, 800));
    r4.paralog_assignment = MakeParalogCall(true);
    records.push_back(r4);

    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.origin_evaluated, 4u);
    EXPECT_EQ(stats.origin_correct, 2u);
    EXPECT_DOUBLE_EQ(stats.origin_accuracy, 0.5);
}

// =============================================================================
// Consistency tests
// =============================================================================

TEST(KillSwitchValidator, InconsistentRecordFailsKillSwitch) {
    KillSwitchValidator validator;

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 100, 200, "IGHG1", 0});
    validator.LoadGroundTruth(truths);

    AlignmentRecord bad_record;
    bad_record.read_id = "read_1";
    bad_record.status = AlignmentStatus::Mapped;

    std::vector<AlignmentRecord> records;
    records.push_back(bad_record);
    validator.AddRecords(records);

    auto stats = validator.Validate();
    EXPECT_EQ(stats.lossless_consistent, 0u);
    EXPECT_EQ(stats.lossless_inconsistent, 1u);
    EXPECT_FALSE(stats.lossless_pass);
}

// =============================================================================
// Synthetic data integration tests
// =============================================================================

TEST(KillSwitchValidator, LoadGroundTruthFromDataset) {
    auto cfg = synthetic::presets::tiny_test(42);
    synthetic::IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    KillSwitchValidator validator;
    validator.LoadGroundTruth(dataset);

    EXPECT_EQ(validator.NumGroundTruth(), dataset.reads.size());
}

TEST(KillSwitchValidator, ValidateSyntheticRunConvenience) {
    auto cfg = synthetic::presets::tiny_test(42);
    synthetic::IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    std::vector<AlignmentRecord> records;
    for (const auto& read : dataset.reads) {
        records.push_back(make_mapped(
            read.id, read.sequence.size(), MakeHit(read.true_start, read.true_end)));
    }

    auto stats = ValidateSyntheticRun(dataset, records);
    EXPECT_TRUE(stats.lossless_pass);
    EXPECT_EQ(stats.missing_records, 0u);
    EXPECT_EQ(stats.mapped, dataset.reads.size());
}

// =============================================================================
// Baseline comparison tests
// =============================================================================

TEST(BaselineComparison, BothMappersMapped) {
    std::vector<AlignmentRecord> llmap_records;
    llmap_records.push_back(make_mapped("read_1", 100, MakeHit(100, 200)));
    llmap_records.push_back(make_mapped("read_2", 100, MakeHit(200, 300)));

    std::unordered_map<std::string, bool> minimap2_mapped;
    minimap2_mapped["read_1"] = true;
    minimap2_mapped["read_2"] = true;

    auto cmp = CompareToBaseline(llmap_records, minimap2_mapped);
    EXPECT_EQ(cmp.llmap_mapped, 2u);
    EXPECT_EQ(cmp.minimap2_mapped, 2u);
    EXPECT_EQ(cmp.both_mapped, 2u);
    EXPECT_EQ(cmp.llmap_only, 0u);
    EXPECT_EQ(cmp.minimap2_only, 0u);
    EXPECT_DOUBLE_EQ(cmp.recall_ratio, 1.0);
}

TEST(BaselineComparison, LLmapBetter) {
    std::vector<AlignmentRecord> llmap_records;
    llmap_records.push_back(make_mapped("read_1", 100, MakeHit(100, 200)));
    llmap_records.push_back(make_mapped("read_2", 100, MakeHit(200, 300)));
    llmap_records.push_back(make_mapped("read_3", 100, MakeHit(300, 400)));

    std::unordered_map<std::string, bool> minimap2_mapped;
    minimap2_mapped["read_1"] = true;
    minimap2_mapped["read_2"] = true;
    minimap2_mapped["read_3"] = false;

    auto cmp = CompareToBaseline(llmap_records, minimap2_mapped);
    EXPECT_EQ(cmp.llmap_mapped, 3u);
    EXPECT_EQ(cmp.minimap2_mapped, 2u);
    EXPECT_EQ(cmp.both_mapped, 2u);
    EXPECT_EQ(cmp.llmap_only, 1u);
    EXPECT_EQ(cmp.minimap2_only, 0u);
    EXPECT_GT(cmp.recall_ratio, 1.0);
}

TEST(BaselineComparison, Minimap2Better) {
    std::vector<AlignmentRecord> llmap_records;
    llmap_records.push_back(make_mapped("read_1", 100, MakeHit(100, 200)));
    llmap_records.push_back(make_unmapped("read_2", 100, RejectionReason::DidNotConverge));
    llmap_records.push_back(make_unmapped("read_3", 100, RejectionReason::NoSeeds));

    std::unordered_map<std::string, bool> minimap2_mapped;
    minimap2_mapped["read_1"] = true;
    minimap2_mapped["read_2"] = true;
    minimap2_mapped["read_3"] = true;

    auto cmp = CompareToBaseline(llmap_records, minimap2_mapped);
    EXPECT_EQ(cmp.llmap_mapped, 1u);
    EXPECT_EQ(cmp.minimap2_mapped, 3u);
    EXPECT_EQ(cmp.both_mapped, 1u);
    EXPECT_EQ(cmp.llmap_only, 0u);
    EXPECT_EQ(cmp.minimap2_only, 2u);
    EXPECT_LT(cmp.recall_ratio, 1.0);
}

// =============================================================================
// Detailed validation tests
// =============================================================================

TEST(KillSwitchValidator, DetailedValidationReturnsAllReads) {
    KillSwitchValidator validator;

    std::vector<GroundTruth> truths;
    truths.push_back({"read_1", "CAN", 100, 200, "IGHG1", 0});
    truths.push_back({"read_2", "DUP", 300, 400, "IGHG2", 0});
    validator.LoadGroundTruth(truths);

    std::vector<AlignmentRecord> records;
    records.push_back(make_mapped("read_1", 100, MakeHit(100, 200)));
    validator.AddRecords(records);

    auto detailed = validator.ValidateDetailed();
    EXPECT_EQ(detailed.size(), 2u);

    int found = 0;
    for (const auto& v : detailed) {
        if (v.read_id == "read_1") {
            EXPECT_TRUE(v.has_ground_truth);
            EXPECT_TRUE(v.has_alignment_record);
            ++found;
        } else if (v.read_id == "read_2") {
            EXPECT_TRUE(v.has_ground_truth);
            EXPECT_FALSE(v.has_alignment_record);
            ++found;
        }
    }
    EXPECT_EQ(found, 2);
}

// =============================================================================
// Summary output tests
// =============================================================================

TEST(ValidationStats, SummaryContainsKeyInformation) {
    ValidationStats stats;
    stats.input_reads = 1000;
    stats.output_records = 1000;
    stats.missing_records = 0;
    stats.mapped = 950;
    stats.tentative = 30;
    stats.unmapped = 20;
    stats.position_evaluated = 950;
    stats.position_correct = 945;
    stats.lossless_consistent = 1000;
    stats.lossless_inconsistent = 0;
    stats.Compute();

    std::string summary = stats.Summary();
    EXPECT_NE(summary.find("Kill-Switch"), std::string::npos);
    EXPECT_NE(summary.find("Lossless Guarantee"), std::string::npos);
    EXPECT_NE(summary.find("PASS"), std::string::npos);
}

TEST(ValidationStats, SummaryShowsFailReason) {
    ValidationStats stats;
    stats.input_reads = 100;
    stats.output_records = 98;
    stats.missing_records = 2;
    stats.Compute();

    std::string summary = stats.Summary();
    EXPECT_NE(summary.find("FAIL"), std::string::npos);
    EXPECT_NE(summary.find("Lossless invariant violated"), std::string::npos);
}
