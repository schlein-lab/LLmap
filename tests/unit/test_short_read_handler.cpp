// LLmap — short-read pipeline (paired-end + equivalence-class) tests.

#include "io/short_read_handler.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace llmap::io;

TEST(ShortReadHandler, BuildEquivalenceClassNormalisesWeights) {
    ShortReadConfig cfg;
    cfg.min_member_weight = 0.05f;
    std::vector<std::pair<std::string, float>> hits = {
        {"ENST.A", 60.0f},
        {"ENST.B", 30.0f},
        {"ENST.C", 10.0f},
    };
    auto mem = BuildEquivalenceClass(hits, cfg);
    ASSERT_EQ(mem.size(), 3u);

    float sum = 0.0f;
    for (const auto& m : mem) sum += m.weight;
    EXPECT_NEAR(sum, 1.0f, 1e-5f);

    // Sorted descending.
    EXPECT_EQ(mem[0].transcript_id, "ENST.A");
    EXPECT_NEAR(mem[0].weight, 0.6f, 1e-5f);
    EXPECT_EQ(mem[1].transcript_id, "ENST.B");
    EXPECT_NEAR(mem[1].weight, 0.3f, 1e-5f);
    EXPECT_EQ(mem[2].transcript_id, "ENST.C");
    EXPECT_NEAR(mem[2].weight, 0.1f, 1e-5f);
}

TEST(ShortReadHandler, EquivalenceClassDropsBelowMinWeight) {
    ShortReadConfig cfg;
    cfg.min_member_weight = 0.10f;   // 10 % cut
    std::vector<std::pair<std::string, float>> hits = {
        {"A", 90.0f},
        {"B", 5.0f},     // weight 0.05 < 0.10 → dropped
        {"C", 5.0f},     // weight 0.05 < 0.10 → dropped
    };
    auto mem = BuildEquivalenceClass(hits, cfg);
    ASSERT_EQ(mem.size(), 1u);
    EXPECT_EQ(mem[0].transcript_id, "A");
    EXPECT_NEAR(mem[0].weight, 1.0f, 1e-5f);
}

TEST(ShortReadHandler, EmptyHitsYieldsEmptyClass) {
    ShortReadConfig cfg;
    EXPECT_TRUE(BuildEquivalenceClass({}, cfg).empty());
    std::vector<std::pair<std::string, float>> zero = {{"A", 0.0f}};
    EXPECT_TRUE(BuildEquivalenceClass(zero, cfg).empty());
}

TEST(ShortReadHandler, ImpliedInsertSizeAbsDifference) {
    EXPECT_EQ(ImpliedInsertSize(100, 400), 300u);
    EXPECT_EQ(ImpliedInsertSize(400, 100), 300u);
    EXPECT_EQ(ImpliedInsertSize(100, 100),   0u);
    EXPECT_EQ(ImpliedInsertSize(0,   400),   0u);
}

TEST(ShortReadHandler, BuildOutcomeMergesMates) {
    ShortReadConfig cfg;
    cfg.emit_equivalence_classes = true;
    cfg.min_member_weight = 0.05f;
    FastqRecordPair rec{.read_id = "r1"};
    auto outcome = BuildOutcome(
        rec,
        /*hits_r1=*/ {{"A", 100.0f}, {"B", 10.0f}},
        /*hits_r2=*/ {{"A",  50.0f}, {"C", 10.0f}},
        /*implied_insert=*/ std::nullopt,
        /*junction_seen=*/ false,
        cfg);
    ASSERT_EQ(outcome.read_id, "r1");
    // A:150, B:10, C:10 — normalised to 0.882 / 0.059 / 0.059
    // After 0.05 cutoff all three survive (weights >0.05).
    ASSERT_EQ(outcome.members.size(), 3u);
    EXPECT_EQ(outcome.members[0].transcript_id, "A");
    EXPECT_GT(outcome.members[0].weight, 0.85f);
}

TEST(ShortReadHandler, BuildOutcomeNoEquivalenceClassPicksBest) {
    ShortReadConfig cfg;
    cfg.emit_equivalence_classes = false;
    FastqRecordPair rec{.read_id = "r1"};
    auto outcome = BuildOutcome(
        rec,
        {{"A", 100.0f}, {"B", 50.0f}},
        {},
        std::nullopt,
        false,
        cfg);
    ASSERT_EQ(outcome.members.size(), 1u);
    EXPECT_EQ(outcome.members[0].transcript_id, "A");
    EXPECT_NEAR(outcome.members[0].weight, 1.0f, 1e-5f);
}

TEST(ShortReadHandler, BuildOutcomeFlagsJunctionOnLongInsert) {
    ShortReadConfig cfg;
    cfg.expected_insert_size = 300;
    FastqRecordPair rec{.read_id = "r1"};
    // implied 5000 bp insert > 1.5 × 300 = 450 → junction-implied.
    auto outcome = BuildOutcome(rec, {{"A", 1.0f}}, {},
                                  /*implied_insert=*/ 5000,
                                  false, cfg);
    EXPECT_TRUE(outcome.junction_inferred);
    EXPECT_EQ(outcome.implied_insert_size, 5000u);
}

TEST(ShortReadHandler, BuildOutcomeNoEvidenceYieldsEmpty) {
    ShortReadConfig cfg;
    FastqRecordPair rec{.read_id = "r1"};
    auto outcome = BuildOutcome(rec, {}, {}, std::nullopt, false, cfg);
    EXPECT_TRUE(outcome.members.empty());
    EXPECT_FALSE(outcome.junction_inferred);
}

TEST(ShortReadHandler, ExplicitJunctionFlagPreserved) {
    ShortReadConfig cfg;
    FastqRecordPair rec{.read_id = "r1"};
    auto outcome = BuildOutcome(rec, {{"A", 1.0f}}, {},
                                  /*implied_insert=*/ 100,  // small
                                  /*junction_seen=*/ true, cfg);
    EXPECT_TRUE(outcome.junction_inferred);
}
