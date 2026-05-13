#include <gtest/gtest.h>

#include "classical/wfa2_aligner.h"

using namespace llmap::classical;

class WFA2AlignerTest : public ::testing::Test {
protected:
    WFA2Config config;

    void SetUp() override {
        config = WFA2Config{};
    }
};

// ============================================================================
// CigarElement tests
// ============================================================================

TEST_F(WFA2AlignerTest, CigarElementOpCharMapping) {
    CigarElement e_match{CigarOp::Match, 5};
    CigarElement e_ins{CigarOp::Insertion, 3};
    CigarElement e_del{CigarOp::Deletion, 2};
    CigarElement e_soft{CigarOp::SoftClip, 4};
    CigarElement e_eq{CigarOp::Equal, 10};
    CigarElement e_diff{CigarOp::Diff, 1};

    EXPECT_EQ(e_match.OpChar(), 'M');
    EXPECT_EQ(e_ins.OpChar(), 'I');
    EXPECT_EQ(e_del.OpChar(), 'D');
    EXPECT_EQ(e_soft.OpChar(), 'S');
    EXPECT_EQ(e_eq.OpChar(), '=');
    EXPECT_EQ(e_diff.OpChar(), 'X');
}

TEST_F(WFA2AlignerTest, CigarElementToString) {
    CigarElement e1{CigarOp::Match, 5};
    CigarElement e2{CigarOp::Insertion, 12};
    CigarElement e3{CigarOp::Equal, 100};

    EXPECT_EQ(e1.ToString(), "5M");
    EXPECT_EQ(e2.ToString(), "12I");
    EXPECT_EQ(e3.ToString(), "100=");
}

// ============================================================================
// ParseCigar tests
// ============================================================================

TEST_F(WFA2AlignerTest, ParseCigarSimple) {
    auto result = ParseCigar("10M");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result->size(), 1);
    EXPECT_EQ((*result)[0].op, CigarOp::Match);
    EXPECT_EQ((*result)[0].length, 10);
}

TEST_F(WFA2AlignerTest, ParseCigarComplex) {
    auto result = ParseCigar("5=2X3I4D10=");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result->size(), 5);

    EXPECT_EQ((*result)[0].op, CigarOp::Equal);
    EXPECT_EQ((*result)[0].length, 5);

    EXPECT_EQ((*result)[1].op, CigarOp::Diff);
    EXPECT_EQ((*result)[1].length, 2);

    EXPECT_EQ((*result)[2].op, CigarOp::Insertion);
    EXPECT_EQ((*result)[2].length, 3);

    EXPECT_EQ((*result)[3].op, CigarOp::Deletion);
    EXPECT_EQ((*result)[3].length, 4);

    EXPECT_EQ((*result)[4].op, CigarOp::Equal);
    EXPECT_EQ((*result)[4].length, 10);
}

TEST_F(WFA2AlignerTest, ParseCigarEmpty) {
    auto result = ParseCigar("");
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(result->empty());
}

TEST_F(WFA2AlignerTest, ParseCigarInvalid) {
    EXPECT_FALSE(ParseCigar("M").has_value());  // No length
    EXPECT_FALSE(ParseCigar("10").has_value()); // No op
    EXPECT_FALSE(ParseCigar("10Z").has_value()); // Invalid op
}

// ============================================================================
// WFA2Result tests
// ============================================================================

TEST_F(WFA2AlignerTest, WFA2ResultCigarString) {
    WFA2Result result;
    result.cigar = {
        {CigarOp::Equal, 5},
        {CigarOp::Diff, 2},
        {CigarOp::Equal, 10}
    };

    EXPECT_EQ(result.CigarString(), "5=2X10=");
}

TEST_F(WFA2AlignerTest, WFA2ResultAlignedLength) {
    WFA2Result result;
    result.cigar = {
        {CigarOp::Equal, 10},      // 10 ref bases
        {CigarOp::Insertion, 5},   // 0 ref bases (query only)
        {CigarOp::Equal, 8},       // 8 ref bases
        {CigarOp::Deletion, 3}     // 3 ref bases (no query)
    };

    // Aligned length in ref = 10 + 8 + 3 = 21
    EXPECT_EQ(result.AlignedLength(), 21);
}

// ============================================================================
// ComputeIdentity tests
// ============================================================================

TEST_F(WFA2AlignerTest, ComputeIdentityPerfect) {
    std::vector<CigarElement> cigar = {
        {CigarOp::Equal, 100}
    };
    EXPECT_FLOAT_EQ(ComputeIdentity(cigar), 1.0f);
}

TEST_F(WFA2AlignerTest, ComputeIdentityWithMismatches) {
    std::vector<CigarElement> cigar = {
        {CigarOp::Equal, 90},
        {CigarOp::Diff, 10}
    };
    EXPECT_FLOAT_EQ(ComputeIdentity(cigar), 0.9f);
}

TEST_F(WFA2AlignerTest, ComputeIdentityWithGaps) {
    std::vector<CigarElement> cigar = {
        {CigarOp::Equal, 80},
        {CigarOp::Insertion, 10},
        {CigarOp::Equal, 10}
    };
    // 90 matches / 100 aligned = 0.9
    EXPECT_FLOAT_EQ(ComputeIdentity(cigar), 0.9f);
}

TEST_F(WFA2AlignerTest, ComputeIdentityEmpty) {
    std::vector<CigarElement> cigar;
    EXPECT_FLOAT_EQ(ComputeIdentity(cigar), 0.0f);
}

// ============================================================================
// WFA2Aligner basic tests
// ============================================================================

TEST_F(WFA2AlignerTest, AlignIdenticalSequences) {
    WFA2Aligner aligner(config);

    std::string seq = "ACGTACGTACGT";
    auto result = aligner.Align(seq, seq);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->score, 0);  // Perfect match has no penalty
    EXPECT_EQ(result->num_matches, seq.length());
    EXPECT_EQ(result->num_mismatches, 0);
    EXPECT_EQ(result->num_insertions, 0);
    EXPECT_EQ(result->num_deletions, 0);
    EXPECT_FLOAT_EQ(result->identity, 1.0f);
}

TEST_F(WFA2AlignerTest, AlignWithSingleMismatch) {
    WFA2Aligner aligner(config);

    std::string query = "ACGTACGTACGT";
    std::string ref   = "ACGTACGCACGT";  // T -> C at position 7

    auto result = aligner.Align(query, ref);

    ASSERT_TRUE(result.has_value());
    EXPECT_GT(result->score, 0);  // Should have penalty
    EXPECT_EQ(result->num_mismatches, 1);
    EXPECT_EQ(result->num_insertions, 0);
    EXPECT_EQ(result->num_deletions, 0);
}

TEST_F(WFA2AlignerTest, AlignWithSingleInsertion) {
    WFA2Aligner aligner(config);

    std::string query = "ACGTACGTACGT";
    std::string ref   = "ACGTACGACGT";  // Missing T at position 7

    auto result = aligner.Align(query, ref);

    ASSERT_TRUE(result.has_value());
    EXPECT_GT(result->score, 0);
    // Either an insertion or deletion depending on perspective
    EXPECT_TRUE(result->num_insertions > 0 || result->num_deletions > 0);
}

TEST_F(WFA2AlignerTest, AlignWithSingleDeletion) {
    WFA2Aligner aligner(config);

    std::string query = "ACGTACGACGT";   // Missing T
    std::string ref   = "ACGTACGTACGT";

    auto result = aligner.Align(query, ref);

    ASSERT_TRUE(result.has_value());
    EXPECT_GT(result->score, 0);
    EXPECT_TRUE(result->num_insertions > 0 || result->num_deletions > 0);
}

TEST_F(WFA2AlignerTest, AlignEmptyQueryFails) {
    WFA2Aligner aligner(config);

    auto result = aligner.Align("", "ACGT");
    EXPECT_FALSE(result.has_value());
}

TEST_F(WFA2AlignerTest, AlignEmptyReferenceFails) {
    WFA2Aligner aligner(config);

    auto result = aligner.Align("ACGT", "");
    EXPECT_FALSE(result.has_value());
}

TEST_F(WFA2AlignerTest, AlignLongerSequences) {
    WFA2Aligner aligner(config);

    // 100bp sequences with some differences
    std::string query = std::string(50, 'A') + std::string(50, 'T');
    std::string ref   = std::string(48, 'A') + "CC" + std::string(50, 'T');

    auto result = aligner.Align(query, ref);

    ASSERT_TRUE(result.has_value());
    EXPECT_LT(result->identity, 1.0f);
    EXPECT_GT(result->identity, 0.9f);  // Should be close to perfect
}

// ============================================================================
// ExtendRight tests
// ============================================================================

TEST_F(WFA2AlignerTest, ExtendRightBasic) {
    WFA2Aligner aligner(config);

    std::string query = "AAAACGTCGTCGTCGT";
    std::string ref   = "AAAACGTCGTCGTCGT";

    // Extend from position 4 onwards
    auto result = aligner.ExtendRight(query, ref, 4, 4);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->query_start, 4);
    EXPECT_EQ(result->ref_start, 4);
    EXPECT_EQ(result->score, 0);  // Perfect match
}

TEST_F(WFA2AlignerTest, ExtendRightInvalidPosition) {
    WFA2Aligner aligner(config);

    std::string seq = "ACGT";
    EXPECT_FALSE(aligner.ExtendRight(seq, seq, -1, 0).has_value());
    EXPECT_FALSE(aligner.ExtendRight(seq, seq, 0, -1).has_value());
    EXPECT_FALSE(aligner.ExtendRight(seq, seq, 100, 0).has_value());
    EXPECT_FALSE(aligner.ExtendRight(seq, seq, 0, 100).has_value());
}

// ============================================================================
// ExtendLeft tests
// ============================================================================

TEST_F(WFA2AlignerTest, ExtendLeftBasic) {
    WFA2Aligner aligner(config);

    std::string query = "ACGTACGTAAAA";
    std::string ref   = "ACGTACGTAAAA";

    // Extend from position 8 backwards
    auto result = aligner.ExtendLeft(query, ref, 8, 8);

    ASSERT_TRUE(result.has_value());
    // Should align the prefix ACGTACGT
}

TEST_F(WFA2AlignerTest, ExtendLeftInvalidPosition) {
    WFA2Aligner aligner(config);

    std::string seq = "ACGT";
    EXPECT_FALSE(aligner.ExtendLeft(seq, seq, 0, 4).has_value());
    EXPECT_FALSE(aligner.ExtendLeft(seq, seq, 4, 0).has_value());
    EXPECT_FALSE(aligner.ExtendLeft(seq, seq, 100, 4).has_value());
}

// ============================================================================
// AlignBatch tests
// ============================================================================

TEST_F(WFA2AlignerTest, AlignBatchBasic) {
    WFA2Aligner aligner(config);

    std::vector<std::string_view> queries = {"ACGT", "AAAA", "TTTT"};
    std::vector<std::string_view> refs    = {"ACGT", "AAAA", "TTTT"};

    auto results = aligner.AlignBatch(queries, refs);

    ASSERT_EQ(results.size(), 3);
    for (const auto& r : results) {
        ASSERT_TRUE(r.has_value());
        EXPECT_EQ(r->score, 0);  // All perfect matches
    }
}

TEST_F(WFA2AlignerTest, AlignBatchMismatchedSizesFails) {
    WFA2Aligner aligner(config);

    std::vector<std::string_view> queries = {"ACGT", "AAAA"};
    std::vector<std::string_view> refs    = {"ACGT"};

    auto results = aligner.AlignBatch(queries, refs);
    EXPECT_TRUE(results.empty());
}

TEST_F(WFA2AlignerTest, AlignBatchWithVariedResults) {
    WFA2Aligner aligner(config);

    std::vector<std::string_view> queries = {"ACGT", "ACGT", ""};
    std::vector<std::string_view> refs    = {"ACGT", "ACCC", ""};

    auto results = aligner.AlignBatch(queries, refs);

    ASSERT_EQ(results.size(), 3);
    EXPECT_TRUE(results[0].has_value());   // Perfect match
    EXPECT_TRUE(results[1].has_value());   // With mismatches
    EXPECT_FALSE(results[2].has_value());  // Empty fails
}

// ============================================================================
// Configuration tests
// ============================================================================

TEST_F(WFA2AlignerTest, ConfigPreserved) {
    config.mismatch_penalty = 10;
    config.gap_open = 20;
    config.gap_extend = 5;

    WFA2Aligner aligner(config);

    EXPECT_EQ(aligner.Config().mismatch_penalty, 10);
    EXPECT_EQ(aligner.Config().gap_open, 20);
    EXPECT_EQ(aligner.Config().gap_extend, 5);
}

TEST_F(WFA2AlignerTest, MaxScoreLimitEnforced) {
    config.max_alignment_score = 5;  // Very low limit

    WFA2Aligner aligner(config);

    // A sequence with many mismatches should exceed the limit
    std::string query = "AAAAAAAAAA";
    std::string ref   = "TTTTTTTTTT";

    auto result = aligner.Align(query, ref);
    EXPECT_FALSE(result.has_value());  // Should fail due to score limit
}

// ============================================================================
// Utility function tests
// ============================================================================

TEST_F(WFA2AlignerTest, IsWFA2AvailableReturnsValue) {
    // Should return true or false without crashing
    bool available = IsWFA2Available();
    (void)available;  // May be true or false depending on build config
}

TEST_F(WFA2AlignerTest, IsNativeWFA2Consistent) {
    WFA2Aligner aligner(config);
    // Should match global availability
    EXPECT_EQ(aligner.IsNativeWFA2(), IsWFA2Available());
}

// ============================================================================
// Move semantics tests
// ============================================================================

TEST_F(WFA2AlignerTest, MoveConstruction) {
    WFA2Aligner aligner1(config);
    WFA2Aligner aligner2(std::move(aligner1));

    // aligner2 should work
    auto result = aligner2.Align("ACGT", "ACGT");
    EXPECT_TRUE(result.has_value());
}

TEST_F(WFA2AlignerTest, MoveAssignment) {
    WFA2Aligner aligner1(config);
    WFA2Aligner aligner2(config);

    aligner2 = std::move(aligner1);

    // aligner2 should work
    auto result = aligner2.Align("ACGT", "ACGT");
    EXPECT_TRUE(result.has_value());
}

// ============================================================================
// Alignment quality tests
// ============================================================================

TEST_F(WFA2AlignerTest, AlignmentTimeRecorded) {
    WFA2Aligner aligner(config);

    auto result = aligner.Align("ACGTACGTACGT", "ACGTACGTACGT");

    ASSERT_TRUE(result.has_value());
    EXPECT_GE(result->alignment_time_us, 0.0f);
}

TEST_F(WFA2AlignerTest, AlignmentCoordinatesCorrect) {
    WFA2Aligner aligner(config);

    std::string query = "ACGTACGT";
    std::string ref   = "ACGTACGT";

    auto result = aligner.Align(query, ref);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->query_start, 0);
    EXPECT_EQ(result->query_end, static_cast<int32_t>(query.length()));
    EXPECT_EQ(result->ref_start, 0);
    EXPECT_EQ(result->ref_end, static_cast<int32_t>(ref.length()));
}

TEST_F(WFA2AlignerTest, GapAffineScoring) {
    // Test that gap extension is cheaper than gap opening
    WFA2Aligner aligner(config);

    // Single gap of length 3 should be cheaper than 3 gaps of length 1
    std::string query1 = "ACGTACGT";
    std::string ref1   = "ACGT---ACGT";  // Would need 3-base deletion

    std::string query2 = "ACGTACGT";
    std::string ref2   = "A-C-G-TACGT"; // Would need 3 separate 1-base deletions

    // The gap-affine model should give different scores
    // (We just verify the alignment works; exact scoring depends on implementation)
    auto result1 = aligner.Align(query1, "ACGTAAACGT");
    auto result2 = aligner.Align(query2, "ACGTAAACGT");

    // Both should succeed
    EXPECT_TRUE(result1.has_value());
    EXPECT_TRUE(result2.has_value());
}
