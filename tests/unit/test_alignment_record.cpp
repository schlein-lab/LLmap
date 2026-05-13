// LLmap Phase 0 — lossless invariant tests.
//
// These tests encode the non-negotiable promise: every input read produces a
// record with status ∈ {Mapped, Tentative, Unmapped}, no silent drops.

#include "core/alignment_record.h"

#include <gtest/gtest.h>

using namespace llmap;

namespace {

AlignmentHit dummy_hit(std::string target = "chr14") {
    AlignmentHit h;
    h.target_id = std::move(target);
    h.start = 105'000'000;
    h.end = 105'010'000;
    h.cigar.ops = "10000M";
    h.score = 9500;
    h.nm = 50;
    return h;
}

TentativeTarget dummy_target(std::string target = "chr14") {
    TentativeTarget t;
    t.target_id = std::move(target);
    t.approx_start = 105'000'000;
    t.approx_end = 105'010'000;
    t.n_seeds = 12;
    t.partial_chain_score = 4200;
    t.sequence_identity_estimate = 0.92f;
    t.final_probability = 0.55f;
    return t;
}

}  // namespace

TEST(AlignmentRecord, MappedFactoryProducesConsistentRecord) {
    auto r = make_mapped("read_001", 10'000, dummy_hit());
    EXPECT_EQ(r.read_id, "read_001");
    EXPECT_EQ(r.read_len, 10'000u);
    EXPECT_EQ(r.status, AlignmentStatus::Mapped);
    EXPECT_TRUE(r.primary.has_value());
    EXPECT_FALSE(r.rejection_reason.has_value());
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, TentativeFactoryProducesConsistentRecord) {
    auto r = make_tentative(
        "read_002", 12'500,
        {dummy_target("chr14"), dummy_target("chr15")},
        RejectionReason::DidNotConverge
    );
    EXPECT_EQ(r.status, AlignmentStatus::Tentative);
    EXPECT_FALSE(r.primary.has_value());
    EXPECT_EQ(r.tentative_targets.size(), 2u);
    ASSERT_TRUE(r.rejection_reason.has_value());
    EXPECT_EQ(*r.rejection_reason, RejectionReason::DidNotConverge);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, UnmappedFactoryProducesConsistentRecord) {
    auto r = make_unmapped("read_003", 8'000, RejectionReason::NoSeeds);
    EXPECT_EQ(r.status, AlignmentStatus::Unmapped);
    EXPECT_FALSE(r.primary.has_value());
    EXPECT_TRUE(r.tentative_targets.empty());
    ASSERT_TRUE(r.rejection_reason.has_value());
    EXPECT_EQ(*r.rejection_reason, RejectionReason::NoSeeds);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, EmptyReadIdViolatesInvariant) {
    auto r = make_unmapped("", 1'000, RejectionReason::NoSeeds);
    EXPECT_FALSE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, MappedWithoutPrimaryViolatesInvariant) {
    AlignmentRecord r;
    r.read_id = "read_x";
    r.status = AlignmentStatus::Mapped;
    // primary intentionally unset
    EXPECT_FALSE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, TentativeWithoutTargetsViolatesInvariant) {
    AlignmentRecord r;
    r.read_id = "read_x";
    r.status = AlignmentStatus::Tentative;
    r.rejection_reason = RejectionReason::DidNotConverge;
    // tentative_targets intentionally empty
    EXPECT_FALSE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, UnmappedWithoutRejectionReasonViolatesInvariant) {
    AlignmentRecord r;
    r.read_id = "read_x";
    r.status = AlignmentStatus::Unmapped;
    // rejection_reason intentionally unset
    EXPECT_FALSE(r.is_lossless_consistent());
}

TEST(AlignmentRecord, LosslessInvariantOverBulkInputs) {
    // Encodes the V1.0 promise: 1000 input reads → 1000 consistent records.
    constexpr std::size_t N = 1000;
    std::vector<AlignmentRecord> records;
    records.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto id = "read_" + std::to_string(i);
        if (i % 3 == 0) {
            records.push_back(make_mapped(id, 10'000, dummy_hit()));
        } else if (i % 3 == 1) {
            records.push_back(make_tentative(
                id, 10'000, {dummy_target()}, RejectionReason::DidNotConverge));
        } else {
            records.push_back(make_unmapped(id, 10'000, RejectionReason::NoSeeds));
        }
    }
    EXPECT_EQ(records.size(), N);
    for (const auto& r : records) {
        EXPECT_TRUE(r.is_lossless_consistent()) << "record " << r.read_id << " inconsistent";
    }
}
