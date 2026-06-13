// LLmap — splice_determinism unit tests.
//
// Verifies the operator's D(pos) examples: constitutive exon → 100, a cassette
// position excluded in 1/10 reads → 90, a 5/5 split → 50, plus junction usage
// counting and the exonic-only exclusion restriction.

#include "mapping/splice_determinism.h"

#include <gtest/gtest.h>

#include <string>

namespace {

using llmap::mapping::Determinism;
using llmap::mapping::DeterminismAccumulator;
using llmap::mapping::PositionCounts;
using llmap::mapping::SortedJunctions;

// ---------------------------------------------------------------------------
// Determinism() formula
// ---------------------------------------------------------------------------

TEST(SpliceDeterminism, FormulaModalFraction) {
    EXPECT_DOUBLE_EQ(Determinism({10, 0}), 100.0);  // constitutive
    EXPECT_DOUBLE_EQ(Determinism({9, 1}), 90.0);    // 9/10 modal (included)
    EXPECT_DOUBLE_EQ(Determinism({1, 9}), 90.0);    // 9/10 modal (excluded)
    EXPECT_DOUBLE_EQ(Determinism({5, 5}), 50.0);    // maximal ambiguity
    EXPECT_DOUBLE_EQ(Determinism({0, 0}), 0.0);     // uncovered
}

// ---------------------------------------------------------------------------
// Accumulation: cassette exon at position 100
// ---------------------------------------------------------------------------

TEST(SpliceDeterminism, CassettePositionScores90) {
    DeterminismAccumulator acc;
    // 9 reads include the exon body (cover [95,105) → includes pos 95..104).
    for (int i = 0; i < 9; ++i) acc.MarkIncluded("chrS", 95, "10M");
    // 1 read skips [100,110) as an intron: [95,100)M [100,110)N [110,115)M.
    acc.MarkIncluded("chrS", 95, "5M10N5M");
    // Pass 2: exclusions (no-op for the pure-M includers; the N read excludes).
    for (int i = 0; i < 9; ++i) acc.MarkExcluded("chrS", 95, "10M");
    acc.MarkExcluded("chrS", 95, "5M10N5M");

    const auto& pm = acc.positions.at("chrS");

    // pos 95: included by all 10, never in an intron → constitutive → 100.
    EXPECT_EQ(pm.at(95).included, 10u);
    EXPECT_EQ(pm.at(95).excluded, 0u);
    EXPECT_DOUBLE_EQ(Determinism(pm.at(95)), 100.0);

    // pos 100: included by the 9, excluded by the 1 → cassette → 90.
    EXPECT_EQ(pm.at(100).included, 9u);
    EXPECT_EQ(pm.at(100).excluded, 1u);
    EXPECT_DOUBLE_EQ(Determinism(pm.at(100)), 90.0);
}

// ---------------------------------------------------------------------------
// Junction usage + ordering
// ---------------------------------------------------------------------------

TEST(SpliceDeterminism, JunctionUsageCounted) {
    DeterminismAccumulator acc;
    // 3 reads splice the same intron [100,110); 1 read a different one [100,120).
    for (int i = 0; i < 3; ++i) acc.MarkIncluded("chrS", 95, "5M10N5M");
    acc.MarkIncluded("chrS", 95, "5M20N5M");

    const auto j = SortedJunctions(acc, "chrS");
    ASSERT_EQ(j.size(), 2u);
    // Sorted by (donor, acceptor): both donor 100, acceptors 110 then 120.
    EXPECT_EQ(j[0].donor, 100u);
    EXPECT_EQ(j[0].acceptor, 110u);
    EXPECT_EQ(j[0].n_reads, 3u);
    EXPECT_EQ(j[1].acceptor, 120u);
    EXPECT_EQ(j[1].n_reads, 1u);
}

// ---------------------------------------------------------------------------
// Exclusion only touches positions that are exonic in some read
// ---------------------------------------------------------------------------

TEST(SpliceDeterminism, ExclusionRestrictedToExonicPositions) {
    DeterminismAccumulator acc;
    // A single read that only skips [100,110) and never includes it: the intron
    // positions must NOT appear (no read makes them exonic → not interesting).
    acc.MarkIncluded("chrS", 95, "5M10N5M");
    acc.MarkExcluded("chrS", 95, "5M10N5M");
    const auto& pm = acc.positions.at("chrS");
    EXPECT_EQ(pm.count(100), 0u);  // intron-only position absent
    EXPECT_EQ(pm.count(95), 1u);   // exonic position present
    EXPECT_EQ(pm.count(110), 1u);  // exonic (downstream exon) present
}

}  // namespace
