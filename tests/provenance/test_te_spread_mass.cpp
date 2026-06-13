// Unit tests for TE-scale WaveCollapse spread-mass.

#include "provenance/te_spread_mass.h"

#include <gtest/gtest.h>

namespace llmap::provenance {
namespace {

TePlacement P(std::int32_t score, bool anchor) {
    TePlacement p;
    p.chrom = "chr1";
    p.align_score = score;
    p.has_unique_anchor = anchor;
    p.te_class = TeClass::Alu;
    return p;
}

TEST(TeSpreadMass, EmptyIsAmbiguousZeroConfidence) {
    const auto r = ResolveTeSpreadMass({});
    EXPECT_TRUE(r.is_ambiguous);
    EXPECT_FLOAT_EQ(r.confidence, 0.0f);
    EXPECT_EQ(r.collapsed, -1);
}

TEST(TeSpreadMass, UniqueAnchorCollapsesConfidently) {
    // Two equal TE copies, but one has a unique flank → collapse there.
    std::vector<TePlacement> c = {P(100, false), P(100, true)};
    const auto r = ResolveTeSpreadMass(c);
    EXPECT_EQ(r.collapsed, 1);
    EXPECT_FALSE(r.is_ambiguous);
    EXPECT_FLOAT_EQ(r.confidence, 1.0f);
    EXPECT_EQ(r.mapq, 60);
    EXPECT_FLOAT_EQ(r.posterior[1], 1.0f);
}

TEST(TeSpreadMass, EqualCopiesKeepSpreadHonestLowMapq) {
    // Two identical Alu copies, no unique anchor → 50/50, honest low MAPQ
    // (NOT a faked 60). This is the whole point vs minimap2.
    std::vector<TePlacement> c = {P(100, false), P(100, false)};
    const auto r = ResolveTeSpreadMass(c);
    EXPECT_TRUE(r.is_ambiguous);
    EXPECT_EQ(r.collapsed, -1);
    EXPECT_NEAR(r.posterior[0], 0.5f, 1e-4);
    EXPECT_NEAR(r.posterior[1], 0.5f, 1e-4);
    EXPECT_NEAR(r.confidence, 0.0f, 1e-4);   // max entropy → zero confidence
    EXPECT_LT(r.mapq, 10);                    // honest: ambiguous, low MAPQ
}

TEST(TeSpreadMass, DominantCopyHigherConfidence) {
    // One copy fits much better → mass concentrates, confidence + MAPQ rise.
    std::vector<TePlacement> c = {P(200, false), P(100, false), P(100, false)};
    const auto r = ResolveTeSpreadMass(c, 30.0f);
    EXPECT_GT(r.posterior[0], r.posterior[1]);
    EXPECT_GT(r.confidence, 0.0f);
    EXPECT_GT(r.mapq, ResolveTeSpreadMass({P(100,false),P(100,false),P(100,false)}, 30.0f).mapq);
}

TEST(TeSpreadMass, PopulationPriorShiftsMass) {
    // Two equal-alignment copies, but the pangenome supports copy 0 far more
    // (it is the consistently-mappable placement across populations) → mass moves
    // to copy 0. The Operator's directive: population reproducibility informs the
    // placement likelihood. Absent support (<0) → even spread (the ambiguity case).
    std::vector<TePlacement> c = {P(100, false), P(100, false)};
    c[0].population_support = 0.95f;   // sharp across the pangenome
    c[1].population_support = 0.05f;   // rarely the right placement
    const auto r = ResolveTeSpreadMass(c);
    EXPECT_GT(r.posterior[0], r.posterior[1]);
    EXPECT_GT(r.posterior[0], 0.7f);   // population prior concentrates the mass

    // Same scores, no population data → honest 50/50 (consistent blurriness).
    std::vector<TePlacement> u = {P(100, false), P(100, false)};
    const auto ru = ResolveTeSpreadMass(u);
    EXPECT_NEAR(ru.posterior[0], 0.5f, 1e-4);
}

TEST(TeSpreadMass, ManyEqualCopiesNearZeroMapq) {
    // The Alu reality: a read fitting ~hundreds of near-identical copies.
    std::vector<TePlacement> c(50, P(100, false));
    const auto r = ResolveTeSpreadMass(c);
    EXPECT_TRUE(r.is_ambiguous);
    EXPECT_LE(r.mapq, 1);                      // ~0 — honest "could be anywhere"
    EXPECT_NEAR(r.confidence, 0.0f, 1e-3);
}

}  // namespace
}  // namespace llmap::provenance
