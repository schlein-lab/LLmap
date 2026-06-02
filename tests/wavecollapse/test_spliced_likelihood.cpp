// LLmap — Spliced likelihood + kernel tests.

#include "wavecollapse/spliced_likelihood.h"

#include "anchor/anchor_record.h"
#include "annot/junction_db.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace llmap::wavecollapse;

TEST(SplicedLikelihood, FactoredEmptySegmentsReturnsFloor) {
    // No segments at all — lossless floor, never 0.
    auto v = SplicedLikelihoodFactored({}, {});
    EXPECT_GT(v, 0.0f);
    EXPECT_LE(v, 1.0e-5f);  // around kFloor
}

TEST(SplicedLikelihood, FactoredSingleSegmentNoJunctions) {
    // One segment 0.9, no junctions → 0.9
    auto v = SplicedLikelihoodFactored({0.9f}, {});
    EXPECT_NEAR(v, 0.9f, 1e-5f);
}

TEST(SplicedLikelihood, FactoredMultipleSegmentsMultiplyWithJunctions) {
    // 2 segs × 1 junction: 0.9 * 0.8 * 0.95
    auto v = SplicedLikelihoodFactored({0.9f, 0.8f}, {0.95f});
    EXPECT_NEAR(v, 0.684f, 1e-3f);
}

TEST(SplicedLikelihood, FactoredShapeMismatchReturnsFloor) {
    // 2 segments + 2 junctions (should be 1) → shape mismatch
    auto v = SplicedLikelihoodFactored({0.9f, 0.8f}, {0.5f, 0.5f});
    EXPECT_LE(v, 1.0e-5f);
}

TEST(SplicedLikelihood, FactoredOnesEqualsOne) {
    auto v = SplicedLikelihoodFactored({1.0f, 1.0f, 1.0f}, {1.0f, 1.0f});
    EXPECT_NEAR(v, 1.0f, 1e-5f);
}

// ===========================================================================
// JunctionProbability
// ===========================================================================

TEST(JunctionProbability, GencodeAnnotatedBoostsAboveBaseline) {
    llmap::anchor::ExonBoundary b;
    b.donor_score = 0.9f;
    b.acceptor_score = 0.9f;

    llmap::annot::JunctionEvidence ev_none{};
    llmap::annot::JunctionEvidence ev_gen{};
    ev_gen.in_gencode = true;

    auto p_none = JunctionProbability(b, ev_none);
    auto p_gen  = JunctionProbability(b, ev_gen);
    EXPECT_GT(p_gen, p_none);
}

TEST(JunctionProbability, MissingPwmScoresGetNeutralFallback) {
    llmap::anchor::ExonBoundary b;  // both scores default 0
    llmap::annot::JunctionEvidence ev{};
    ev.in_gencode = true;
    auto p = JunctionProbability(b, ev);
    // 0.5 * 0.5 * 1.05 = 0.2625 ; tissue_compat=1 default
    EXPECT_NEAR(p, 0.2625f, 0.01f);
}

TEST(JunctionProbability, BackSpliceCapsAt095) {
    llmap::anchor::ExonBoundary b;
    b.donor_score = 1.0f;
    b.acceptor_score = 1.0f;
    b.spliceosome_class = 3;  // back-splice
    llmap::annot::JunctionEvidence ev{};
    ev.in_gencode = true;
    auto p = JunctionProbability(b, ev);
    EXPECT_LE(p, 0.95f);
}

TEST(JunctionProbability, TissueCompatibilityScalesDown) {
    llmap::anchor::ExonBoundary b;
    b.donor_score = 0.9f;
    b.acceptor_score = 0.9f;
    llmap::annot::JunctionEvidence ev{};
    ev.in_gencode = true;
    auto p_full = JunctionProbability(b, ev, /*tissue_compat=*/1.0f);
    auto p_half = JunctionProbability(b, ev, /*tissue_compat=*/0.5f);
    EXPECT_GT(p_full, p_half);
}

// ===========================================================================
// KSpliced kernel
// ===========================================================================

TEST(KSpliced, DistanceZeroYieldsAtMostUnity) {
    auto k = KSpliced(/*delta=*/0, /*sigma=*/10000.0f,
                       /*junction=*/1.0f, /*tissue=*/1.0f);
    EXPECT_LE(k, 1.0f);
    EXPECT_GT(k, 0.99f);
}

TEST(KSpliced, FarDistanceDecaysToFloor) {
    // 10 σ out — Gaussian ≈ exp(-50) ≈ 0; floor takes over.
    auto k = KSpliced(/*delta=*/100000, /*sigma=*/10000.0f,
                       1.0f, 1.0f);
    EXPECT_GT(k, 0.0f);   // never 0 (lossless floor)
    EXPECT_LE(k, 1.0e-5f);
}

TEST(KSpliced, JunctionFactorMultiplies) {
    auto k_high = KSpliced(0, 10000.0f, 1.0f, 1.0f);
    auto k_low  = KSpliced(0, 10000.0f, 0.1f, 1.0f);
    EXPECT_GT(k_high, k_low);
}

TEST(KSpliced, BadSigmaReturnsFloor) {
    auto k = KSpliced(0, /*sigma=*/0.0f, 1.0f, 1.0f);
    EXPECT_LE(k, 1.0e-5f);
}

// ===========================================================================
// TissueCompatibility
// ===========================================================================

TEST(TissueCompatibility, NegativeTpmIsNeutral) {
    EXPECT_NEAR(TissueCompatibility(-1.0f), 1.0f, 1e-5f);
}

TEST(TissueCompatibility, HighTpmApproachesUnity) {
    EXPECT_GT(TissueCompatibility(1000.0f), 0.85f);
}

TEST(TissueCompatibility, ZeroTpmFloorsAt020) {
    EXPECT_GE(TissueCompatibility(0.0f), 0.20f);
    EXPECT_LE(TissueCompatibility(0.0f), 0.30f);
}
