// LLmap — Multi-Signal Fusion Engine tests.
//
// Verify each of the 10 factors individually + the Product() floor +
// the disable-mask + a smoke "10-factor synthetic" integration.

#include "fusion/likelihood_factors.h"

#include "anchor/anchor_record.h"
#include "core/transcript_kind.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

using namespace llmap::fusion;
using llmap::anchor::AnchorRecord;
using llmap::anchor::AnchorSource;
using llmap::anchor::ExonBoundary;
using llmap::core::TranscriptKind;

namespace {

AnchorRecord MakeAnchor(TranscriptKind k = TranscriptKind::MatureMrna,
                         std::string seq = std::string(500, 'A'),
                         std::vector<std::string> tags = {}) {
    AnchorRecord r;
    r.anchor_id = "GENCODE:ENST.1:exon1";
    r.source = AnchorSource::Gencode;
    r.kind = k;
    r.sequence = std::move(seq);
    r.tags = std::move(tags);
    return r;
}

ReadContext MakeRead(std::uint32_t length = 1000,
                      std::int32_t mapq = 30,
                      std::string_view platform = "hifi") {
    ReadContext c;
    c.read_id = "r1";
    c.read_length = length;
    c.platform = platform;
    c.mapq = mapq;
    return c;
}

}  // namespace

// ===========================================================================
// Lossless floor — Product never falls below kFloor, even when every
// factor is at the per-factor floor.
// ===========================================================================
TEST(LikelihoodFactors, ProductHonoursLosslessFloor) {
    LikelihoodFactors f;
    f.L_sequence = kFloor;
    f.L_modification = kFloor;
    f.L_depth_coverage = kFloor;
    f.L_expression_prior = kFloor;
    f.L_phasing = kFloor;
    f.L_pseudogene_compatibility = kFloor;
    f.L_junction = kFloor;
    f.L_barcode_context = kFloor;
    f.L_mapq_signal = kFloor;
    f.L_length_plausibility = kFloor;
    EXPECT_EQ(f.Product(), kFloor);
}

TEST(LikelihoodFactors, ProductOnesEquals1) {
    LikelihoodFactors f;
    EXPECT_EQ(f.Product(), 1.0f);
}

// ===========================================================================
// L_sequence — platform-dependent.
// ===========================================================================
TEST(LikelihoodFactors, SequencePlatformVarianceHiFiHigherThanOntLegacy) {
    auto a = MakeAnchor();
    auto r_hifi   = MakeRead(500, 30, "hifi");
    auto r_legacy = MakeRead(500, 30, "ont_legacy");
    auto f_hifi   = ComputeFactors(r_hifi,   a, {}, {});
    auto f_legacy = ComputeFactors(r_legacy, a, {}, {});
    EXPECT_GT(f_hifi.L_sequence, f_legacy.L_sequence);
}

// ===========================================================================
// L_mapq_signal — sigmoid, not threshold.
// ===========================================================================
TEST(LikelihoodFactors, MapqSignalIsContinuousNotThresholded) {
    auto a = MakeAnchor();
    auto f0  = ComputeFactors(MakeRead(500,  0), a, {}, {});
    auto f10 = ComputeFactors(MakeRead(500, 10), a, {}, {});
    auto f30 = ComputeFactors(MakeRead(500, 30), a, {}, {});
    EXPECT_NEAR(f10.L_mapq_signal, 0.5f, 0.01f);
    EXPECT_GT(f30.L_mapq_signal, 0.95f);
    // CRITICAL: MAPQ=0 is NOT 0 — it's the lossless floor (0.05).
    EXPECT_GE(f0.L_mapq_signal, 0.05f);
    EXPECT_LT(f0.L_mapq_signal, 0.15f);
}

TEST(LikelihoodFactors, NegativeMapqIsNeutral) {
    // Unknown MAPQ (=-1) ⇒ factor = 1.0 (neutral, no penalty).
    auto a = MakeAnchor();
    auto f = ComputeFactors(MakeRead(500, -1), a, {}, {});
    EXPECT_EQ(f.L_mapq_signal, 1.0f);
}

// ===========================================================================
// L_length_plausibility — anchor-kind-aware.
// ===========================================================================
TEST(LikelihoodFactors, LengthPlausibilityWindowMirna) {
    auto a = MakeAnchor(TranscriptKind::Mirna, std::string(22, 'A'));
    auto f_in  = ComputeFactors(MakeRead(22, 30), a, {}, {});
    auto f_out = ComputeFactors(MakeRead(5000, 30), a, {}, {});
    EXPECT_EQ(f_in.L_length_plausibility,  1.0f);
    EXPECT_LE(f_out.L_length_plausibility, 0.2f);
}

TEST(LikelihoodFactors, LengthPlausibilityWindowMatureMrna) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna);
    auto f_in  = ComputeFactors(MakeRead(2000, 30), a, {}, {});
    auto f_tiny = ComputeFactors(MakeRead(20, 30), a, {}, {});
    EXPECT_EQ(f_in.L_length_plausibility,  1.0f);
    EXPECT_LE(f_tiny.L_length_plausibility, 0.2f);
}

// ===========================================================================
// L_pseudogene_compatibility — biotype-driven.
// ===========================================================================
TEST(LikelihoodFactors, PseudogeneCompatibilityFunctionalIsOne) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(500, 'A'),
                         {"biotype:protein_coding"});
    auto f = ComputeFactors(MakeRead(), a, {}, {});
    EXPECT_EQ(f.L_pseudogene_compatibility, 1.0f);
}

TEST(LikelihoodFactors, PseudogeneCompatibilityIgPseudogeneIsLow) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(500, 'A'),
                         {"biotype:IG_C_pseudogene"});
    auto f = ComputeFactors(MakeRead(), a, {}, {});
    EXPECT_EQ(f.L_pseudogene_compatibility, 0.40f);
}

TEST(LikelihoodFactors, PseudogeneCompatibilityTranscribedHigh) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(500, 'A'),
                         {"biotype:transcribed_processed_pseudogene"});
    auto f = ComputeFactors(MakeRead(), a, {}, {});
    EXPECT_EQ(f.L_pseudogene_compatibility, 0.90f);
}

// ===========================================================================
// L_phasing — hap consistency.
// ===========================================================================
TEST(LikelihoodFactors, PhasingMatchingHaplotypeOne) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(500, 'A'),
                         {"PANGEN_HG00329_hap1"});
    auto r = MakeRead();
    r.haplotype = 0;  // hap1
    auto f = ComputeFactors(r, a, {}, {});
    EXPECT_EQ(f.L_phasing, 1.0f);
}

TEST(LikelihoodFactors, PhasingCrossHaplotypeLowButNonZero) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(500, 'A'),
                         {"PANGEN_HG00329_hap2"});
    auto r = MakeRead();
    r.haplotype = 0;  // hap1 but anchor hap2
    auto f = ComputeFactors(r, a, {}, {});
    EXPECT_EQ(f.L_phasing, 0.05f);
}

TEST(LikelihoodFactors, PhasingUnknownReadHapNeutral) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(500, 'A'),
                         {"PANGEN_HG00329_hap1"});
    auto r = MakeRead();
    // haplotype not set
    auto f = ComputeFactors(r, a, {}, {});
    EXPECT_EQ(f.L_phasing, 1.0f);
}

// ===========================================================================
// L_junction — splice-PWM aggregation.
// ===========================================================================
TEST(LikelihoodFactors, JunctionNoBoundariesIsNeutral) {
    auto a = MakeAnchor();
    auto f = ComputeFactors(MakeRead(), a, {}, {});
    EXPECT_EQ(f.L_junction, 1.0f);
}

TEST(LikelihoodFactors, JunctionHighScoringBoundariesProductMultiplies) {
    auto a = MakeAnchor();
    ExonBoundary b;
    b.donor_score = 0.9f;
    b.acceptor_score = 0.9f;
    a.exon_boundaries.push_back(b);
    a.exon_boundaries.push_back(b);
    auto f = ComputeFactors(MakeRead(), a, {}, {});
    // 0.81 * 0.81 = 0.6561
    EXPECT_NEAR(f.L_junction, 0.6561f, 0.01f);
}

// ===========================================================================
// L_modification — m6A signal via DRACH context.
// ===========================================================================
TEST(LikelihoodFactors, ModificationNoCallsIsNeutral) {
    auto a = MakeAnchor();
    auto f = ComputeFactors(MakeRead(), a, {}, {});
    EXPECT_EQ(f.L_modification, 1.0f);
}

TEST(LikelihoodFactors, ModificationDrachM6AIsRewarded) {
    auto a = MakeAnchor();
    auto r = MakeRead();
    // Use a 100-bp read whose pos 10 sits in a DRACH context.
    // Context "GGACT" at pos-2..pos+2 of the m6A site.
    std::string seq = std::string(100, 'T');
    // place GGACT at positions 8..12 so pos=10 → middle 'A'
    seq[8] = 'G'; seq[9] = 'G'; seq[10] = 'A'; seq[11] = 'C'; seq[12] = 'T';
    r.read_sequence = seq;
    ObservedModificationCalls m;
    m.m6a_calls.push_back({10u, 1.0f});   // high-confidence call
    auto f = ComputeFactors(r, a, m, {});
    EXPECT_GE(f.L_modification, 0.94f);  // reward applied
    EXPECT_LE(f.L_modification, 0.96f);
}

// ===========================================================================
// FactorDisableMask — disable returns 1.0 (neutral).
// ===========================================================================
TEST(LikelihoodFactors, DisableMaskTurnsFactorOff) {
    auto a = MakeAnchor();
    auto r_low_mapq = MakeRead(500, 0);
    FactorDisableMask dis;
    dis.mapq_signal = true;
    auto f = ComputeFactors(r_low_mapq, a, {}, {}, dis);
    EXPECT_EQ(f.L_mapq_signal, 1.0f);
}

// ===========================================================================
// 10-factor smoke test — Product reflects multi-signal multiplication.
// ===========================================================================
TEST(LikelihoodFactors, IntegrationProductFromAllFactorsStaysInBounds) {
    auto a = MakeAnchor(TranscriptKind::MatureMrna,
                         std::string(1000, 'A'),
                         {"biotype:protein_coding"});
    ExonBoundary b;
    b.donor_score = 0.95f;
    b.acceptor_score = 0.85f;
    a.exon_boundaries.push_back(b);

    auto r = MakeRead(1000, 30, "hifi");
    r.haplotype = 0;
    TissueContext t{.label = "lymph"};

    auto f = ComputeFactors(r, a, {}, t);
    const float p = f.Product();
    EXPECT_GE(p, kFloor);
    EXPECT_LE(p, 1.0f);
    EXPECT_GT(p, 0.5f);  // a well-matched scenario should pass solidly above floor
}
