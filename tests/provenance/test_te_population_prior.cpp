// Unit tests for the pangenome M(pos) prior → TE spread-mass bridge.

#include "provenance/te_population_prior.h"

#include <gtest/gtest.h>

namespace llmap::provenance {
namespace {

TePlacement P(std::uint64_t pos, std::int32_t score = 100) {
    TePlacement p;
    p.chrom = "chr1";
    p.pos = pos;
    p.align_score = score;
    p.has_unique_anchor = false;
    p.te_class = TeClass::Alu;
    return p;
}

TEST(TePopulationPrior, FillsSupportFromTrackAndShiftsMass) {
    // Two equal-alignment Alu copies. The pangenome track says copy A (pos 1000)
    // is reproducibly sharp across samples, copy B (pos 5000) reproducibly blurry.
    core::PangenomeMappability track(100);
    for (int s = 0; s < 10; ++s) {
        track.AddSample("chr1", 1000, 0.95f);   // sharp everywhere
        track.AddSample("chr1", 5000, 0.05f);   // blurry everywhere
    }
    std::vector<TePlacement> cands = {P(1000), P(5000)};
    ApplyPopulationPrior(cands, track);

    EXPECT_GT(cands[0].population_support, 0.9f);
    EXPECT_LT(cands[1].population_support, 0.1f);

    const auto r = ResolveTeSpreadMass(cands);
    EXPECT_GT(r.posterior[0], r.posterior[1]);   // mass moved to the sharp copy
    EXPECT_GT(r.posterior[0], 0.7f);
}

TEST(TePopulationPrior, UnseenWindowKeepsSentinelByDefault) {
    // No pangenome data at these loci → leave the te_spread_mass "no baseline"
    // sentinel (<0) so behaviour falls back to reference-only honest 50/50.
    core::PangenomeMappability track(100);   // empty
    std::vector<TePlacement> cands = {P(1000), P(5000)};
    ApplyPopulationPrior(cands, track);

    EXPECT_LT(cands[0].population_support, 0.0f);
    EXPECT_LT(cands[1].population_support, 0.0f);

    const auto r = ResolveTeSpreadMass(cands);
    EXPECT_NEAR(r.posterior[0], 0.5f, 1e-4);     // honest even spread, no fake prior
}

TEST(TePopulationPrior, ResolveWithPopulationPriorConvenience) {
    core::PangenomeMappability track(100);
    for (int s = 0; s < 8; ++s) {
        track.AddSample("chr1", 1000, 0.9f);
        track.AddSample("chr1", 5000, 0.1f);
    }
    const auto r = ResolveWithPopulationPrior({P(1000), P(5000)}, track);
    EXPECT_GT(r.posterior[0], r.posterior[1]);
}

TEST(TePopulationPrior, FlagsBlurryOnlyHereAnomaly) {
    // The pangenome says pos 1000 is reproducibly SHARP (M=0.95). This read,
    // however, did not anchor there (no unique anchor) → sample-specific anomaly,
    // NOT a real ambiguity. pos 5000 is reproducibly blurry → expected, not flagged.
    core::PangenomeMappability track(100);
    for (int s = 0; s < 10; ++s) {
        track.AddSample("chr1", 1000, 0.95f);
        track.AddSample("chr1", 5000, 0.10f);
    }
    std::vector<TePlacement> cands = {P(1000), P(5000)};
    const auto anomalies = FlagBlurryOnlyHere(cands, track, 0.8f);
    ASSERT_EQ(anomalies.size(), 1u);
    EXPECT_EQ(anomalies[0], 0);                  // only the sharp-locus mismatch
}

TEST(TePopulationPrior, AnchoredPlacementNotFlagged) {
    // Same sharp locus, but the read DOES anchor uniquely → legitimate, not an
    // anomaly even though M is high.
    core::PangenomeMappability track(100);
    for (int s = 0; s < 10; ++s) track.AddSample("chr1", 1000, 0.95f);
    auto c = P(1000);
    c.has_unique_anchor = true;
    const auto anomalies = FlagBlurryOnlyHere({c}, track, 0.8f);
    EXPECT_TRUE(anomalies.empty());
}

}  // namespace
}  // namespace llmap::provenance
