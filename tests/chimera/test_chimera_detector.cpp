// LLmap — Chimera detector + VDJ-mask tests.

#include "chimera/chimera_detector.h"

#include <gtest/gtest.h>

#include <vector>

using namespace llmap::chimera;

namespace {

AlignedPart MakePart(std::string chrom,
                      std::uint64_t r_start, std::uint64_t r_end,
                      std::uint32_t q_off, std::uint32_t q_len,
                      bool switch_region = false) {
    AlignedPart p;
    p.ref_chrom = std::move(chrom);
    p.ref_start = r_start;
    p.ref_end   = r_end;
    p.read_offset = q_off;
    p.read_length = q_len;
    p.score = 100.0f;
    p.is_switch_region = switch_region;
    return p;
}

VdjLocusMask MakeMask() {
    VdjLocusMask m;
    m.LoadGrch38Defaults();
    return m;
}

}  // namespace

// ===========================================================================
// VdjLocusMask
// ===========================================================================

TEST(VdjLocusMask, GrCh38DefaultsLoadAllSeven) {
    VdjLocusMask m;
    m.LoadGrch38Defaults();
    EXPECT_EQ(m.size(), 7u);
}

TEST(VdjLocusMask, IghPositionResolvesToIGH) {
    auto m = MakeMask();
    EXPECT_EQ(m.LocusAt("chr14", 105'700'000), "IGH");
    EXPECT_EQ(m.LocusAt("chr14", 105'586'437), "IGH");  // start inclusive
    EXPECT_EQ(m.LocusAt("chr14", 106'879'843), "IGH");  // end-1 inclusive
}

TEST(VdjLocusMask, OutsideIghReturnsEmpty) {
    auto m = MakeMask();
    EXPECT_EQ(m.LocusAt("chr14", 105'586'436), "");  // pre-IGH
    EXPECT_EQ(m.LocusAt("chr14", 106'879'844), "");  // post-IGH (end excl)
    EXPECT_EQ(m.LocusAt("chr1",  100'000'000), "");
}

TEST(VdjLocusMask, BothInSameIghLocus) {
    auto m = MakeMask();
    EXPECT_TRUE(m.BothInSameVdjLocus(
        "chr14", 105'700'000, "chr14", 106'800'000));
    EXPECT_FALSE(m.BothInSameVdjLocus(
        "chr14", 105'700'000, "chr14", 50'000'000));  // pos2 outside IGH
}

// ===========================================================================
// Analyze — chimerism kinds
// ===========================================================================

TEST(Analyze, FewerThanTwoPartsReturnsEmpty) {
    auto m = MakeMask();
    std::vector<AlignedPart> parts = {
        MakePart("chr14", 1000, 1500, 0, 500),
    };
    EXPECT_TRUE(Analyze(parts, m, {}).empty());
}

TEST(Analyze, IntraRegionSameChromCloseDistance) {
    auto m = MakeMask();
    std::vector<AlignedPart> parts = {
        MakePart("chr1", 100'000, 100'500,   0, 500),
        MakePart("chr1", 200'000, 200'500, 500, 500),  // 100 kb gap
    };
    auto out = Analyze(parts, m, {});
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].kind, 'I');
    EXPECT_FALSE(out[0].vdj_class_switch_detected);
    EXPECT_GT(out[0].genomic_distance_bp, 99'999u);
}

TEST(Analyze, VdjClassSwitchYieldsKindV) {
    auto m = MakeMask();
    std::vector<AlignedPart> parts = {
        MakePart("chr14", 105'860'000, 105'860'400,   0, 400),  // IGHM-CH1
        MakePart("chr14", 105'625'700, 105'625'994, 400, 294,
                  /*switch_region=*/true),                       // S-region between
    };
    auto out = Analyze(parts, m, {});
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].kind, 'V');
    EXPECT_TRUE(out[0].vdj_class_switch_detected);
}

TEST(Analyze, CrossChromYieldsKindX) {
    auto m = MakeMask();
    std::vector<AlignedPart> parts = {
        MakePart("chr1",  1000, 1500,   0, 500),
        MakePart("chr22", 5000, 5500, 500, 500),
    };
    auto out = Analyze(parts, m, {});
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].kind, 'X');
    EXPECT_EQ(out[0].genomic_distance_bp, 0u);  // not meaningful for cross-chrom
}

TEST(Analyze, CrossChromSuppressibleViaConfig) {
    auto m = MakeMask();
    std::vector<AlignedPart> parts = {
        MakePart("chr1",  1000, 1500,   0, 500),
        MakePart("chr22", 5000, 5500, 500, 500),
    };
    ChimeraConfig cfg;
    cfg.emit_low_prior_chimeric = false;
    auto out = Analyze(parts, m, cfg);
    EXPECT_TRUE(out.empty());
}

TEST(Analyze, BelowMinPartCoverageDropped) {
    auto m = MakeMask();
    // First part 95 %, second part 2 % — second below default 10 %
    std::vector<AlignedPart> parts = {
        MakePart("chr1", 1000, 2000,   0, 950),
        MakePart("chr1", 3000, 3050, 950,  20),
    };
    auto out = Analyze(parts, m, {});
    EXPECT_TRUE(out.empty()) << "minor part below 10% coverage → not chimeric";
}

TEST(Analyze, BelowMinCombinedCoverageDropped) {
    auto m = MakeMask();
    // Two parts, each 30 %, combined 60 % < default 90 %
    std::vector<AlignedPart> parts = {
        MakePart("chr1", 1000, 1300,   0, 300),
        MakePart("chr1", 5000, 5300, 700, 300),
    };
    auto out = Analyze(parts, m, {});
    EXPECT_TRUE(out.empty());
}

// ===========================================================================
// DistancePrior
// ===========================================================================

TEST(DistancePrior, CrossChromIsLowest) {
    auto p = DistancePrior(0, /*cross_chrom=*/true);
    EXPECT_LT(p, -7.0f);
}

TEST(DistancePrior, IntraNearHighesterThanFar) {
    auto p_near = DistancePrior(500'000, false);
    auto p_far  = DistancePrior(40'000'000, false);
    EXPECT_GT(p_near, p_far);
}
