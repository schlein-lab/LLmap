// LLmap — pangenome_mappability tests (M(pos) cross-sample mapping prior).

#include "core/pangenome_mappability.h"

#include <gtest/gtest.h>

#include <filesystem>

namespace {

using llmap::core::PangenomeMappability;

TEST(PangenomeMappability, CrossSampleMeanPerWindow) {
    PangenomeMappability m(100);  // 100 bp windows
    // Three samples in the same window [1000,1100): 0.9, 0.8, 1.0 → M=0.9.
    m.AddSample("chr1", 1050, 0.9f);
    m.AddSample("chr1", 1050, 0.8f);
    m.AddSample("chr1", 1070, 1.0f);
    EXPECT_NEAR(m.PopulationSupport("chr1", 1050), 0.9, 1e-6);
    EXPECT_NEAR(m.PopulationSupport("chr1", 1099), 0.9, 1e-6);  // same window
    EXPECT_EQ(m.SampleCount("chr1", 1010), 3u);
}

TEST(PangenomeMappability, UnseenLocusReturnsNeutral) {
    PangenomeMappability m(100);
    m.AddSample("chr1", 1050, 0.9f);
    EXPECT_FLOAT_EQ(m.PopulationSupport("chr1", 999'999, 0.5f), 0.5f);
    EXPECT_FLOAT_EQ(m.PopulationSupport("chr2", 1050, 0.5f), 0.5f);  // other ref
    EXPECT_EQ(m.SampleCount("chr1", 999'999), 0u);
}

TEST(PangenomeMappability, ReproduciblyBlurryIsLow) {
    PangenomeMappability m(100);
    // A segdup/TE locus blurry across the pangenome → low M (honest, reliable).
    m.AddSample("chr1", 5000, 0.10f);
    m.AddSample("chr1", 5050, 0.15f);
    m.AddSample("chr1", 5099, 0.05f);
    EXPECT_LT(m.PopulationSupport("chr1", 5025), 0.20f);
    EXPECT_EQ(m.SampleCount("chr1", 5025), 3u);
}

TEST(PangenomeMappability, SaveLoadRoundTrip) {
    namespace fs = std::filesystem;
    PangenomeMappability m(100);
    m.AddSample("chr1", 1050, 0.9f);
    m.AddSample("chr1", 1050, 0.8f);
    m.AddSample("chr1", 1070, 1.0f);
    m.AddSample("chr1", 5000, 0.10f);

    const fs::path p = fs::path(::testing::TempDir()) / "mappability.bedgraph";
    ASSERT_TRUE(m.Save(p.string()));

    PangenomeMappability m2(100);
    ASSERT_TRUE(m2.Load(p.string()));
    EXPECT_NEAR(m2.PopulationSupport("chr1", 1050), 0.9, 1e-6);
    EXPECT_EQ(m2.SampleCount("chr1", 1050), 3u);
    EXPECT_NEAR(m2.PopulationSupport("chr1", 5000), 0.10, 1e-6);
    fs::remove(p);
}

TEST(PangenomeMappability, ConfidenceClampedToUnit) {
    PangenomeMappability m(100);
    m.AddSample("chr1", 100, 5.0f);   // clamped to 1.0
    m.AddSample("chr1", 100, -2.0f);  // clamped to 0.0
    EXPECT_NEAR(m.PopulationSupport("chr1", 100), 0.5, 1e-6);  // (1+0)/2
}

}  // namespace
