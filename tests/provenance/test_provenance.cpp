// Unit tests for the three-layer read provenance model + spectrum.

#include "provenance/provenance_class.h"
#include "provenance/contamination_spectrum.h"

#include <gtest/gtest.h>

namespace llmap::provenance {
namespace {

TEST(ProvenanceClass, TagParseRoundTrip) {
    for (std::size_t i = 0; i < static_cast<std::size_t>(ProvenanceClass::Count); ++i) {
        const auto c = static_cast<ProvenanceClass>(i);
        const auto parsed = ParseProvenanceClass(ProvenanceClassTag(c));
        ASSERT_TRUE(parsed.has_value());
        EXPECT_EQ(*parsed, c);
    }
    EXPECT_FALSE(ParseProvenanceClass("nonsense").has_value());
    EXPECT_STREQ(ProvenanceClassTag(ProvenanceClass::Host), "host");
    EXPECT_STREQ(ProvenanceClassTag(ProvenanceClass::RefArtefact), "refartefact");
    EXPECT_STREQ(ProvenanceClassTag(ProvenanceClass::CrossIndividual), "xindiv");
}

TEST(ProvenanceClass, OriginBucketFamily) {
    // A: competing-reference origin buckets.
    EXPECT_TRUE(IsOriginBucket(ProvenanceClass::Exogenous));
    EXPECT_TRUE(IsOriginBucket(ProvenanceClass::CrossIndividual));
    EXPECT_TRUE(IsOriginBucket(ProvenanceClass::Numt));
    EXPECT_TRUE(IsOriginBucket(ProvenanceClass::MobileElement));
    // partition members that are NOT "another sequence" buckets.
    EXPECT_FALSE(IsOriginBucket(ProvenanceClass::Host));
    EXPECT_FALSE(IsOriginBucket(ProvenanceClass::RefArtefact));
    EXPECT_FALSE(IsOriginBucket(ProvenanceClass::Chimera));
    EXPECT_FALSE(IsOriginBucket(ProvenanceClass::Multiplicity));
}

TEST(ContaminationSpectrum, Layer1PartitionLosslessSigma) {
    ContaminationSpectrum sp;
    for (int i = 0; i < 96; ++i) sp.Add({ProvenanceClass::Host, 1.0f, 150, "", 0});
    for (int i = 0; i < 3; ++i)  sp.Add({ProvenanceClass::Exogenous, 0.92f, 150, "exo:phix", 0});
    sp.Add({ProvenanceClass::Numt, 0.80f, 150, "numt", 0});

    EXPECT_EQ(sp.TotalReads(), 100u);
    EXPECT_TRUE(sp.CheckLossless(100));      // Σ over origin classes == N
    EXPECT_FALSE(sp.CheckLossless(101));

    const ClassStat exo = sp.Stat(ProvenanceClass::Exogenous);
    EXPECT_EQ(exo.n_reads, 3u);
    EXPECT_DOUBLE_EQ(exo.fraction, 0.03);
    EXPECT_EQ(exo.bases, 450u);
    EXPECT_NEAR(exo.mean_posterior, 0.92, 1e-5);
    EXPECT_DOUBLE_EQ(sp.Stat(ProvenanceClass::Host).fraction, 0.96);
}

TEST(ContaminationSpectrum, Layer3OverlayDoesNotBreakSigma) {
    ContaminationSpectrum sp;
    // 10 host reads; 4 of them carry V(D)J biology (overlay), 1 also mt-heteroplasmy.
    for (int i = 0; i < 6; ++i) sp.Add({ProvenanceClass::Host, 1.0f, 150, "", 0});
    for (int i = 0; i < 3; ++i)
        sp.Add({ProvenanceClass::Host, 1.0f, 150, "",
                static_cast<std::uint16_t>(BioConfounder::Vdj)});
    sp.Add({ProvenanceClass::Host, 1.0f, 150, "",
            static_cast<std::uint16_t>(BioConfounder::Vdj) |
            static_cast<std::uint16_t>(BioConfounder::MtHeteroplasmy)});

    // Partition: all 10 are Host — overlay flags do NOT move them out.
    EXPECT_EQ(sp.TotalReads(), 10u);
    EXPECT_TRUE(sp.CheckLossless(10));
    EXPECT_EQ(sp.Stat(ProvenanceClass::Host).n_reads, 10u);
    // Overlay counts (orthogonal, may overlap Host).
    EXPECT_EQ(sp.BioConfounderReads(BioConfounder::Vdj), 4u);
    EXPECT_EQ(sp.BioConfounderReads(BioConfounder::MtHeteroplasmy), 1u);
    EXPECT_EQ(sp.BioConfounderReads(BioConfounder::Shm), 0u);
}

TEST(ContaminationSpectrum, ToStringHasBothLayers) {
    ContaminationSpectrum sp;
    sp.Add({ProvenanceClass::Host, 1.0f, 100, "",
            static_cast<std::uint16_t>(BioConfounder::Vdj)});
    const std::string s = sp.ToString();
    EXPECT_NE(s.find("layer1_origin_partition"), std::string::npos);
    EXPECT_NE(s.find("host"), std::string::npos);
    EXPECT_NE(s.find("layer3_bioconfounder_overlay"), std::string::npos);
    EXPECT_NE(s.find("bio:vdj"), std::string::npos);
}

}  // namespace
}  // namespace llmap::provenance
