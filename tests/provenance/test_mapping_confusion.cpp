// LLmap — mapping_confusion provenance detector tests (Layer-1 origin classes).
//
// Class C: paralog/segdup mis-placement, NUMT artifact (fake heteroplasmy)
// vs genuine mt heteroplasmy (Layer-3, kept separate), processed pseudogene,
// rDNA/satellite. Pure classifier over an explicit evidence struct.

#include "provenance/mapping_confusion.h"

#include <gtest/gtest.h>

namespace {

using llmap::provenance::ClassifyMappingConfusion;
using llmap::provenance::IsRealMtHeteroplasmy;
using llmap::provenance::MappingConfusion;
using llmap::provenance::MappingConfusionEvidence;
using llmap::provenance::MappingConfusionTag;

TEST(MappingConfusion, PseudogeneIntronlessAtParentLocus) {
    MappingConfusionEvidence e;
    e.at_pseudogene_parent_locus = true;
    e.read_is_spliced = false;  // intronless where the parent has introns
    EXPECT_EQ(ClassifyMappingConfusion(e).kind, MappingConfusion::Pseudogene);

    // The spliced (parent) read at the same locus is NOT a pseudogene call.
    e.read_is_spliced = true;
    EXPECT_NE(ClassifyMappingConfusion(e).kind, MappingConfusion::Pseudogene);
}

TEST(MappingConfusion, NumtArtifactWhenNuclearWins) {
    MappingConfusionEvidence e;
    e.at_mt_homologous_locus = true;
    e.identity_to_mt = 0.95f;
    e.identity_to_nuclear_numt = 0.99f;  // nuclear NUMT clearly better → artifact
    auto c = ClassifyMappingConfusion(e);
    EXPECT_EQ(c.kind, MappingConfusion::Numt);
    EXPECT_GT(c.confidence, 0.0f);
    EXPECT_FALSE(IsRealMtHeteroplasmy(e));  // mutually exclusive
}

TEST(MappingConfusion, RealMtHeteroplasmyNotNumt) {
    MappingConfusionEvidence e;
    e.at_mt_homologous_locus = true;
    e.identity_to_mt = 0.99f;            // genuine mt read
    e.identity_to_nuclear_numt = 0.94f;
    EXPECT_TRUE(IsRealMtHeteroplasmy(e));
    // Layer-1 must NOT bucket it as Numt — it is host mt biology (Layer-3 mthet).
    EXPECT_NE(ClassifyMappingConfusion(e).kind, MappingConfusion::Numt);
}

TEST(MappingConfusion, ParalogWhenPsvAmbiguous) {
    MappingConfusionEvidence e;
    e.has_psv = true;
    e.psv_posterior = 0.55f;  // < 0.90 threshold → ambiguous copy
    auto c = ClassifyMappingConfusion(e);
    EXPECT_EQ(c.kind, MappingConfusion::Paralog);
    EXPECT_NEAR(c.confidence, 0.45f, 1e-5f);

    e.psv_posterior = 0.99f;  // confident → not paralog-confused
    EXPECT_EQ(ClassifyMappingConfusion(e).kind, MappingConfusion::None);
}

TEST(MappingConfusion, RdnaRepeatLowMapq) {
    MappingConfusionEvidence e;
    e.in_repeat_array = true;
    e.mapq = 0;
    EXPECT_EQ(ClassifyMappingConfusion(e).kind, MappingConfusion::Rdna);
}

TEST(MappingConfusion, CleanReadIsNone) {
    MappingConfusionEvidence e;  // all defaults → uniquely placed
    EXPECT_EQ(ClassifyMappingConfusion(e).kind, MappingConfusion::None);
}

TEST(MappingConfusion, Tags) {
    EXPECT_STREQ(MappingConfusionTag(MappingConfusion::Paralog), "para");
    EXPECT_STREQ(MappingConfusionTag(MappingConfusion::Numt), "numt");
    EXPECT_STREQ(MappingConfusionTag(MappingConfusion::Pseudogene), "pseudo");
    EXPECT_STREQ(MappingConfusionTag(MappingConfusion::Rdna), "rdna");
    EXPECT_STREQ(MappingConfusionTag(MappingConfusion::None), "none");
}

}  // namespace
