// LLmap — damage_profile substitution classifier unit tests.
//
// Operator taxonomy class B: 8-oxoG (G>T), FFPE/ancient deamination (C>T),
// ADAR A-to-I (A>G), APOBEC C-to-U (C>T in TC). Per-substitution evidence;
// flags, never filters.

#include "rnamod/damage_profile.h"

#include <gtest/gtest.h>

#include <string>

namespace {

using llmap::rnamod::ClassifySubstitution;
using llmap::rnamod::SubstitutionContext;
using llmap::rnamod::SubstitutionProvenance;
using llmap::rnamod::SubstitutionProvenanceName;

SubstitutionContext Ctx(char ref, char alt, bool rna, bool rev,
                        std::string_view c5 = "NNNNN",
                        std::uint32_t dist = 1000) {
    SubstitutionContext c;
    c.ref = ref;
    c.alt = alt;
    c.is_rna = rna;
    c.read_reverse = rev;
    c.ctx5 = c5;
    c.dist_from_read_end = dist;
    return c;
}

TEST(DamageProfile, OxoGFromGtoTAndCtoA) {
    auto a = ClassifySubstitution(Ctx('G', 'T', /*rna=*/false, /*rev=*/false));
    EXPECT_EQ(a.kind, SubstitutionProvenance::Damage8oxoG);
    EXPECT_TRUE(a.strand_biased);
    // Complementary lesion on the other strand: C>A.
    auto b = ClassifySubstitution(Ctx('C', 'A', false, false));
    EXPECT_EQ(b.kind, SubstitutionProvenance::Damage8oxoG);
}

TEST(DamageProfile, FfpeDeaminationMidRead) {
    auto a = ClassifySubstitution(Ctx('C', 'T', /*rna=*/false, false, "NNNNN",
                                      /*dist=*/50));
    EXPECT_EQ(a.kind, SubstitutionProvenance::DamageFfpeDeam);
    EXPECT_TRUE(a.strand_biased);
}

TEST(DamageProfile, AncientDeaminationAtReadEnd) {
    auto a = ClassifySubstitution(Ctx('C', 'T', /*rna=*/false, false, "NNNNN",
                                      /*dist=*/2));
    EXPECT_EQ(a.kind, SubstitutionProvenance::DamageAncientDeam);
}

TEST(DamageProfile, AdarA2IForwardAndReverse) {
    // Forward sense A>G.
    auto fwd = ClassifySubstitution(Ctx('A', 'G', /*rna=*/true, /*rev=*/false));
    EXPECT_EQ(fwd.kind, SubstitutionProvenance::EditAdarA2I);
    // '-'-strand read: genomic T>C is sense A>G after complementing.
    auto rev = ClassifySubstitution(Ctx('T', 'C', /*rna=*/true, /*rev=*/true));
    EXPECT_EQ(rev.kind, SubstitutionProvenance::EditAdarA2I);
}

TEST(DamageProfile, ApobecC2UInTcContext) {
    // ctx5 centred on the C (index 2), 5' neighbour T (index 1) → TC.
    auto a = ClassifySubstitution(Ctx('C', 'T', /*rna=*/true, false, "ATCGA"));
    EXPECT_EQ(a.kind, SubstitutionProvenance::EditApobecC2U);
    // Without the TC context it is not APOBEC (falls through; no signature here).
    auto b = ClassifySubstitution(Ctx('C', 'T', /*rna=*/true, false, "AGCGA"));
    EXPECT_NE(b.kind, SubstitutionProvenance::EditApobecC2U);
}

TEST(DamageProfile, TransversionIsVariantCandidate) {
    auto a = ClassifySubstitution(Ctx('A', 'C', /*rna=*/false, false));
    EXPECT_EQ(a.kind, SubstitutionProvenance::VariantCandidate);
    auto same = ClassifySubstitution(Ctx('A', 'A', false, false));
    EXPECT_EQ(same.kind, SubstitutionProvenance::VariantCandidate);
}

TEST(DamageProfile, TagNamesMatchProvenanceScheme) {
    EXPECT_STREQ(SubstitutionProvenanceName(SubstitutionProvenance::Damage8oxoG),
                 "dmg:8oxoG");
    EXPECT_STREQ(SubstitutionProvenanceName(SubstitutionProvenance::DamageFfpeDeam),
                 "dmg:ffpe");
    EXPECT_STREQ(
        SubstitutionProvenanceName(SubstitutionProvenance::DamageAncientDeam),
        "dmg:deam");
    EXPECT_STREQ(SubstitutionProvenanceName(SubstitutionProvenance::EditAdarA2I),
                 "edit:adar");
    EXPECT_STREQ(SubstitutionProvenanceName(SubstitutionProvenance::EditApobecC2U),
                 "edit:apobec");
}

}  // namespace
