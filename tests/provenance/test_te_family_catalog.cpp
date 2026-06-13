// Unit tests for the TE family catalog (Block 2).

#include "provenance/te_family_catalog.h"

#include <gtest/gtest.h>

namespace llmap::provenance {
namespace {

TEST(TeFamilyCatalog, ClassTags) {
    EXPECT_STREQ(TeClassTag(TeClass::Alu), "alu");
    EXPECT_STREQ(TeClassTag(TeClass::Line1), "l1");
    EXPECT_STREQ(TeClassTag(TeClass::Sva), "sva");
    EXPECT_STREQ(TeClassTag(TeClass::Herv), "herv");
}

TEST(TeFamilyCatalog, BuiltinStarterLookup) {
    TeFamilyCatalog cat;
    cat.LoadBuiltinStarter();
    EXPECT_GT(cat.Size(), 0u);

    // Inside the chr1 AluYa5 locus [145000000,145000300).
    const TeLocus* a = cat.Lookup("chr1", 145000150);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->cls, TeClass::Alu);
    EXPECT_EQ(a->subfamily, "AluYa5");

    // Inside the chr1 L1HS locus.
    const TeLocus* l1 = cat.Lookup("chr1", 145103000);
    ASSERT_NE(l1, nullptr);
    EXPECT_EQ(l1->cls, TeClass::Line1);
    EXPECT_EQ(l1->subfamily, "L1HS");

    // Other chromosome's Alu.
    const TeLocus* a7 = cat.Lookup("chr7", 142000100);
    ASSERT_NE(a7, nullptr);
    EXPECT_EQ(a7->cls, TeClass::Alu);
}

TEST(TeFamilyCatalog, GapAndWrongChromAreUnique) {
    TeFamilyCatalog cat;
    cat.LoadBuiltinStarter();
    EXPECT_EQ(cat.Lookup("chr1", 145000500), nullptr);   // gap between loci
    EXPECT_EQ(cat.Lookup("chr1", 1000), nullptr);        // before any locus
    EXPECT_EQ(cat.Lookup("chrX", 145000150), nullptr);   // chromosome with no TE
    EXPECT_EQ(cat.Lookup("chr1", 145000300), nullptr);   // exactly at end (half-open)
}

TEST(TeFamilyCatalog, BoundaryHalfOpen) {
    TeFamilyCatalog cat;
    cat.LoadBuiltinStarter();
    EXPECT_NE(cat.Lookup("chr1", 145000000), nullptr);   // start is inclusive
    EXPECT_NE(cat.Lookup("chr1", 145000299), nullptr);   // end-1 inside
    EXPECT_EQ(cat.Lookup("chr1", 145000300), nullptr);   // end is exclusive
}

}  // namespace
}  // namespace llmap::provenance
