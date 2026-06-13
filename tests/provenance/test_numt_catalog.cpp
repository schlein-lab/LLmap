// LLmap — numt_catalog tests.

#include "provenance/numt_catalog.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

namespace {

using llmap::provenance::NumtCatalog;

TEST(NumtCatalog, BuiltinStarterLookup) {
    NumtCatalog cat;
    cat.LoadBuiltinStarter();
    EXPECT_GE(cat.Size(), 4u);

    const auto* p = cat.Lookup("chr1", 566'000);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->chrom, "chr1");
    EXPECT_FALSE(p->verified);  // starter coords are approximate

    EXPECT_EQ(cat.Lookup("chr1", 1'000'000), nullptr);  // outside any NUMT
    EXPECT_EQ(cat.Lookup("chr2", 566'000), nullptr);     // wrong chrom
}

TEST(NumtCatalog, LoadBedAuthoritative) {
    namespace fs = std::filesystem;
    const fs::path p = fs::path(::testing::TempDir()) / "numts.bed";
    {
        std::ofstream out(p);
        out << "# chrom start end mt_region\n";
        out << "chr9\t1000\t2000\tMT:ND4\n";
        out << "chr9\t5000\t6000\tMT:CYTB\n";
    }
    NumtCatalog cat;
    ASSERT_TRUE(cat.LoadBed(p.string()));
    EXPECT_EQ(cat.Size(), 2u);

    const auto* a = cat.Lookup("chr9", 1500);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->mt_region, "MT:ND4");
    EXPECT_TRUE(a->verified);  // curated catalog is authoritative

    EXPECT_EQ(cat.Lookup("chr9", 3000), nullptr);  // gap
    fs::remove(p);
}

TEST(NumtCatalog, LoadBedMissingFileFails) {
    NumtCatalog cat;
    EXPECT_FALSE(cat.LoadBed("/nonexistent/numts.bed"));
}

}  // namespace
