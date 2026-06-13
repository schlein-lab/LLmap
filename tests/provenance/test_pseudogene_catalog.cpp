// LLmap — pseudogene_catalog tests.

#include "provenance/pseudogene_catalog.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace {

using llmap::provenance::PseudogeneCatalog;
using llmap::provenance::PseudogeneRole;

TEST(PseudogeneCatalog, BuiltinStarterLookup) {
    PseudogeneCatalog cat;
    cat.LoadBuiltinStarter();
    EXPECT_EQ(cat.Size(), 6u);

    // Inside the GBA1 parent locus.
    const auto* p = cat.Lookup("chr1", 155'240'000);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->parent, "GBA1");
    EXPECT_EQ(p->role, PseudogeneRole::Parent);
    EXPECT_FALSE(p->verified);  // starter coords are approximate

    // Inside the GBAP1 pseudogene locus.
    const auto* q = cat.Lookup("chr1", 155'220'000);
    ASSERT_NE(q, nullptr);
    EXPECT_EQ(q->pseudogene, "GBAP1");
    EXPECT_EQ(q->role, PseudogeneRole::Pseudogene);

    // Outside any pair.
    EXPECT_EQ(cat.Lookup("chr1", 1'000'000), nullptr);
    EXPECT_EQ(cat.Lookup("chr2", 155'240'000), nullptr);  // wrong chrom
}

TEST(PseudogeneCatalog, LoadBedAuthoritative) {
    namespace fs = std::filesystem;
    const fs::path p = fs::path(::testing::TempDir()) / "pseudo_pairs.bed";
    {
        std::ofstream out(p);
        out << "# chrom start end parent pseudogene role\n";
        out << "chrX\t1000\t2000\tFOO\tFOOP1\tparent\n";
        out << "chrX\t5000\t6000\tFOO\tFOOP1\tpseudogene\n";
    }
    PseudogeneCatalog cat;
    ASSERT_TRUE(cat.LoadBed(p.string()));
    EXPECT_EQ(cat.Size(), 2u);

    const auto* a = cat.Lookup("chrX", 1500);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->role, PseudogeneRole::Parent);
    EXPECT_TRUE(a->verified);  // curated catalog is authoritative

    const auto* b = cat.Lookup("chrX", 5500);
    ASSERT_NE(b, nullptr);
    EXPECT_EQ(b->role, PseudogeneRole::Pseudogene);

    EXPECT_EQ(cat.Lookup("chrX", 3000), nullptr);  // gap between the two
    fs::remove(p);
}

TEST(PseudogeneCatalog, LoadBedMissingFileFails) {
    PseudogeneCatalog cat;
    EXPECT_FALSE(cat.LoadBed("/nonexistent/pseudo.bed"));
}

}  // namespace
