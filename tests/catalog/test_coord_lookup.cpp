// LLmap — Unit tests for SegDupCatalog coordinate lookup.

#include <gtest/gtest.h>

#include "catalog/segdup_catalog.h"

#include <filesystem>

namespace llmap::catalog {
namespace {

class CoordLookupTest : public ::testing::Test {
protected:
    SegDupCatalog cat_;

    void SetUp() override {
        auto status = cat_.load_curated_dir(
            std::filesystem::path(LLMAP_CATALOG_CURATED_DIR));
        ASSERT_TRUE(status.ok) << status.error;
    }
};

TEST_F(CoordLookupTest, IGHG4_locus_lookup_at_CH1) {
    // Position 105625900 falls inside three nested curated entries:
    //   IGHG_canondup_nahr_block (~300 kb)
    //   IGHG4_chimdup_tandem      (~19.5 kb)
    //   IGHG4_chimdup_canonical_arch (294 bp CH1 only)
    // lookup_by_coords returns the most specific (smallest span). At
    // 105625900 that is canonical_arch.
    auto hit = cat_.lookup_by_coords("GRCh38", "chr14", 105625900);
    ASSERT_TRUE(hit.has_value());
    EXPECT_EQ(hit->locus_id, "IGHG4_chimdup_canonical_arch");

    // lookup_overlapping must surface all three.
    auto all = cat_.lookup_overlapping("GRCh38", "chr14",
                                       105625900, 105625901);
    bool saw_tandem = false, saw_nahr = false, saw_canon = false;
    for (const auto& h : all) {
        if (h.locus_id == "IGHG4_chimdup_tandem")           saw_tandem = true;
        if (h.locus_id == "IGHG_canondup_nahr_block")       saw_nahr   = true;
        if (h.locus_id == "IGHG4_chimdup_canonical_arch")   saw_canon  = true;
    }
    EXPECT_TRUE(saw_tandem);
    EXPECT_TRUE(saw_nahr);
    EXPECT_TRUE(saw_canon);
}

TEST_F(CoordLookupTest, IGHG4_chimdup_tandem_unique_region) {
    // 105640000 is inside tandem (start 105625772, end 105645272) but
    // past canonical_arch.end (105626066). Most-specific match must be
    // tandem.
    auto hit = cat_.lookup_by_coords("GRCh38", "chr14", 105640000);
    ASSERT_TRUE(hit.has_value());
    EXPECT_EQ(hit->locus_id, "IGHG4_chimdup_tandem");
}

TEST_F(CoordLookupTest, IGHG4_chimdup_tandem_boundary) {
    // 105625772 is start (inclusive) of both tandem and canonical_arch.
    // Most specific = canonical_arch (shorter span).
    auto hit_start = cat_.lookup_by_coords("GRCh38", "chr14", 105625772);
    ASSERT_TRUE(hit_start.has_value());
    EXPECT_EQ(hit_start->locus_id, "IGHG4_chimdup_canonical_arch");
    // 105645272 is end (exclusive) of tandem — must NOT return tandem.
    auto hit_end = cat_.lookup_by_coords("GRCh38", "chr14", 105645272);
    if (hit_end.has_value()) {
        EXPECT_NE(hit_end->locus_id, "IGHG4_chimdup_tandem");
        EXPECT_NE(hit_end->locus_id, "IGHG4_chimdup_canonical_arch");
    }
}

TEST_F(CoordLookupTest, OutsideAssemblyMisses) {
    auto hit = cat_.lookup_by_coords("GRCh38", "chr1", 1000);
    EXPECT_FALSE(hit.has_value());
}

TEST_F(CoordLookupTest, UnknownAssemblyReturnsNone) {
    auto hit = cat_.lookup_by_coords("hg19", "chr14", 105625900);
    EXPECT_FALSE(hit.has_value());
}

TEST_F(CoordLookupTest, OverlappingReturnsMultiple) {
    // 105625800..105645000 spans IGHG4_chimdup_tandem; if the larger
    // NAHR block overlaps this window it should appear too. We assert
    // at least 1 hit (the tandem itself).
    auto hits = cat_.lookup_overlapping("GRCh38", "chr14",
                                        105625800, 105645000);
    ASSERT_FALSE(hits.empty());
    bool saw_tandem = false;
    for (const auto& h : hits) {
        if (h.locus_id == "IGHG4_chimdup_tandem") saw_tandem = true;
    }
    EXPECT_TRUE(saw_tandem);
}

TEST_F(CoordLookupTest, SyntheticAddEntryIsIndexed) {
    SegDupCatalog c;
    SegDupCatalogEntry e;
    e.locus_id = "SYNTH_locus_A";
    e.structural_architecture = "single_copy";
    e.tier = SegDupCatalogEntry::Tier::T2_Bulk;
    GenomicCoords g;
    g.assembly = "GRCh38";
    g.chrom = "chrX";
    g.start = 1000;
    g.end   = 2000;
    e.coords_by_assembly.emplace("GRCh38", g);
    c.add_entry(e);

    auto hit = c.lookup_by_coords("GRCh38", "chrX", 1500);
    ASSERT_TRUE(hit.has_value());
    EXPECT_EQ(hit->locus_id, "SYNTH_locus_A");

    auto miss = c.lookup_by_coords("GRCh38", "chrX", 2000);
    EXPECT_FALSE(miss.has_value());
}

TEST_F(CoordLookupTest, LookupByLocusIdRoundTrip) {
    auto hits = cat_.lookup_by_locus_id("IGHG4_chimdup_canonical_arch");
    ASSERT_EQ(hits.size(), 1u);
    EXPECT_EQ(hits.front().locus_id, "IGHG4_chimdup_canonical_arch");

    auto missing = cat_.lookup_by_locus_id("does_not_exist");
    EXPECT_TRUE(missing.empty());
}

}  // namespace
}  // namespace llmap::catalog
