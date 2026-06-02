// LLmap — Unit tests for curated (T1) SegDup catalog loader.
//
// Loads catalog/curated/*.json from the repo and asserts that the
// well-known IGHG T1 entries are parsed correctly (locus_id,
// structural_architecture, haplotype_class, n_snps).

#include <gtest/gtest.h>

#include "catalog/segdup_catalog.h"

#include <filesystem>

namespace llmap::catalog {
namespace {

class CuratedLoaderTest : public ::testing::Test {
protected:
    SegDupCatalog cat_;

    void SetUp() override {
        auto status = cat_.load_curated_dir(
            std::filesystem::path(LLMAP_CATALOG_CURATED_DIR));
        ASSERT_TRUE(status.ok) << status.error;
        ASSERT_GT(status.records_loaded, 0u);
    }
};

TEST_F(CuratedLoaderTest, LoadsAllCuratedFiles) {
    // The repo currently ships >= 10 curated T1 JSONs (rDNA, HLA, IGHG,
    // FCGR, Y-palindromes, ...). The test pins to "at least 3" so it
    // remains stable while the curated set grows.
    EXPECT_GE(cat_.size(), 3u);
}

TEST_F(CuratedLoaderTest, IGHG4_chimdup_tandem_fields) {
    auto hits = cat_.lookup_by_locus_id("IGHG4_chimdup_tandem");
    ASSERT_EQ(hits.size(), 1u);
    const auto& e = hits.front();
    EXPECT_EQ(e.structural_architecture, "tandem_dup_small");
    EXPECT_EQ(e.haplotype_class, "IGHG4_chimdup_homozygous");
    EXPECT_EQ(e.n_snps(), 5u);
    EXPECT_TRUE(e.is_curated());
    ASSERT_TRUE(e.coords_by_assembly.count("GRCh38"));
    const auto& g = e.coords_by_assembly.at("GRCh38");
    EXPECT_EQ(g.chrom, "chr14");
    EXPECT_EQ(g.start, 105625772);
    EXPECT_EQ(g.end,   105645272);
    EXPECT_EQ(g.strand, '-');
    EXPECT_EQ(e.mapping_primary.kmer_size, 25);
    EXPECT_EQ(e.mapping_primary.max_mismatch, 2);
    EXPECT_GE(e.fallback_chain.size(), 1u);
}

TEST_F(CuratedLoaderTest, IGHG4_chimdup_canonical_arch_fields) {
    auto hits = cat_.lookup_by_locus_id("IGHG4_chimdup_canonical_arch");
    ASSERT_EQ(hits.size(), 1u);
    const auto& e = hits.front();
    EXPECT_EQ(e.structural_architecture, "single_copy");
    EXPECT_EQ(e.haplotype_class, "IGHG4_chimdup_single_copy");
    EXPECT_EQ(e.n_snps(), 2u);
}

TEST_F(CuratedLoaderTest, IGHG_canondup_nahr_block_fields) {
    auto hits = cat_.lookup_by_locus_id("IGHG_canondup_nahr_block");
    ASSERT_EQ(hits.size(), 1u);
    const auto& e = hits.front();
    EXPECT_EQ(e.structural_architecture, "tandem_dup_large");
    EXPECT_EQ(e.haplotype_class, "IGHG_canonical_homozygous");
    EXPECT_EQ(e.n_snps(), 2u);
}

TEST_F(CuratedLoaderTest, RejectsMalformedJson) {
    SegDupCatalog tmp;
    std::string err;
    auto e = parse_curated_json("not json at all", "<inline>", err);
    EXPECT_FALSE(e.has_value());
    EXPECT_FALSE(err.empty());
}

TEST_F(CuratedLoaderTest, RejectsMissingLocusId) {
    SegDupCatalog tmp;
    std::string err;
    const char* doc = R"({
        "structural_architecture": "single_copy",
        "coords": {"GRCh38": {"chrom":"chr1","start":0,"end":100}}
    })";
    auto e = parse_curated_json(doc, "<inline>", err);
    EXPECT_FALSE(e.has_value());
    EXPECT_NE(err.find("locus_id"), std::string::npos);
}

TEST_F(CuratedLoaderTest, SnpFieldsParsed) {
    auto hits = cat_.lookup_by_locus_id("IGHG4_chimdup_tandem");
    ASSERT_EQ(hits.size(), 1u);
    const auto& snps = hits.front().discriminating_snps;
    ASSERT_EQ(snps.size(), 5u);
    // Check the IGHG4_CH1_pos69 SNP (DUP_fixed marker).
    bool found = false;
    for (const auto& s : snps) {
        if (s.id == "IGHG4_CH1_pos69") {
            found = true;
            EXPECT_EQ(s.transcript, "IGHG4-201");
            EXPECT_EQ(s.cds_offset, 69);
            EXPECT_EQ(s.grch38_pos, 105625998);
            EXPECT_EQ(s.ref, "C");
            EXPECT_EQ(s.alt, "G");
            EXPECT_EQ(s.snp_class, "DUP_fixed");
            EXPECT_NEAR(s.freq_dup, 1.0, 1e-6);
            EXPECT_NEAR(s.freq_canonical, 0.0, 1e-6);
        }
    }
    EXPECT_TRUE(found);
}

}  // namespace
}  // namespace llmap::catalog
