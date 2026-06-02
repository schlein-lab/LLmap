// LLmap — Unit tests for bulk (T2) JSONL SegDup catalog loader.
//
// Writes a small synthetic v2026.Q2-style JSONL into a tmp dir and
// asserts record counts per structural_architecture. If the repo ships
// the canonical catalog/bulk/v2026.Q2.grch38.jsonl, it is additionally
// smoke-loaded.

#include <gtest/gtest.h>

#include "catalog/segdup_catalog.h"

#include <filesystem>
#include <fstream>

namespace llmap::catalog {
namespace {

class BulkLoaderTest : public ::testing::Test {
protected:
    std::filesystem::path tmp_dir_;

    void SetUp() override {
        tmp_dir_ = std::filesystem::temp_directory_path() /
                   "llmap_catalog_bulk_test";
        std::filesystem::create_directories(tmp_dir_);
    }
    void TearDown() override {
        std::error_code ec;
        std::filesystem::remove_all(tmp_dir_, ec);
    }

    std::filesystem::path write_jsonl(const std::string& name,
                                      const std::string& content) {
        auto p = tmp_dir_ / name;
        std::ofstream out(p);
        out << content;
        return p;
    }
};

TEST_F(BulkLoaderTest, LoadsSyntheticJsonl) {
    const std::string jsonl = R"(# v2026.Q2 synthetic bulk catalog
{"locus_id":"SD_001","structural_architecture":"tandem_dup_small","coords":{"GRCh38":{"chrom":"chr1","start":1000,"end":2000}}}
{"locus_id":"SD_002","structural_architecture":"tandem_dup_small","coords":{"GRCh38":{"chrom":"chr1","start":3000,"end":4000}}}
{"locus_id":"SD_003","structural_architecture":"tandem_dup_large","coords":{"GRCh38":{"chrom":"chr2","start":5000,"end":15000}}}
{"locus_id":"SD_004","structural_architecture":"intrachromosomal_interspersed","coords":{"GRCh38":{"chrom":"chr3","start":7000,"end":7500}}}

{"locus_id":"SD_005","structural_architecture":"single_copy","coords":{"GRCh38":{"chrom":"chr4","start":100,"end":200}}}
)";
    auto path = write_jsonl("v2026.Q2.grch38.jsonl", jsonl);

    SegDupCatalog c;
    auto st = c.load_bulk_jsonl(path);
    ASSERT_TRUE(st.ok) << st.error;
    EXPECT_EQ(st.records_loaded, 5u);
    EXPECT_EQ(st.records_skipped, 0u);
    EXPECT_EQ(c.size(), 5u);

    EXPECT_EQ(c.count_architecture("tandem_dup_small"), 2u);
    EXPECT_EQ(c.count_architecture("tandem_dup_large"), 1u);
    EXPECT_EQ(c.count_architecture("intrachromosomal_interspersed"), 1u);
    EXPECT_EQ(c.count_architecture("single_copy"), 1u);

    auto hit = c.lookup_by_coords("GRCh38", "chr1", 1500);
    ASSERT_TRUE(hit.has_value());
    EXPECT_EQ(hit->locus_id, "SD_001");
    EXPECT_TRUE(hit->is_bulk());
}

TEST_F(BulkLoaderTest, SkipsMalformedLinesAndCommentLines) {
    const std::string jsonl = R"(# header comment
{"locus_id":"OK_A","structural_architecture":"single_copy","coords":{"GRCh38":{"chrom":"chr1","start":0,"end":10}}}
not-json
{"locus_id":"","structural_architecture":"single_copy","coords":{"GRCh38":{"chrom":"chr1","start":0,"end":10}}}
{"locus_id":"OK_B","structural_architecture":"single_copy","coords":{"GRCh38":{"chrom":"chr2","start":0,"end":20}}}
)";
    auto path = write_jsonl("malformed.jsonl", jsonl);

    SegDupCatalog c;
    auto st = c.load_bulk_jsonl(path);
    EXPECT_TRUE(st.ok);
    EXPECT_EQ(st.records_loaded, 2u);
    EXPECT_EQ(st.records_skipped, 2u);  // not-json + empty locus_id
}

TEST_F(BulkLoaderTest, LoadsCanonicalBulkIfPresent) {
    // Smoke: if catalog/bulk/v2026.Q2.grch38.jsonl exists in the repo,
    // it must at least open and parse without crashing. We do not assert
    // exact counts because the bulk file is data-driven.
    std::filesystem::path canon =
        std::filesystem::path(LLMAP_CATALOG_BULK_DIR) / "v2026.Q2.grch38.jsonl";
    if (!std::filesystem::exists(canon)) {
        GTEST_SKIP() << "no canonical bulk JSONL shipped at " << canon;
    }
    SegDupCatalog c;
    auto st = c.load_bulk_jsonl(canon);
    EXPECT_TRUE(st.ok);
    // If it has any record, lookup_by_coords must not crash on a sample
    // record.
    if (c.size() > 0) {
        const auto& e0 = c.entries().front();
        ASSERT_FALSE(e0.coords_by_assembly.empty());
        const auto& [asm_name, coords] = *e0.coords_by_assembly.begin();
        auto hit = c.lookup_by_coords(asm_name, coords.chrom, coords.start);
        EXPECT_TRUE(hit.has_value());
    }
}

}  // namespace
}  // namespace llmap::catalog
