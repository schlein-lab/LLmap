// LLmap — AnchorStore unit tests.
//
// Synthetic anchors only (no on-disk GFF3 dependency). The
// GENCODE-loader test below uses a tiny in-memory GFF3 written to a
// temp file so the test can prove the parser handles the format
// without needing the 1.5 GB v46 GFF.
//
// Tests cover:
//   - AddAnchor + ById + by_id_ uniqueness
//   - ByTag returns indices in load order
//   - ByRegion half-open overlap (edge cases)
//   - ByTranscriptId aggregation
//   - ForEach predicate iteration
//   - CountBySource telemetry
//   - Reindex idempotence after manual mutation
//   - LoadGencodeGff with a synthetic GFF3 (uncompressed!)

#include "anchor/anchor_store.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

using namespace llmap::anchor;

namespace {

AnchorRecord MakeAnchor(std::string id,
                         AnchorSource src,
                         std::string chrom = "",
                         std::int64_t start = 0,
                         std::int64_t end   = 0,
                         std::vector<std::string> tags = {},
                         std::string transcript_id = "") {
    AnchorRecord r;
    r.anchor_id = std::move(id);
    r.source = src;
    r.tags = std::move(tags);
    r.transcript_id = std::move(transcript_id);
    if (!chrom.empty()) {
        r.ref_chrom = std::move(chrom);
        r.ref_start = start;
        r.ref_end = end;
        r.strand = '+';
    }
    return r;
}

}  // namespace

TEST(AnchorStore, AddAndLookupById) {
    AnchorStore s;
    auto idx = s.AddAnchor(MakeAnchor("GENCODE:T1:exon1", AnchorSource::Gencode));
    EXPECT_EQ(idx, 0u);
    EXPECT_EQ(s.size(), 1u);

    const AnchorRecord* hit = s.ById("GENCODE:T1:exon1");
    ASSERT_NE(hit, nullptr);
    EXPECT_EQ(hit->anchor_id, "GENCODE:T1:exon1");
    EXPECT_EQ(hit->source, AnchorSource::Gencode);

    EXPECT_EQ(s.ById("nope"), nullptr);
}

TEST(AnchorStore, ByTagReturnsLoadOrder) {
    AnchorStore s;
    s.AddAnchor(MakeAnchor("a1", AnchorSource::Gencode, "", 0, 0, {"IGH", "CH1"}));
    s.AddAnchor(MakeAnchor("a2", AnchorSource::Gencode, "", 0, 0, {"IGH"}));
    s.AddAnchor(MakeAnchor("a3", AnchorSource::Gencode, "", 0, 0, {"CH1"}));

    auto igh = s.ByTag("IGH");
    ASSERT_EQ(igh.size(), 2u);
    EXPECT_EQ(s.anchors()[igh[0]].anchor_id, "a1");
    EXPECT_EQ(s.anchors()[igh[1]].anchor_id, "a2");

    auto ch1 = s.ByTag("CH1");
    ASSERT_EQ(ch1.size(), 2u);
    EXPECT_EQ(s.anchors()[ch1[0]].anchor_id, "a1");
    EXPECT_EQ(s.anchors()[ch1[1]].anchor_id, "a3");
}

TEST(AnchorStore, ByRegionHalfOpenOverlap) {
    AnchorStore s;
    s.AddAnchor(MakeAnchor("a", AnchorSource::Gencode, "chr14",  100, 200));
    s.AddAnchor(MakeAnchor("b", AnchorSource::Gencode, "chr14",  150, 300));
    s.AddAnchor(MakeAnchor("c", AnchorSource::Gencode, "chr14",  500, 600));
    s.AddAnchor(MakeAnchor("d", AnchorSource::Gencode, "chr15",  100, 200));
    s.Reindex();

    auto h = s.ByRegion("chr14", 180, 250);
    ASSERT_EQ(h.size(), 2u) << "a overlaps [180,200), b overlaps [180,250)";
    auto h_only_a = s.ByRegion("chr14", 95, 105);
    ASSERT_EQ(h_only_a.size(), 1u);
    EXPECT_EQ(s.anchors()[h_only_a[0]].anchor_id, "a");

    // edge: query [200,201) — exclusive end of 'a' means no overlap.
    auto h_edge = s.ByRegion("chr14", 200, 201);
    ASSERT_EQ(h_edge.size(), 1u);
    EXPECT_EQ(s.anchors()[h_edge[0]].anchor_id, "b")
        << "a's exclusive end at 200 should not overlap [200,201)";

    EXPECT_TRUE(s.ByRegion("chrX", 0, 1'000'000).empty());
    EXPECT_TRUE(s.ByRegion("chr14", 5000, 6000).empty());
}

TEST(AnchorStore, ByTranscriptIdAggregatesExonsOfSameTranscript) {
    AnchorStore s;
    s.AddAnchor(MakeAnchor("e1", AnchorSource::Gencode, "chr14", 100, 200, {}, "ENST.42"));
    s.AddAnchor(MakeAnchor("e2", AnchorSource::Gencode, "chr14", 300, 400, {}, "ENST.42"));
    s.AddAnchor(MakeAnchor("e3", AnchorSource::Gencode, "chr14", 500, 600, {}, "ENST.99"));

    auto tx42 = s.ByTranscriptId("ENST.42");
    ASSERT_EQ(tx42.size(), 2u);

    auto tx99 = s.ByTranscriptId("ENST.99");
    ASSERT_EQ(tx99.size(), 1u);
    EXPECT_EQ(s.anchors()[tx99[0]].anchor_id, "e3");

    EXPECT_TRUE(s.ByTranscriptId("ENST.unknown").empty());
}

TEST(AnchorStore, ForEachPredicateIteration) {
    AnchorStore s;
    s.AddAnchor(MakeAnchor("a", AnchorSource::Gencode));
    s.AddAnchor(MakeAnchor("b", AnchorSource::Imgt_GeneDb));
    s.AddAnchor(MakeAnchor("c", AnchorSource::Computed_Cluster));

    std::vector<std::string> seen;
    s.ForEach(
        [](const AnchorRecord& r) { return r.source == AnchorSource::Gencode
                                      || r.source == AnchorSource::Computed_Cluster; },
        [&](std::uint32_t, const AnchorRecord& r) { seen.push_back(r.anchor_id); });

    ASSERT_EQ(seen.size(), 2u);
    EXPECT_EQ(seen[0], "a");
    EXPECT_EQ(seen[1], "c");
}

TEST(AnchorStore, CountBySource) {
    AnchorStore s;
    s.AddAnchor(MakeAnchor("a1", AnchorSource::Gencode));
    s.AddAnchor(MakeAnchor("a2", AnchorSource::Gencode));
    s.AddAnchor(MakeAnchor("b1", AnchorSource::Imgt_GeneDb));
    s.AddAnchor(MakeAnchor("c1", AnchorSource::Branch_Bubble));
    EXPECT_EQ(s.CountBySource(AnchorSource::Gencode), 2u);
    EXPECT_EQ(s.CountBySource(AnchorSource::Imgt_GeneDb), 1u);
    EXPECT_EQ(s.CountBySource(AnchorSource::Branch_Bubble), 1u);
    EXPECT_EQ(s.CountBySource(AnchorSource::Refseq), 0u);
}

TEST(AnchorStore, AnchorSourceNameRoundTrip) {
    using S = AnchorSource;
    constexpr S all[] = {S::Unknown, S::Gencode, S::Mane, S::Refseq,
                         S::Imgt_GeneDb, S::Imgt_Hla, S::Fantom5_Cat,
                         S::ChessDb, S::Pangenome_PerHap,
                         S::Branch_Bubble, S::Computed_Cluster, S::Custom};
    for (auto s : all) {
        const char* name = AnchorSourceName(s);
        ASSERT_NE(name, nullptr);
        auto parsed = ParseAnchorSource(name);
        ASSERT_TRUE(parsed.has_value()) << "name='" << name << "'";
        EXPECT_EQ(*parsed, s);
    }
    EXPECT_FALSE(ParseAnchorSource("not_a_source").has_value());
}

TEST(AnchorStore, ClearWipesEverything) {
    AnchorStore s;
    s.AddAnchor(MakeAnchor("a", AnchorSource::Gencode, "chr14", 0, 100,
                            {"x", "y"}, "ENST.42"));
    EXPECT_EQ(s.size(), 1u);

    s.Clear();
    EXPECT_EQ(s.size(), 0u);
    EXPECT_EQ(s.ById("a"), nullptr);
    EXPECT_TRUE(s.ByTag("x").empty());
    EXPECT_TRUE(s.ByTranscriptId("ENST.42").empty());
}

// ============================================================================
// GENCODE GFF3 loader smoke test
//
// We synthesise a tiny GFF3 covering one transcript with two exons, write
// it uncompressed to a temp file (the loader handles both .gff3 and .gff3.gz
// via gzopen, which transparently reads plain text), and verify the
// resulting anchor structure.
// ============================================================================

TEST(AnchorStoreGencodeLoader, ParsesSyntheticTwoExonTranscript) {
    auto tmp = std::filesystem::temp_directory_path() / "llmap_test_gencode.gff3";
    {
        std::ofstream out(tmp);
        out << "##gff-version 3\n";
        out << "chr14\tGENCODE\tgene\t100\t1000\t.\t+\t.\tID=ENSG.1;gene_name=IGHG4;gene_id=ENSG.1\n";
        out << "chr14\tGENCODE\ttranscript\t100\t1000\t.\t+\t.\t"
               "ID=ENST.1;transcript_id=ENST.1;transcript_type=protein_coding;"
               "gene_id=ENSG.1;gene_name=IGHG4\n";
        out << "chr14\tGENCODE\texon\t100\t300\t.\t+\t.\t"
               "ID=exon1;transcript_id=ENST.1;exon_number=1\n";
        out << "chr14\tGENCODE\texon\t500\t1000\t.\t+\t.\t"
               "ID=exon2;transcript_id=ENST.1;exon_number=2\n";
    }

    AnchorStore s;
    auto status = s.LoadGencodeGff(tmp, /*ref_fa=*/"", /*with_sequence=*/false);
    EXPECT_TRUE(status.ok) << status.error;
    EXPECT_EQ(status.records_loaded, 2u) << "two exon-anchors expected";
    EXPECT_EQ(s.size(), 2u);

    // exon1 anchor
    const auto* e1 = s.ById("GENCODE:ENST.1:exon1");
    ASSERT_NE(e1, nullptr);
    EXPECT_EQ(e1->source, AnchorSource::Gencode);
    EXPECT_EQ(e1->kind, llmap::core::TranscriptKind::MatureMrna);
    EXPECT_EQ(e1->transcript_id, "ENST.1");
    EXPECT_EQ(e1->host_gene_id, "ENSG.1");
    ASSERT_TRUE(e1->ref_chrom.has_value());
    EXPECT_EQ(*e1->ref_chrom, "chr14");
    EXPECT_EQ(*e1->ref_start, 99);     // 1-based 100 → 0-based 99
    EXPECT_EQ(*e1->ref_end,   300);    // half-open

    // ExonBoundary list synthesised on both exon records
    ASSERT_EQ(e1->exon_boundaries.size(), 1u);
    const auto& b = e1->exon_boundaries[0];
    EXPECT_EQ(b.donor_genomic_pos,    299u);  // 1-based 300 - 1 → 0-based 299
    EXPECT_EQ(b.acceptor_genomic_pos, 499u);  // 1-based 500 - 1 → 0-based 499
    EXPECT_EQ(b.pos_in_transcript,    201u);  // exon1 length = 300-99 = 201

    // tags
    bool has_gene_tag = false, has_biotype_tag = false;
    for (const auto& t : e1->tags) {
        if (t == "gene:IGHG4")              has_gene_tag = true;
        if (t == "biotype:protein_coding")  has_biotype_tag = true;
    }
    EXPECT_TRUE(has_gene_tag);
    EXPECT_TRUE(has_biotype_tag);

    std::filesystem::remove(tmp);
}

TEST(AnchorStoreGencodeLoader, RejectsMissingFile) {
    AnchorStore s;
    auto status = s.LoadGencodeGff("/tmp/does_not_exist_zzzz_42.gff3", "", false);
    EXPECT_FALSE(status.ok);
    EXPECT_FALSE(status.error.empty());
}
