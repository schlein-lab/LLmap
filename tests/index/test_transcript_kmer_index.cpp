// LLmap — TranscriptKmerIndex tests.
//
// Cover all six origin tables + the lookup-direction asymmetry that
// motivates the segregation:
//   - IntraExon — bidirectional (matches genome AND reads)
//   - JunctionSpanning — read-only (must NOT match genome)
//   - BackSpliceSpanning — read-only, circRNA reads
//   - SterileIntronic — read-only, sterile germline transcripts
//   - ShortRna — k=21 fallback for miRNA/piRNA/siRNA
//   - PreMrnaIntronic — opt-in, scaffold-only

#include "index/transcript_kmer_index.h"

#include "anchor/anchor_store.h"
#include "annot/splice_site_db.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

using namespace llmap::index;
using llmap::anchor::AnchorRecord;
using llmap::anchor::AnchorSource;
using llmap::anchor::AnchorStore;
using llmap::anchor::ExonBoundary;
using llmap::annot::SpliceSiteDb;
using llmap::core::TranscriptKind;

namespace {

AnchorRecord MakeExonAnchor(std::string id,
                             std::string seq,
                             std::string chrom = "chr14",
                             std::int64_t start = 100,
                             std::vector<ExonBoundary> boundaries = {}) {
    AnchorRecord r;
    r.anchor_id = std::move(id);
    r.source = AnchorSource::Gencode;
    r.kind = TranscriptKind::MatureMrna;
    r.sequence = std::move(seq);
    r.ref_chrom = std::move(chrom);
    r.ref_start = start;
    r.ref_end = start + static_cast<std::int64_t>(r.sequence.size());
    r.strand = '+';
    r.exon_boundaries = std::move(boundaries);
    return r;
}

SpliceSiteDb MakeSplice() {
    SpliceSiteDb d;
    d.LoadDefaults();
    return d;
}

TranscriptKmerIndexConfig SmallCfg() {
    TranscriptKmerIndexConfig c;
    c.k_intra = 11;       // tiny for synthetic tests
    c.k_junction = 8;
    c.k_circular = 8;
    c.alt_k = 8;
    return c;
}

}  // namespace

TEST(TranscriptKmerIndex, EmptyIndexHasZeroTotals) {
    TranscriptKmerIndex idx;
    EXPECT_EQ(idx.TotalKmers(), 0u);
    EXPECT_EQ(idx.TableSize(KmerOrigin::IntraExon), 0u);
    EXPECT_TRUE(idx.QueryReadKmer("ACGTACGTACGTAC").empty());
    EXPECT_TRUE(idx.QueryGenomeWindow("ACGTACGTACGT").empty());
}

TEST(TranscriptKmerIndex, IntraExonKmerFoundInReadAndGenome) {
    AnchorStore store;
    // 30-bp anchor sequence (long enough for k_intra=11 → many k-mers)
    store.AddAnchor(MakeExonAnchor("a1",
        "ACGTACGTACGTACGTACGTACGTACGTAC"));   // 30 bp
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());

    // Pick a known 11-mer inside the anchor: positions 0..10 = ACGTACGTACG
    auto hits = idx.QueryReadKmer("ACGTACGTACG");
    ASSERT_FALSE(hits.empty());
    EXPECT_EQ(hits[0].anchor_id_idx, 0u);
    EXPECT_EQ(hits[0].origin, KmerOrigin::IntraExon);

    // The same k-mer must ALSO appear when queried via genome window.
    auto gh = idx.QueryGenomeWindow("ACGTACGTACGTACGTACGTACGT");
    ASSERT_FALSE(gh.empty());
}

TEST(TranscriptKmerIndex, JunctionSpanningKmerOnlyMatchesReadNotGenome) {
    // Construct an anchor with two exons concatenated.
    // exon1 length = 12 bp, exon2 length = 12 bp.
    // junction at pos_in_transcript = 12.
    // A k=8 junction-spanning window centred on pos 12 covers positions 8..15
    // i.e. the 8-mer ACGTACGT|GTCAGTCA crossing the boundary.
    std::string seq = "ACGTACGTACGT" "GTCAGTCAGTCA";  // 12+12=24 bp
    ExonBoundary b;
    b.pos_in_transcript = 12;
    b.spliceosome_class = 0;
    auto anchor = MakeExonAnchor("jx", seq, "chr14", 100, {b});

    AnchorStore store;
    store.AddAnchor(std::move(anchor));
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());

    // Pick a junction-spanning 8-mer that straddles position 12.
    // bases 8..15: "ACGT" + "GTCA" = "ACGTGTCA"
    auto read_hits = idx.QueryReadKmer("ACGTGTCA");
    bool seen_junction = false;
    for (const auto& h : read_hits) {
        if (h.origin == KmerOrigin::JunctionSpanning) seen_junction = true;
    }
    EXPECT_TRUE(seen_junction)
        << "junction-spanning k-mer must appear in QueryReadKmer";

    // Now query via genome window — junction-spanning k-mers MUST NOT
    // be in the genome-direction result. We pass a genomic-style window
    // that just happens to contain "ACGTGTCA"; only IntraExon hits
    // should come back.
    auto gh = idx.QueryGenomeWindow("XXXACGTGTCAYYY");  // X/Y filtered out
    for (const auto& h : gh) {
        EXPECT_NE(h.origin, KmerOrigin::JunctionSpanning)
            << "genome-direction queries must never return JunctionSpanning";
    }
}

TEST(TranscriptKmerIndex, BackSpliceSpanningOnlyForCircularAnchor) {
    // Build a circular-flagged anchor (spliceosome_class == 3) of 20 bp.
    ExonBoundary b;
    b.pos_in_transcript = 10;
    b.spliceosome_class = 3;  // back-splice
    auto anchor = MakeExonAnchor("circ", "AAAACCCCGGGGTTTTAAAA", "chr14",
                                  100, {b});

    AnchorStore store;
    store.AddAnchor(std::move(anchor));
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());

    // BackSpliceSpanning table must be populated.
    EXPECT_GT(idx.TableSize(KmerOrigin::BackSpliceSpanning), 0u);
}

TEST(TranscriptKmerIndex, SterileIntronicForSwitchRegionAnchor) {
    AnchorRecord a = MakeExonAnchor("sg", "ACGTACGTACGTACGTACGT");
    a.tags = {"switch_region_Sgamma4"};
    a.host_gene_id = "IGHG4";
    a.kind = TranscriptKind::SterileGermline;
    AnchorStore store;
    store.AddAnchor(std::move(a));
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());

    // SterileIntronic table populated for switch-region anchor.
    EXPECT_GT(idx.TableSize(KmerOrigin::SterileIntronic), 0u);
}

TEST(TranscriptKmerIndex, ShortRnaTableForSmallRnaAnchor) {
    AnchorRecord mir = MakeExonAnchor("mir1", "ACGTACGTACGTACGTACGT");  // 20 bp
    mir.kind = TranscriptKind::Mirna;
    AnchorStore store;
    store.AddAnchor(std::move(mir));
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());

    EXPECT_GT(idx.TableSize(KmerOrigin::ShortRna), 0u);
    EXPECT_EQ(idx.TableSize(KmerOrigin::IntraExon), 0u)
        << "small-RNA anchors should NOT populate the intra-exon table";
}

TEST(TranscriptKmerIndex, OriginNameRoundTripStable) {
    EXPECT_STREQ(KmerOriginName(KmerOrigin::IntraExon),          "intra_exon");
    EXPECT_STREQ(KmerOriginName(KmerOrigin::JunctionSpanning),   "junction_spanning");
    EXPECT_STREQ(KmerOriginName(KmerOrigin::BackSpliceSpanning), "back_splice_spanning");
    EXPECT_STREQ(KmerOriginName(KmerOrigin::SterileIntronic),    "sterile_intronic");
    EXPECT_STREQ(KmerOriginName(KmerOrigin::PreMrnaIntronic),    "pre_mrna_intronic");
    EXPECT_STREQ(KmerOriginName(KmerOrigin::ShortRna),           "short_rna");
}

TEST(TranscriptKmerIndex, RebuildClearsPreviousState) {
    AnchorStore store;
    store.AddAnchor(MakeExonAnchor("a", "ACGTACGTACGTACGT"));
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());
    EXPECT_GT(idx.TotalKmers(), 0u);

    // Now rebuild from an empty store — should clear.
    AnchorStore empty;
    idx.BuildFromAnchorStore(empty, splice, SmallCfg());
    EXPECT_EQ(idx.TotalKmers(), 0u);
}

TEST(TranscriptKmerIndex, NonAcgtBasesSkipSilently) {
    AnchorStore store;
    // Anchor sequence contains 'N' — surrounding k-mers skipped, others kept.
    store.AddAnchor(MakeExonAnchor(
        "a", "ACGTACGTNCGTACGTACGTACGT"));
    store.Reindex();

    TranscriptKmerIndex idx;
    auto splice = MakeSplice();
    idx.BuildFromAnchorStore(store, splice, SmallCfg());

    // The table should still have SOME k-mers (the windows that don't
    // overlap the 'N').
    EXPECT_GT(idx.TableSize(KmerOrigin::IntraExon), 0u);
    // ...but querying the 'N'-window itself returns nothing.
    EXPECT_TRUE(idx.QueryReadKmer("CGTACGTN").empty());
}
