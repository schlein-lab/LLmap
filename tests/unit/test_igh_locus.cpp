// LLmap — unit tests for the IGH locus post-hoc re-sort module.

#include "igh_locus/igh_anchor_catalog.h"
#include "igh_locus/igh_region.h"
#include "igh_locus/igh_resort.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

using namespace llmap;
using namespace llmap::igh_locus;

namespace {

// Distinct synthetic exon sequences (no cross-substring collisions).
const std::string kG4_CH1 = "ACACACGTGTACACGTACGTACACGTGTACACGTACGTACACGTGTACAC";
const std::string kG4_CH2 = "TTGGTTGGCCAACCAATTGGCCAATTGGCCAATTGGTTGGCCAACCAATT";
const std::string kG2_CH1 = "GAGAGTCTCTGAGTCTCAGAGTCTCTGAGTCTCAGAGTCTCTGAGTCTCA";
const std::string kG2_CH2 = "CCTTCCTTAAGGAAGGCCTTAAGGCCTTAAGGCCTTCCTTAAGGAAGGCC";

std::string RevComp(const std::string& s) {
    std::string out(s.size(), 'N');
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        out[i] = (c == 'A') ? 'T' : (c == 'C') ? 'G' : (c == 'G') ? 'C'
               : (c == 'T') ? 'A' : 'N';
    }
    return out;
}

std::string WriteCatalogFasta() {
    auto path = std::filesystem::temp_directory_path() /
                "llmap_igh_test_anchors.fa";
    std::ofstream out(path);
    out << ">IGHG4_C1_hap1_CH1_50bp\n" << kG4_CH1 << "\n";
    out << ">IGHG4_C1_hap1_CH2_50bp\n" << kG4_CH2 << "\n";
    out << ">IGHG2_C2_hap1_CH1_50bp\n" << kG2_CH1 << "\n";
    out << ">IGHG2_C2_hap1_CH2_50bp\n" << kG2_CH2 << "\n";
    return path.string();
}

// Two IGHG4 copies sharing an identical CH2 but with distinct CH1: a read over
// only the shared CH2 is copy-ambiguous; a read over CH1 resolves to one copy.
const std::string kC1_CH1 = "ACACACGTGTACACGTACGTACACGTGTACACGTACGTACACGTGTACAC";
const std::string kC2_CH1 = "GAGAGTCTCTGAGTCTCAGAGTCTCTGAGTCTCAGAGTCTCTGAGTCTCA";
const std::string kShared_CH2 = "TTGGTTGGCCAACCAATTGGCCAATTGGCCAATTGGTTGGCCAACCAATT";

std::string WriteAmbigCatalogFasta() {
    auto path = std::filesystem::temp_directory_path() /
                "llmap_igh_ambig_anchors.fa";
    std::ofstream out(path);
    out << ">IGHG4_C1_hap1_CH1_50bp\n" << kC1_CH1 << "\n";
    out << ">IGHG4_C1_hap1_CH2_50bp\n" << kShared_CH2 << "\n";
    out << ">IGHG4_C2_hap2_CH1_50bp\n" << kC2_CH1 << "\n";
    out << ">IGHG4_C2_hap2_CH2_50bp\n" << kShared_CH2 << "\n";  // identical CH2
    return path.string();
}

}  // namespace

TEST(IghAnchorCatalog, AmbiguousCopyOnSharedExon) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteAmbigCatalogFasta());
    ASSERT_TRUE(cat.has_value());
    // Read covers ONLY the shared CH2 -> both copies tie at 1 exon.
    IghMatch m = cat->Match("GG" + kShared_CH2 + "CC");
    EXPECT_TRUE(m.matched);
    EXPECT_EQ(m.gene, "IGHG4");
    EXPECT_TRUE(m.ambiguous_copy);
    EXPECT_EQ(m.tied_copies.size(), 2u);
}

TEST(IghAnchorCatalog, ResolvedCopyOnDiscriminatingExon) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteAmbigCatalogFasta());
    ASSERT_TRUE(cat.has_value());
    // Read covers C1's distinct CH1 + the shared CH2 -> C1 wins (2 vs 1 exon).
    IghMatch m = cat->Match("GG" + kC1_CH1 + "AT" + kShared_CH2 + "CC");
    EXPECT_TRUE(m.matched);
    EXPECT_EQ(m.copy_label, "IGHG4_C1_hap1");
    EXPECT_FALSE(m.ambiguous_copy);
    EXPECT_EQ(m.n_distinct_exons, 2);
}

TEST(IghRegion, GenomicMembershipSuffixContig) {
    IghRegion r = IghRegion::Default();
    EXPECT_TRUE(r.Contains("chr14", 105'600'000ULL));
    EXPECT_TRUE(r.Contains("HG002#1#chr14", 105'600'000ULL));
    EXPECT_FALSE(r.Contains("chr14", 1'000'000ULL));
    EXPECT_FALSE(r.Contains("chr7", 105'600'000ULL));
}

TEST(IghAnchorCatalog, LoadAndGeneOfTarget) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());
    EXPECT_EQ(cat->size(), 4u);
    EXPECT_EQ(cat->GeneOfTarget("IGHG4_C1_hap1").value_or(""), "IGHG4");
    EXPECT_EQ(cat->GeneOfTarget("IGHG2").value_or(""), "IGHG2");
    EXPECT_FALSE(cat->GeneOfTarget("chr7").has_value());
}

TEST(IghAnchorCatalog, MatchForwardTwoExons) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());
    std::string read = "GGGGG" + kG4_CH1 + "ATATAT" + kG4_CH2 + "CCCCC";
    IghMatch m = cat->Match(read);
    EXPECT_TRUE(m.matched);
    EXPECT_EQ(m.gene, "IGHG4");
    EXPECT_EQ(m.copy_label, "IGHG4_C1_hap1");
    EXPECT_EQ(m.n_distinct_exons, 2);
    EXPECT_FALSE(m.is_reverse);
    EXPECT_FALSE(m.ambiguous_gene);
}

TEST(IghAnchorCatalog, MatchReverseStrand) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());
    std::string read = "TT" + RevComp(kG2_CH2) + "AA" + RevComp(kG2_CH1) + "GG";
    IghMatch m = cat->Match(read);
    EXPECT_TRUE(m.matched);
    EXPECT_EQ(m.gene, "IGHG2");
    EXPECT_TRUE(m.is_reverse);
    EXPECT_EQ(m.n_distinct_exons, 2);
}

TEST(IghAnchorCatalog, NoMatchForUnrelatedRead) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());
    IghMatch m = cat->Match("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    EXPECT_FALSE(m.matched);
}

TEST(IghAnchorCatalog, MaxMismatchToleratesError) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta(), /*max_mismatch=*/2);
    ASSERT_TRUE(cat.has_value());
    // Introduce 1 substitution into each G4 exon.
    std::string e1 = kG4_CH1; e1[25] = (e1[25] == 'A') ? 'C' : 'A';
    std::string e2 = kG4_CH2; e2[10] = (e2[10] == 'T') ? 'G' : 'T';
    std::string read = "GG" + e1 + "ATAT" + e2 + "CC";
    IghMatch m = cat->Match(read);
    EXPECT_TRUE(m.matched);
    EXPECT_EQ(m.gene, "IGHG4");
    EXPECT_EQ(m.n_distinct_exons, 2);
}

TEST(IghResort, RelabelsTranscriptomicMisassignment) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());

    // Read truly from IGHG4 but mapper put it on the IGHG2 copy.
    AlignmentHit hit;
    hit.target_id = "IGHG2_C2_hap1";
    hit.start = 0;
    hit.end = 100;
    AlignmentRecord rec = make_mapped("read1", 120, hit);

    std::vector<AlignmentRecord> records{rec};
    std::vector<std::string> seqs{"GG" + kG4_CH1 + "AT" + kG4_CH2 + "CC"};

    ResortStats st = ApplyResort(*cat, records, seqs);
    EXPECT_EQ(st.n_in_locus, 1u);
    EXPECT_EQ(st.n_matched, 1u);
    EXPECT_EQ(st.n_resorted, 1u);
    ASSERT_TRUE(records[0].primary.has_value());
    EXPECT_EQ(records[0].primary->target_id, "IGHG4_C1_hap1");
    ASSERT_TRUE(records[0].paralog_assignment.has_value());
    ASSERT_FALSE(records[0].paralog_assignment->inter_paralog.empty());
    EXPECT_EQ(records[0].paralog_assignment->inter_paralog.front().first, "IGHG4");
    EXPECT_TRUE(records[0].is_lossless_consistent());
}

TEST(IghResort, LeavesCorrectAssignmentButSetsParalog) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());

    AlignmentHit hit;
    hit.target_id = "IGHG4_C1_hap1";
    AlignmentRecord rec = make_mapped("read2", 120, hit);
    std::vector<AlignmentRecord> records{rec};
    std::vector<std::string> seqs{"GG" + kG4_CH1 + "AT" + kG4_CH2 + "CC"};

    ResortStats st = ApplyResort(*cat, records, seqs);
    EXPECT_EQ(st.n_resorted, 0u);          // already correct
    EXPECT_EQ(st.n_paralog_set, 1u);
    EXPECT_EQ(records[0].primary->target_id, "IGHG4_C1_hap1");
}

TEST(IghResort, IgnoresNonIghReads) {
    auto cat = IghAnchorCatalog::LoadFasta(WriteCatalogFasta());
    ASSERT_TRUE(cat.has_value());

    AlignmentHit hit;
    hit.target_id = "chr7";
    hit.start = 1'000'000;
    AlignmentRecord rec = make_mapped("read3", 120, hit);
    std::vector<AlignmentRecord> records{rec};
    std::vector<std::string> seqs{"GG" + kG4_CH1 + "AT" + kG4_CH2 + "CC"};

    ResortStats st = ApplyResort(*cat, records, seqs);
    EXPECT_EQ(st.n_in_locus, 0u);
    EXPECT_EQ(st.n_resorted, 0u);
    EXPECT_EQ(records[0].primary->target_id, "chr7");
}
