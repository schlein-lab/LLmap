// LLmap — Tests for IMGT FASTA, MANE TSV, BRANCH bubble BED loaders.

#include "anchor/anchor_store.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

using namespace llmap::anchor;

namespace {

std::filesystem::path WriteTmp(const std::string& name,
                                const std::string& body) {
    auto p = std::filesystem::temp_directory_path() / name;
    std::ofstream out(p);
    out << body;
    return p;
}

}  // namespace

// ===========================================================================
// LoadImgtGeneDb
// ===========================================================================

TEST(ImgtGeneDb, LoadsTwoRecordsFromMinimalFasta) {
    auto p = WriteTmp("llmap_test_imgt.fa",
        ">IGHV1-2*01\n"
        "CAGGTGCAGCTGGTGCAGTCTGGGGCT\n"
        ">IGHC*01\n"
        "GCCTCCACCAAGGGCCCATCCGTCTTC\n");
    AnchorStore s;
    auto status = s.LoadImgtGeneDb(p);
    EXPECT_TRUE(status.ok) << status.error;
    EXPECT_EQ(status.records_loaded, 2u);
    EXPECT_EQ(s.size(), 2u);

    const auto* v = s.ById("IMGT:IGHV1-2*01");
    ASSERT_NE(v, nullptr);
    EXPECT_EQ(v->source, AnchorSource::Imgt_GeneDb);
    EXPECT_FALSE(v->sequence.empty());

    bool has_igh_tag = false, has_v_tag = false;
    for (const auto& t : v->tags) {
        if (t == "IGH")    has_igh_tag = true;
        if (t == "V_gene") has_v_tag = true;
    }
    EXPECT_TRUE(has_igh_tag);
    EXPECT_TRUE(has_v_tag);

    std::filesystem::remove(p);
}

TEST(ImgtGeneDb, MissingFileFails) {
    AnchorStore s;
    auto status = s.LoadImgtGeneDb("/tmp/does_not_exist_xy_42.fa");
    EXPECT_FALSE(status.ok);
}

// ===========================================================================
// ImportBranchBubbles
// ===========================================================================

TEST(BranchBubbles, ImportsBedRowsAsAnchors) {
    auto p = WriteTmp("llmap_test_branch.bed",
        "chr14\t105000000\t105005000\tbubble_42\t0.95\t10,5,2\t0.02\t17\n"
        "chr14\t105100000\t105105000\tbubble_99\t0.85\t8\t0.10\t8\n");
    AnchorStore s;
    auto status = s.ImportBranchBubbles(p);
    EXPECT_TRUE(status.ok) << status.error;
    EXPECT_EQ(status.records_loaded, 2u);

    const auto* a = s.ById("BRANCH:bubble_42");
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->source, AnchorSource::Branch_Bubble);
    EXPECT_EQ(a->kind, llmap::core::TranscriptKind::NovelUnclassified);
    bool has_branch_tag = false, has_vaf_tag = false;
    for (const auto& t : a->tags) {
        if (t == "branch_bubble") has_branch_tag = true;
        if (t == "vaf:0.02")      has_vaf_tag = true;
    }
    EXPECT_TRUE(has_branch_tag);
    EXPECT_TRUE(has_vaf_tag);

    std::filesystem::remove(p);
}

TEST(BranchBubbles, MissingBedFails) {
    AnchorStore s;
    auto status = s.ImportBranchBubbles("/tmp/does_not_exist_xy_42.bed");
    EXPECT_FALSE(status.ok);
}

// ===========================================================================
// LoadMane
// ===========================================================================

TEST(ManeSelect, TagsExistingGencodeAnchors) {
    // Pre-populate one GENCODE anchor.
    AnchorRecord pre;
    pre.anchor_id = "GENCODE:ENST00000390557.4:exon1";
    pre.source = AnchorSource::Gencode;
    pre.transcript_id = "ENST00000390557";
    AnchorStore s;
    s.AddAnchor(std::move(pre));
    s.Reindex();

    auto p = WriteTmp("llmap_test_mane.tsv",
        "GENCODE\tRefSeq\tGene\n"
        "ENST00000390557\tNM_000533\tIGHG4\n");
    auto status = s.LoadMane(p);
    EXPECT_TRUE(status.ok);
    EXPECT_EQ(status.records_loaded, 1u);

    const auto* a = s.ById("GENCODE:ENST00000390557.4:exon1");
    ASSERT_NE(a, nullptr);
    bool has_mane = false;
    for (const auto& t : a->tags) {
        if (t == "mane_select") has_mane = true;
    }
    EXPECT_TRUE(has_mane);
    std::filesystem::remove(p);
}

TEST(ManeSelect, MissingFileFails) {
    AnchorStore s;
    auto status = s.LoadMane("/tmp/does_not_exist_xy_42.tsv");
    EXPECT_FALSE(status.ok);
}

// ===========================================================================
// TranscriptMode round-trip
// ===========================================================================

#include "core/transcript_mode.h"
using llmap::core::TranscriptMode;
using llmap::core::TranscriptModeName;
using llmap::core::ParseTranscriptMode;

TEST(TranscriptMode, NameRoundTrip) {
    constexpr TranscriptMode all[] = {
        TranscriptMode::Auto,
        TranscriptMode::Transcript,
        TranscriptMode::GenomeReads,
        TranscriptMode::Assembly,
        TranscriptMode::ReadsVsAssembly,
    };
    for (auto m : all) {
        const char* name = TranscriptModeName(m);
        ASSERT_NE(name, nullptr);
        auto parsed = ParseTranscriptMode(name);
        ASSERT_TRUE(parsed.has_value()) << "name='" << name << "'";
        EXPECT_EQ(*parsed, m);
    }
    EXPECT_FALSE(ParseTranscriptMode("nonsense").has_value());
    // 'genome_reads' should map to GenomeReads alongside 'reads'.
    EXPECT_EQ(ParseTranscriptMode("genome_reads"), TranscriptMode::GenomeReads);
}
