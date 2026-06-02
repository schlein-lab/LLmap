// LLmap — TranscriptKindClassifier tests.
//
// Synthetic anchors covering every detector branch:
//   - mature mRNA (loader-set kind preserved)
//   - sterile germline (IGH I-promoter→S→C, no V-D-J)
//   - circRNA (back-splice exon boundary)
//   - pre-mRNA (intron-retained tag)
//   - piRNA (~28 nt)
//   - siRNA (~22 nt)
//   - snoRNA C/D box + H/ACA box (motif-driven)
//   - novel-unclassified (no signature)

#include "annot/transcript_kind_classifier.h"
#include "annot/splice_site_db.h"
#include "anchor/anchor_record.h"

#include <gtest/gtest.h>

using namespace llmap::annot;
using namespace llmap::anchor;
using llmap::core::TranscriptKind;

namespace {

SpliceSiteDb MakeSplice() {
    SpliceSiteDb d;
    d.LoadDefaults();
    return d;
}

AnchorRecord MakeAnchor(std::string id,
                         TranscriptKind kind = TranscriptKind::Unknown,
                         std::string seq = "",
                         std::vector<std::string> tags = {},
                         std::string host_gene = "") {
    AnchorRecord r;
    r.anchor_id = std::move(id);
    r.kind = kind;
    r.sequence = std::move(seq);
    r.tags = std::move(tags);
    r.host_gene_id = std::move(host_gene);
    return r;
}

ExonBoundary BackSplice() {
    ExonBoundary b;
    b.spliceosome_class = 3;
    return b;
}

}  // namespace

TEST(TranscriptKindClassifier, PreservesLoaderKindWhenNothingSpecific) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("a", TranscriptKind::MatureMrna);
    auto [k, tag] = c.Classify(a, splice, /*tissue*/ {});
    EXPECT_EQ(k, TranscriptKind::MatureMrna);
    EXPECT_FALSE(tag.has_value());
}

TEST(TranscriptKindClassifier, DetectsCircularRna) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("circ");
    a.exon_boundaries.push_back(BackSplice());
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::CircularRna);
}

TEST(TranscriptKindClassifier, DetectsSterileGermlineTranscript) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("sg", TranscriptKind::MatureMrna,
                         /*seq*/ "",
                         /*tags*/ {"switch_region_Sgamma4"},
                         /*host*/ "IGHG4");
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::SterileGermline);
}

TEST(TranscriptKindClassifier, SterileGermlineRejectedWhenVdjPresent) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("not_sg", TranscriptKind::MatureMrna, "",
                         {"switch_region_Sgamma4", "V_gene_IGHV3-23"},
                         "IGHG4");
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::MatureMrna)  // loader kind preserved
        << "presence of V_gene tag must veto sterile classification";
}

TEST(TranscriptKindClassifier, DetectsPreMrna) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("pre", TranscriptKind::MatureMrna, "",
                         {"intron_retained"});
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::PreMrna);
}

TEST(TranscriptKindClassifier, DetectsPirnaByLength) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("pir");
    a.sequence = std::string(28, 'A');  // 28 nt → piRNA range
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::Pirna);
}

TEST(TranscriptKindClassifier, DetectsSirnaByLength) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("si");
    a.sequence = std::string(22, 'A');  // 22 nt → siRNA range
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::Sirna);
}

TEST(TranscriptKindClassifier, RefinesCDBoxSnornaByMotif) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    // 5' AUGAUGA + 3' CUGA + body length 80
    std::string seq = "GTGATGA" + std::string(70, 'A') + "CTGA";
    auto a = MakeAnchor("sno", TranscriptKind::Snorna_CDbox, seq);
    auto k = c.RefineSnornaSubclass(a);
    EXPECT_EQ(k, TranscriptKind::Snorna_CDbox);
}

TEST(TranscriptKindClassifier, RefinesHACABoxSnornaByMotif) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    // 5' ACANNA + 3' ACA terminus, length 140
    std::string seq = std::string(50, 'G') + "ACANNA"
                       + std::string(80, 'G') + "ACA";
    auto a = MakeAnchor("snoH", TranscriptKind::Snorna_HacaBox, seq);
    auto k = c.RefineSnornaSubclass(a);
    EXPECT_EQ(k, TranscriptKind::Snorna_HacaBox);
}

TEST(TranscriptKindClassifier, OpenEndedAssignsCustomKindTag) {
    TranscriptKindClassifier c;
    auto splice = MakeSplice();
    auto a = MakeAnchor("x");  // empty
    auto [k, tag] = c.Classify(a, splice, {});
    EXPECT_EQ(k, TranscriptKind::NovelUnclassified);
    ASSERT_TRUE(tag.has_value());
    EXPECT_TRUE(tag->label.starts_with("novel_"))
        << "label='" << tag->label << "'";
    EXPECT_FALSE(tag->reason_signature.empty());
}
