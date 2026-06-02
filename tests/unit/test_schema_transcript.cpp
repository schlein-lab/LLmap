// LLmap — Transcript-Mode output schema + bridges tests.

#include "output/schema_transcript.h"
#include "output/pangenome_bridge.h"
#include "output/branch_bridge.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>

using namespace llmap;
using namespace llmap::output;
using namespace llmap::output::transcript_schema;

namespace {
std::string ReadFile(const std::filesystem::path& p) {
    std::ifstream in(p);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}
}  // namespace

// ===========================================================================
// BAM tag builders
// ===========================================================================

TEST(SchemaTranscript, XsTagPassThroughOrQuestionMark) {
    EXPECT_EQ(XsTag('+'), '+');
    EXPECT_EQ(XsTag('-'), '-');
    EXPECT_EQ(XsTag('.'), '?');
    EXPECT_EQ(XsTag('X'), '?');
}

TEST(SchemaTranscript, JiTagEncodesPairs) {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> j = {
        {299, 499}, {599, 999},
    };
    EXPECT_EQ(JiTag(j), "299,499,599,999");
    EXPECT_TRUE(JiTag({}).empty());
}

TEST(SchemaTranscript, JmTagFormatsTwoDecimals) {
    std::vector<float> c = {0.95f, 0.85f};
    EXPECT_EQ(JmTag(c), "0.95,0.85");
    EXPECT_TRUE(JmTag({}).empty());
}

TEST(SchemaTranscript, XkTagMatchesStatusNames) {
    EXPECT_STREQ(XkTag(AlignmentStatus::Mapped),               "MAPPED");
    EXPECT_STREQ(XkTag(AlignmentStatus::ChimericIntraRegion),  "CHIM_INTRA");
    EXPECT_STREQ(XkTag(AlignmentStatus::DarkNovel),            "DARK");
}

TEST(SchemaTranscript, XcTagDashForZeroCluster) {
    EXPECT_EQ(XcTag(0), "-");
    EXPECT_EQ(XcTag(42), "42");
}

TEST(SchemaTranscript, XaTagCommaJoins) {
    std::vector<std::string> srcs = {"GENCODE", "MANE", "IMGT"};
    EXPECT_EQ(XaTag(srcs), "GENCODE,MANE,IMGT");
    EXPECT_TRUE(XaTag({}).empty());
}

TEST(SchemaTranscript, XqTagMinus1ForUnknown) {
    EXPECT_EQ(XqTag(-1), "-1");
    EXPECT_EQ(XqTag(30), "30");
    EXPECT_EQ(XqTag(60), "60");
}

TEST(SchemaTranscript, XmTagFormatCallList) {
    std::vector<ModCallView> c = {
        {"m6a", 100u, 0.95f},
        {"psi", 250u, 0.80f},
    };
    EXPECT_EQ(XmTag(c), "m6a:100:0.95,psi:250:0.80");
}

TEST(SchemaTranscript, XfTagPassThrough) {
    EXPECT_EQ(XfTag("canonical"), "canonical");
    EXPECT_EQ(XfTag("intron_retained"), "intron_retained");
}

// ===========================================================================
// Pangenome GAF
// ===========================================================================

TEST(PangenomeBridge, GafRowTabFormat) {
    GafRow r;
    r.query_name = "read_001";
    r.query_length = 12345;
    r.query_start = 0;
    r.query_end = 1000;
    r.strand = '+';
    r.path = ">chr14:105000000";
    r.path_length = 5000;
    r.path_start = 100;
    r.path_end = 600;
    r.residue_matches = 950;
    r.block_length = 1000;
    r.mapping_quality = 60;

    auto s = EmitGafRow(r);
    EXPECT_NE(s.find("read_001\t12345"), std::string::npos);
    EXPECT_NE(s.find(">chr14:105000000"), std::string::npos);
    EXPECT_NE(s.find("\t60"), std::string::npos);
    // 12 tab-separated columns ⇒ 11 tabs
    EXPECT_EQ(std::count(s.begin(), s.end(), '\t'), 11);
}

TEST(PangenomeBridge, WriteGafFileRoundTrip) {
    auto p = std::filesystem::temp_directory_path() / "llmap_test_pan.gaf";
    GafRow r;
    r.query_name = "rA";
    r.query_length = 500;
    r.query_end = 500;
    r.path = ">chr1:1";
    r.path_length = 500;
    r.path_end = 500;
    r.residue_matches = 500;
    r.block_length = 500;
    EXPECT_TRUE(WriteGafFile(p, std::vector{r}));
    auto body = ReadFile(p);
    EXPECT_NE(body.find("rA"), std::string::npos);
    std::filesystem::remove(p);
}

// ===========================================================================
// BRANCH GAF bridge
// ===========================================================================

TEST(BranchBridge, GafRowFormat) {
    BranchGafRow r;
    r.read_id = "rB";
    r.cluster_id = 42;
    r.haplotype = '1';
    r.path = "@bubble42";
    r.score = 950;
    r.n_supporting_reads = 1;
    auto s = EmitBranchGafRow(r);
    EXPECT_NE(s.find("rB\t42\t1\t@bubble42\t950\t1\tLLMAP"), std::string::npos);
}

TEST(BranchBridge, GafRowDashForZeroCluster) {
    BranchGafRow r;
    r.read_id = "rC";
    r.cluster_id = 0;
    r.haplotype = '.';
    r.path = "|chr14:105000000";
    auto s = EmitBranchGafRow(r);
    EXPECT_NE(s.find("rC\t-\t."), std::string::npos);
}

TEST(BranchBridge, WriteBranchGafFile) {
    auto p = std::filesystem::temp_directory_path() / "llmap_test_branch.gaf";
    BranchGafRow r;
    r.read_id = "rD";
    r.path = "@bx";
    EXPECT_TRUE(WriteBranchGafFile(p, std::vector{r}));
    auto body = ReadFile(p);
    EXPECT_NE(body.find("rD"), std::string::npos);
    std::filesystem::remove(p);
}
