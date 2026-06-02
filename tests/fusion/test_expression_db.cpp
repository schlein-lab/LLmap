// LLmap — ExpressionDb tests.

#include "fusion/expression_db.h"
#include "fusion/likelihood_factors.h"
#include "anchor/anchor_record.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

using namespace llmap::fusion;

namespace {

std::filesystem::path WriteTmp(const std::string& name,
                                const std::string& body) {
    auto p = std::filesystem::temp_directory_path() / name;
    std::ofstream out(p);
    out << body;
    return p;
}

}  // namespace

TEST(ExpressionDb, EmptyLookupReturnsMinusOne) {
    ExpressionDb d;
    EXPECT_LT(d.ExpectedTpm("ENST.X", "lymph"), 0.0f);
    EXPECT_LT(d.ExpectedCellTypeFraction("ENST.X", "lymph", "b_cell"), 0.0f);
    EXPECT_EQ(d.TotalEntries(), 0u);
}

TEST(ExpressionDb, LoadGtexBulkGct) {
    auto p = WriteTmp("llmap_expr_gct.gct",
        "#1.2\n"
        "2\t3\n"
        "Name\tDescription\tBlood\tSpleen\tLymph\n"
        "ENST.A\tIGHG4\t10.0\t25.0\t100.0\n"
        "ENST.B\tIGHG1\t5.0\t8.0\t30.0\n");
    ExpressionDb d;
    ASSERT_TRUE(d.LoadGtexBulkGct(p));
    EXPECT_NEAR(d.ExpectedTpm("ENST.A", "Blood"), 10.0f, 0.001f);
    EXPECT_NEAR(d.ExpectedTpm("ENST.A", "Lymph"), 100.0f, 0.001f);
    EXPECT_NEAR(d.ExpectedTpm("ENST.B", "Spleen"), 8.0f, 0.001f);
    EXPECT_LT(d.ExpectedTpm("ENST.A", "Heart"), 0.0f);
    EXPECT_EQ(d.DistinctTranscripts(), 2u);
    EXPECT_EQ(d.DistinctTissues(), 3u);
    std::filesystem::remove(p);
}

TEST(ExpressionDb, LoadTabulaSapiensTsvFeedsCellTypeAndTissue) {
    auto p = WriteTmp("llmap_expr_tabula.tsv",
        "transcript_id\torgan\tcell_type\tmean_expression\n"
        "ENST.A\tlymph\tb_cell\t0.85\n"
        "ENST.A\tlymph\tt_cell\t0.10\n"
        "ENST.A\tblood\tneutrophil\t0.20\n");
    ExpressionDb d;
    ASSERT_TRUE(d.LoadTabulaSapiensTsv(p));
    EXPECT_NEAR(d.ExpectedCellTypeFraction("ENST.A", "lymph", "b_cell"),
                 0.85f, 0.01f);
    // Fold-in: lymph organ max across cell types = 0.85
    EXPECT_NEAR(d.ExpectedTpm("ENST.A", "lymph"), 0.85f, 0.01f);
    std::filesystem::remove(p);
}

TEST(ExpressionDb, MissingFileReturnsFalse) {
    ExpressionDb d;
    EXPECT_FALSE(d.LoadGtexBulkGct("/tmp/does_not_exist_xx_42.gct"));
    EXPECT_FALSE(d.LoadGtexLongReadTsv("/tmp/does_not_exist_xx_42.tsv"));
    EXPECT_FALSE(d.LoadTabulaSapiensTsv("/tmp/does_not_exist_xx_42.tsv"));
    EXPECT_FALSE(d.LoadHcaTsv("/tmp/does_not_exist_xx_42.tsv"));
    EXPECT_FALSE(d.LoadHpaTsv("/tmp/does_not_exist_xx_42.tsv"));
    EXPECT_FALSE(d.LoadRecount3Tsv("/tmp/does_not_exist_xx_42.tsv"));
    EXPECT_EQ(d.TotalEntries(), 0u);
}

TEST(ExpressionDb, EntropyZeroForUnknownTranscript) {
    ExpressionDb d;
    EXPECT_EQ(d.ExpressionEntropy("ENST.unknown"), 0.0f);
}

TEST(ExpressionDb, EntropyHigherForBroadlyExpressedTranscript) {
    // 3 equal-expression tissues → entropy log2(3) ≈ 1.585
    // 1 dominant tissue (90 %)   → entropy ≈ 0.469
    auto p = WriteTmp("llmap_expr_entropy.gct",
        "#1.2\n"
        "2\t3\n"
        "Name\tDescription\tT1\tT2\tT3\n"
        "ENST.BROAD\tg1\t10.0\t10.0\t10.0\n"
        "ENST.NARROW\tg2\t90.0\t5.0\t5.0\n");
    ExpressionDb d;
    ASSERT_TRUE(d.LoadGtexBulkGct(p));
    const float broad = d.ExpressionEntropy("ENST.BROAD");
    const float narrow = d.ExpressionEntropy("ENST.NARROW");
    EXPECT_GT(broad, narrow);
    EXPECT_NEAR(broad, 1.585f, 0.01f);
    std::filesystem::remove(p);
}

// ===========================================================================
// Wiring: ComputeFactorsWithExpression actually consults the DB.
// ===========================================================================

TEST(ExpressionDbWiring, ExpressionPriorChangesWithTpm) {
    // Load a tiny DB: one tx with TPM=1000 in lymph, TPM=0 elsewhere.
    auto p = WriteTmp("llmap_expr_wiring.gct",
        "#1.2\n"
        "1\t2\n"
        "Name\tDescription\tlymph\tmuscle\n"
        "ENST.HIGH\tg\t1000.0\t0.0\n");
    ExpressionDb d;
    ASSERT_TRUE(d.LoadGtexBulkGct(p));

    llmap::anchor::AnchorRecord a;
    a.anchor_id = "GENCODE:ENST.HIGH:exon1";
    a.source = llmap::anchor::AnchorSource::Gencode;
    a.kind = llmap::core::TranscriptKind::MatureMrna;
    a.sequence = std::string(500, 'A');
    a.transcript_id = "ENST.HIGH";

    ReadContext r;
    r.read_id = "r1";
    r.read_length = 500;
    r.platform = "hifi";
    r.mapq = 30;

    TissueContext t_lymph{.label = "lymph"};
    TissueContext t_muscle{.label = "muscle"};

    auto f_lymph = ComputeFactorsWithExpression(r, a, {}, t_lymph, &d);
    auto f_muscle = ComputeFactorsWithExpression(r, a, {}, t_muscle, &d);

    // Plan-Block 4.5 table: TPM=1000 → 0.95 (within ±0.02 of sigmoid
    // anchor), TPM=0 → 0.18.
    EXPECT_GT(f_lymph.L_expression_prior, f_muscle.L_expression_prior);
    EXPECT_GT(f_lymph.L_expression_prior, 0.93f);
    EXPECT_LT(f_muscle.L_expression_prior, 0.25f);
    std::filesystem::remove(p);
}

TEST(ExpressionDbWiring, NullDbFallsBackToNeutralOne) {
    llmap::anchor::AnchorRecord a;
    a.anchor_id = "X";
    a.source = llmap::anchor::AnchorSource::Gencode;
    a.kind = llmap::core::TranscriptKind::MatureMrna;
    a.sequence = std::string(100, 'A');
    a.transcript_id = "ENST.X";

    ReadContext r;
    r.read_length = 100;
    r.platform = "hifi";
    r.mapq = 30;

    TissueContext t{.label = "lymph"};
    auto f = ComputeFactorsWithExpression(r, a, {}, t, /*expr_db=*/nullptr);
    EXPECT_EQ(f.L_expression_prior, 1.0f);
    EXPECT_EQ(f.L_depth_coverage,   1.0f);
}

TEST(ExpressionDbWiring, BarcodeContextUsesCellTypeFraction) {
    auto p = WriteTmp("llmap_expr_bc.tsv",
        "transcript_id\torgan\tcell_type\tmean_expression\n"
        "ENST.B\tlymph\tb_cell\t0.95\n"
        "ENST.B\tlymph\tt_cell\t0.05\n");
    ExpressionDb d;
    ASSERT_TRUE(d.LoadTabulaSapiensTsv(p));

    llmap::anchor::AnchorRecord a;
    a.anchor_id = "X";
    a.source = llmap::anchor::AnchorSource::Gencode;
    a.kind = llmap::core::TranscriptKind::MatureMrna;
    a.sequence = std::string(100, 'A');
    a.transcript_id = "ENST.B";

    ReadContext r;
    r.read_length = 100;
    r.platform = "hifi";
    r.mapq = 30;

    TissueContext t{.label = "lymph"};

    r.cell_type = "b_cell";
    auto f_b = ComputeFactorsWithExpression(r, a, {}, t, &d);
    r.cell_type = "t_cell";
    auto f_t = ComputeFactorsWithExpression(r, a, {}, t, &d);

    EXPECT_GT(f_b.L_barcode_context, f_t.L_barcode_context);
    std::filesystem::remove(p);
}
