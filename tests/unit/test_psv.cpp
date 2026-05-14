// LLmap — PSV (Paralog-Specific Variant) module tests.

#include "psv/psv_assigner.h"
#include "psv/psv_catalog.h"
#include "psv/psv_types.h"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

namespace llmap::psv {
namespace {

// === PSV Types tests ===

TEST(PsvTypesTest, AlleleToNumValid) {
    EXPECT_EQ(AlleleToNum('A'), 0);
    EXPECT_EQ(AlleleToNum('C'), 1);
    EXPECT_EQ(AlleleToNum('G'), 2);
    EXPECT_EQ(AlleleToNum('T'), 3);
}

TEST(PsvTypesTest, AlleleToNumLowercase) {
    EXPECT_EQ(AlleleToNum('a'), 0);
    EXPECT_EQ(AlleleToNum('c'), 1);
    EXPECT_EQ(AlleleToNum('g'), 2);
    EXPECT_EQ(AlleleToNum('t'), 3);
}

TEST(PsvTypesTest, AlleleToNumInvalid) {
    EXPECT_EQ(AlleleToNum('N'), 4);
    EXPECT_EQ(AlleleToNum('X'), 4);
    EXPECT_EQ(AlleleToNum('-'), 4);
}

TEST(PsvTypesTest, NumToAllele) {
    EXPECT_EQ(NumToAllele(0), 'A');
    EXPECT_EQ(NumToAllele(1), 'C');
    EXPECT_EQ(NumToAllele(2), 'G');
    EXPECT_EQ(NumToAllele(3), 'T');
    EXPECT_EQ(NumToAllele(4), 'N');
    EXPECT_EQ(NumToAllele(5), 'N');
}

TEST(PsvTypesTest, IsValidAllele) {
    EXPECT_TRUE(IsValidAllele('A'));
    EXPECT_TRUE(IsValidAllele('C'));
    EXPECT_TRUE(IsValidAllele('G'));
    EXPECT_TRUE(IsValidAllele('T'));
    EXPECT_FALSE(IsValidAllele('N'));
    EXPECT_FALSE(IsValidAllele('X'));
}

// === PSV Catalog tests ===

TEST(PsvCatalogTest, AddSiteAndLookupById) {
    PsvCatalog catalog;

    PsvSite site;
    site.psv_id = 42;
    site.chrom = "chr1";
    site.position = 1000;
    site.ref_allele = 'A';
    site.paralog_alleles["GENE1"] = 'A';
    site.paralog_alleles["GENE2"] = 'G';
    site.informativeness = 1.0f;

    catalog.AddSite(site);

    EXPECT_EQ(catalog.Size(), 1);
    EXPECT_FALSE(catalog.Empty());

    const auto* found = catalog.GetSiteById(42);
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->chrom, "chr1");
    EXPECT_EQ(found->position, 1000);
    EXPECT_EQ(found->ref_allele, 'A');
}

TEST(PsvCatalogTest, LookupByPosition) {
    PsvCatalog catalog;

    PsvSite site;
    site.psv_id = 1;
    site.chrom = "chr2";
    site.position = 5000;
    site.ref_allele = 'C';
    site.paralog_alleles["P1"] = 'C';
    site.paralog_alleles["P2"] = 'T';

    catalog.AddSite(site);
    catalog.BuildIndex();

    const auto* found = catalog.GetSiteAtPosition("chr2", 5000);
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->psv_id, 1);
    EXPECT_EQ(found->ref_allele, 'C');

    const auto* not_found = catalog.GetSiteAtPosition("chr2", 5001);
    EXPECT_EQ(not_found, nullptr);
}

TEST(PsvCatalogTest, GetSitesInRegion) {
    PsvCatalog catalog;

    for (int i = 0; i < 10; ++i) {
        PsvSite site;
        site.psv_id = i;
        site.chrom = "chr1";
        site.position = 100 + i * 10;
        site.ref_allele = 'A';
        site.paralog_alleles["P1"] = 'A';
        site.paralog_alleles["P2"] = 'G';
        catalog.AddSite(site);
    }

    catalog.BuildIndex();

    auto sites = catalog.GetSitesInRegion("chr1", 120, 160);
    EXPECT_EQ(sites.size(), 5);

    sites = catalog.GetSitesInRegion("chr1", 0, 50);
    EXPECT_EQ(sites.size(), 0);

    sites = catalog.GetSitesInRegion("chr2", 100, 200);
    EXPECT_EQ(sites.size(), 0);
}

TEST(PsvCatalogTest, GetParalogs) {
    PsvCatalog catalog;

    PsvSite site1;
    site1.psv_id = 1;
    site1.chrom = "chr1";
    site1.position = 100;
    site1.paralog_alleles["GENE_A"] = 'A';
    site1.paralog_alleles["GENE_B"] = 'G';
    catalog.AddSite(site1);

    PsvSite site2;
    site2.psv_id = 2;
    site2.chrom = "chr1";
    site2.position = 200;
    site2.paralog_alleles["GENE_B"] = 'T';
    site2.paralog_alleles["GENE_C"] = 'C';
    catalog.AddSite(site2);

    auto paralogs = catalog.GetParalogs();
    EXPECT_EQ(paralogs.size(), 3);
}

TEST(PsvCatalogTest, GetChromosomes) {
    PsvCatalog catalog;

    PsvSite site1;
    site1.psv_id = 1;
    site1.chrom = "chr1";
    site1.position = 100;
    site1.paralog_alleles["P1"] = 'A';
    catalog.AddSite(site1);

    PsvSite site2;
    site2.psv_id = 2;
    site2.chrom = "chr2";
    site2.position = 100;
    site2.paralog_alleles["P1"] = 'G';
    catalog.AddSite(site2);

    auto chroms = catalog.GetChromosomes();
    EXPECT_EQ(chroms.size(), 2);
}

TEST(PsvCatalogTest, Clear) {
    PsvCatalog catalog;

    PsvSite site;
    site.psv_id = 1;
    site.chrom = "chr1";
    site.position = 100;
    site.paralog_alleles["P1"] = 'A';
    catalog.AddSite(site);

    EXPECT_EQ(catalog.Size(), 1);

    catalog.Clear();
    EXPECT_TRUE(catalog.Empty());
    EXPECT_EQ(catalog.Size(), 0);
}

// === Informativeness tests ===

TEST(InformativenessTest, AllSameAllele) {
    PsvSite site;
    site.paralog_alleles["P1"] = 'A';
    site.paralog_alleles["P2"] = 'A';
    site.paralog_alleles["P3"] = 'A';

    EXPECT_FLOAT_EQ(ComputeInformativeness(site), 0.0f);
}

TEST(InformativenessTest, TwoDifferentAlleles) {
    PsvSite site;
    site.paralog_alleles["P1"] = 'A';
    site.paralog_alleles["P2"] = 'G';

    EXPECT_FLOAT_EQ(ComputeInformativeness(site), 1.0f);
}

TEST(InformativenessTest, ThreeUniqueAlleles) {
    PsvSite site;
    site.paralog_alleles["P1"] = 'A';
    site.paralog_alleles["P2"] = 'G';
    site.paralog_alleles["P3"] = 'C';

    EXPECT_FLOAT_EQ(ComputeInformativeness(site), 1.0f);
}

TEST(InformativenessTest, MixedAlleles) {
    PsvSite site;
    site.paralog_alleles["P1"] = 'A';
    site.paralog_alleles["P2"] = 'A';
    site.paralog_alleles["P3"] = 'G';

    float info = ComputeInformativeness(site);
    EXPECT_GT(info, 0.0f);
    EXPECT_LT(info, 1.0f);
}

TEST(InformativenessTest, EmptySite) {
    PsvSite site;
    EXPECT_FLOAT_EQ(ComputeInformativeness(site), 0.0f);
}

// === PSV Catalog Builder tests ===

TEST(PsvCatalogBuilderTest, BuildFromSequences) {
    PsvCatalogBuilder builder;
    builder.SetReferenceSequence("ACGTACGT");
    builder.AddParalogSequence("P1", "ACGTACGT");
    builder.AddParalogSequence("P2", "GCGTACGT");

    auto catalog = builder.Build("chr1", 0, 0.0f);

    EXPECT_GE(catalog.Size(), 1);
    const auto* site = catalog.GetSiteAtPosition("chr1", 0);
    ASSERT_NE(site, nullptr);
    EXPECT_EQ(site->ref_allele, 'A');
}

TEST(PsvCatalogBuilderTest, BuildWithMinInformativeness) {
    PsvCatalogBuilder builder;
    builder.SetReferenceSequence("ACGT");
    builder.AddParalogSequence("P1", "ACGT");
    builder.AddParalogSequence("P2", "GCGT");

    auto catalog_low = builder.Build("chr1", 0, 0.0f);
    auto catalog_high = builder.Build("chr1", 0, 0.99f);

    EXPECT_GE(catalog_low.Size(), catalog_high.Size());
}

// === PSV Catalog I/O tests ===

TEST(PsvCatalogIOTest, SaveAndLoadBed) {
    PsvCatalog catalog;

    PsvSite site1;
    site1.psv_id = 1;
    site1.chrom = "chr1";
    site1.position = 100;
    site1.ref_allele = 'A';
    site1.paralog_alleles["GENE1"] = 'A';
    site1.paralog_alleles["GENE2"] = 'G';
    site1.informativeness = 1.0f;
    catalog.AddSite(site1);

    PsvSite site2;
    site2.psv_id = 2;
    site2.chrom = "chr1";
    site2.position = 200;
    site2.ref_allele = 'C';
    site2.paralog_alleles["GENE1"] = 'C';
    site2.paralog_alleles["GENE2"] = 'T';
    site2.informativeness = 1.0f;
    catalog.AddSite(site2);

    auto temp_dir = std::filesystem::temp_directory_path();
    auto path = temp_dir / "test_psv_catalog.bed";

    ASSERT_TRUE(SavePsvCatalogToBed(catalog, path));

    auto loaded = LoadPsvCatalogFromBed(path);
    ASSERT_TRUE(loaded.has_value());
    EXPECT_EQ(loaded->Size(), 2);

    std::filesystem::remove(path);
}

// === PSV Assigner tests ===

class PsvAssignerTest : public ::testing::Test {
protected:
    void SetUp() override {
        PsvSite site1;
        site1.psv_id = 0;
        site1.chrom = "chr1";
        site1.position = 100;
        site1.ref_allele = 'A';
        site1.paralog_alleles["GENE1"] = 'A';
        site1.paralog_alleles["GENE2"] = 'G';
        site1.informativeness = 1.0f;
        catalog_.AddSite(site1);

        PsvSite site2;
        site2.psv_id = 1;
        site2.chrom = "chr1";
        site2.position = 200;
        site2.ref_allele = 'C';
        site2.paralog_alleles["GENE1"] = 'C';
        site2.paralog_alleles["GENE2"] = 'T';
        site2.informativeness = 1.0f;
        catalog_.AddSite(site2);

        catalog_.BuildIndex();
    }

    PsvCatalog catalog_;
};

TEST_F(PsvAssignerTest, AssignWithMatchingObservations) {
    PsvAssignmentConfig config;
    config.min_base_quality = 20;
    config.confidence_threshold = 0.9f;

    PsvAssigner assigner(catalog_, config);

    std::vector<PsvObservation> obs;
    obs.push_back({0, 'A', 30, 10});
    obs.push_back({1, 'C', 30, 20});

    auto result = assigner.Assign("read1", obs);

    EXPECT_EQ(result.read_id, "read1");
    EXPECT_EQ(result.total_psvs_observed, 2);
    EXPECT_EQ(result.informative_psvs, 2);
    EXPECT_EQ(result.best_paralog, "GENE1");
    EXPECT_GT(result.best_posterior, 0.9f);
    EXPECT_TRUE(result.is_confident);
}

TEST_F(PsvAssignerTest, AssignWithConflictingObservations) {
    PsvAssignmentConfig config;
    config.min_base_quality = 20;
    config.confidence_threshold = 0.95f;

    PsvAssigner assigner(catalog_, config);

    std::vector<PsvObservation> obs;
    obs.push_back({0, 'A', 30, 10});
    obs.push_back({1, 'T', 30, 20});

    auto result = assigner.Assign("read2", obs);

    EXPECT_EQ(result.informative_psvs, 2);
    EXPECT_FALSE(result.is_confident);
}

TEST_F(PsvAssignerTest, AssignEmptyObservations) {
    PsvAssigner assigner(catalog_);

    std::vector<PsvObservation> obs;
    auto result = assigner.Assign("read3", obs);

    EXPECT_EQ(result.total_psvs_observed, 0);
    EXPECT_EQ(result.informative_psvs, 0);
    EXPECT_FALSE(result.is_confident);
}

TEST_F(PsvAssignerTest, AssignLowQualityFiltered) {
    PsvAssignmentConfig config;
    config.min_base_quality = 30;

    PsvAssigner assigner(catalog_, config);

    std::vector<PsvObservation> obs;
    obs.push_back({0, 'A', 10, 10});
    obs.push_back({1, 'C', 10, 20});

    auto result = assigner.Assign("read4", obs);

    EXPECT_EQ(result.informative_psvs, 0);
    EXPECT_FALSE(result.is_confident);
}

TEST_F(PsvAssignerTest, GetStats) {
    PsvAssigner assigner(catalog_);

    std::vector<PsvObservation> obs1 = {{0, 'A', 30, 10}};
    std::vector<PsvObservation> obs2 = {{0, 'G', 30, 10}};

    (void)assigner.Assign("read1", obs1);
    (void)assigner.Assign("read2", obs2);

    const auto& stats = assigner.GetStats();
    EXPECT_EQ(stats.reads_processed, 2);
    EXPECT_EQ(stats.reads_with_psvs, 2);
}

TEST_F(PsvAssignerTest, ResetStats) {
    PsvAssigner assigner(catalog_);

    std::vector<PsvObservation> obs = {{0, 'A', 30, 10}};
    (void)assigner.Assign("read1", obs);

    EXPECT_EQ(assigner.GetStats().reads_processed, 1);

    assigner.ResetStats();
    EXPECT_EQ(assigner.GetStats().reads_processed, 0);
}

// === ResultToParalogCall tests ===

TEST(ResultToParalogCallTest, ConvertResult) {
    PsvAssignmentResult result;
    result.read_id = "read1";
    result.informative_psvs = 5;

    ParalogLikelihood pl1;
    pl1.paralog_id = "GENE1";
    pl1.posterior = 0.7f;
    result.likelihoods.push_back(pl1);

    ParalogLikelihood pl2;
    pl2.paralog_id = "GENE2";
    pl2.posterior = 0.3f;
    result.likelihoods.push_back(pl2);

    auto call = ResultToParalogCall(result);

    EXPECT_EQ(call.n_discriminating_psvs, 5);
    EXPECT_EQ(call.inter_paralog.size(), 2);
    EXPECT_EQ(call.inter_paralog[0].first, "GENE1");
    EXPECT_FLOAT_EQ(call.inter_paralog[0].second, 0.7f);
}

// === MergeAssignments tests ===

TEST(MergeAssignmentsTest, MergeTwoAssignments) {
    ParalogCall existing;
    existing.n_discriminating_psvs = 3;
    existing.inter_paralog.emplace_back("GENE1", 0.8f);
    existing.inter_paralog.emplace_back("GENE2", 0.2f);

    PsvAssignmentResult psv_result;
    psv_result.informative_psvs = 2;
    psv_result.likelihoods.push_back({"GENE1", 0.0f, 0.6f, 0, 0});
    psv_result.likelihoods.push_back({"GENE2", 0.0f, 0.4f, 0, 0});

    auto merged = MergeAssignments(existing, psv_result);

    EXPECT_EQ(merged.n_discriminating_psvs, 5);
    EXPECT_EQ(merged.inter_paralog.size(), 2);

    float sum = 0.0f;
    for (const auto& [_, prob] : merged.inter_paralog) {
        sum += prob;
    }
    EXPECT_NEAR(sum, 1.0f, 0.01f);
}

// === Extract observations tests ===

TEST_F(PsvAssignerTest, ExtractObservationsFromHit) {
    PsvAssigner assigner(catalog_);

    AlignmentHit hit;
    hit.target_id = "chr1";
    hit.start = 90;
    hit.end = 210;
    hit.cigar.ops = "120M";

    std::string read_seq(120, 'A');
    read_seq[10] = 'A';
    read_seq[110] = 'C';

    auto obs = assigner.ExtractObservations(hit, read_seq);

    EXPECT_GE(obs.size(), 2);
}

// === Update AlignmentRecord tests ===

TEST_F(PsvAssignerTest, UpdateRecordWithPSV) {
    PsvAssigner assigner(catalog_);

    AlignmentRecord record;
    record.read_id = "read1";
    record.read_len = 120;
    record.status = AlignmentStatus::Mapped;

    AlignmentHit hit;
    hit.target_id = "chr1";
    hit.start = 90;
    hit.end = 210;
    hit.cigar.ops = "120M";
    record.primary = hit;

    std::string read_seq(120, 'A');
    read_seq[10] = 'A';
    read_seq[110] = 'C';

    assigner.UpdateRecord(record, read_seq);

    EXPECT_TRUE(record.paralog_assignment.has_value());
}

}  // namespace
}  // namespace llmap::psv
