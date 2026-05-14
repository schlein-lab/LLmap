// LLmap — Unit tests for single-cell paralog QC metrics.

#include "singlecell/sc_paralog_qc.h"
#include "core/alignment_record.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <random>

namespace llmap::singlecell::test {

class ScParalogQcTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::random_device rd;
        test_dir_ = std::filesystem::temp_directory_path() /
                    ("sc_qc_test_" + std::to_string(rd()));
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    CellParalogMatrix CreateTestMatrix() {
        CellParalogMatrix matrix;

        for (int cell = 0; cell < 10; ++cell) {
            for (int paralog = 0; paralog < 3; ++paralog) {
                AlignmentRecord rec;
                rec.cell_barcode = "CELL" + std::to_string(cell);
                ParalogCall pc;
                std::string paralog_id = "PARALOG" + std::to_string(paralog);
                float prob = 0.5f + 0.05f * static_cast<float>(cell);
                pc.inter_paralog.emplace_back(paralog_id, prob);
                rec.paralog_assignment = pc;
                rec.confidence_scores.push_back(prob);
                for (int i = 0; i < cell + 1; ++i) {
                    matrix.AddRecord(rec);
                }
            }
        }

        CellParalogConfig config;
        config.aggregation = CellParalogConfig::AggregationMethod::Mean;
        matrix.Finalize(config);
        return matrix;
    }

    std::filesystem::path test_dir_;
};

TEST_F(ScParalogQcTest, ComputeEntropy_Uniform) {
    std::vector<float> probs = {0.25f, 0.25f, 0.25f, 0.25f};
    float entropy = ComputeEntropy(probs);
    EXPECT_NEAR(entropy, 2.0f, 0.01f);
}

TEST_F(ScParalogQcTest, ComputeEntropy_Certain) {
    std::vector<float> probs = {1.0f, 0.0f, 0.0f};
    float entropy = ComputeEntropy(probs);
    EXPECT_NEAR(entropy, 0.0f, 0.01f);
}

TEST_F(ScParalogQcTest, ComputeEntropy_Binary) {
    std::vector<float> probs = {0.5f, 0.5f};
    float entropy = ComputeEntropy(probs);
    EXPECT_NEAR(entropy, 1.0f, 0.01f);
}

TEST_F(ScParalogQcTest, ComputeEntropy_Empty) {
    std::vector<float> probs;
    float entropy = ComputeEntropy(probs);
    EXPECT_FLOAT_EQ(entropy, 0.0f);
}

TEST_F(ScParalogQcTest, ComputeDominance_Single) {
    std::vector<float> probs = {0.9f};
    float dom = ComputeDominance(probs);
    EXPECT_FLOAT_EQ(dom, 1.0f);
}

TEST_F(ScParalogQcTest, ComputeDominance_ClearWinner) {
    std::vector<float> probs = {0.8f, 0.2f};
    float dom = ComputeDominance(probs);
    EXPECT_FLOAT_EQ(dom, 4.0f);
}

TEST_F(ScParalogQcTest, ComputeDominance_Tie) {
    std::vector<float> probs = {0.5f, 0.5f};
    float dom = ComputeDominance(probs);
    EXPECT_FLOAT_EQ(dom, 1.0f);
}

TEST_F(ScParalogQcTest, ComputeDominance_Empty) {
    std::vector<float> probs;
    float dom = ComputeDominance(probs);
    EXPECT_FLOAT_EQ(dom, 0.0f);
}

TEST_F(ScParalogQcTest, ConfidenceDistribution_Basic) {
    std::vector<float> confs = {0.1f, 0.2f, 0.3f, 0.7f, 0.8f, 0.9f};
    auto dist = ComputeConfidenceDistribution(confs, 0.5f, 0.1f);

    EXPECT_EQ(dist.below_threshold, 3);
    EXPECT_EQ(dist.above_threshold, 3);
    EXPECT_NEAR(dist.mean, 0.5f, 0.01f);
    EXPECT_EQ(dist.histogram.size(), 10);
}

TEST_F(ScParalogQcTest, ConfidenceDistribution_Empty) {
    std::vector<float> confs;
    auto dist = ComputeConfidenceDistribution(confs);

    EXPECT_EQ(dist.below_threshold, 0);
    EXPECT_EQ(dist.above_threshold, 0);
    EXPECT_FLOAT_EQ(dist.mean, 0.0f);
}

TEST_F(ScParalogQcTest, ConfidenceDistribution_Histogram) {
    std::vector<float> confs = {0.05f, 0.15f, 0.95f};
    auto dist = ComputeConfidenceDistribution(confs, 0.5f, 0.1f);

    EXPECT_EQ(dist.histogram[0], 1);
    EXPECT_EQ(dist.histogram[1], 1);
    EXPECT_EQ(dist.histogram[9], 1);
}

TEST_F(ScParalogQcTest, ComputeCellQcMetrics_Basic) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    EXPECT_EQ(metrics.size(), 10);
    for (const auto& m : metrics) {
        EXPECT_FALSE(m.cell_barcode.empty());
        EXPECT_EQ(m.paralogs_detected, 3);
        EXPECT_GT(m.total_reads, 0);
    }
}

TEST_F(ScParalogQcTest, ComputeCellQcMetrics_AssignmentRate) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    for (const auto& m : metrics) {
        EXPECT_GE(m.assignment_rate, 0.0f);
        EXPECT_LE(m.assignment_rate, 1.0f);
    }
}

TEST_F(ScParalogQcTest, ComputeParalogQcMetrics_Basic) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeParalogQcMetrics(matrix);

    EXPECT_EQ(metrics.size(), 3);
    for (const auto& m : metrics) {
        EXPECT_FALSE(m.paralog_id.empty());
        EXPECT_EQ(m.total_cells, 10);
        EXPECT_GT(m.total_reads, 0);
    }
}

TEST_F(ScParalogQcTest, ComputeParalogQcMetrics_DetectionRate) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeParalogQcMetrics(matrix);

    for (const auto& m : metrics) {
        EXPECT_FLOAT_EQ(m.detection_rate, 1.0f);
    }
}

TEST_F(ScParalogQcTest, ComputeGlobalQcSummary_Basic) {
    auto matrix = CreateTestMatrix();
    auto cell_metrics = ComputeCellQcMetrics(matrix);
    auto paralog_metrics = ComputeParalogQcMetrics(matrix);
    auto summary = ComputeGlobalQcSummary(matrix, cell_metrics, paralog_metrics);

    EXPECT_EQ(summary.total_cells, 10);
    EXPECT_EQ(summary.total_paralogs, 3);
    EXPECT_FLOAT_EQ(summary.mean_paralogs_per_cell, 3.0f);
}

TEST_F(ScParalogQcTest, GenerateQcReport_Basic) {
    auto matrix = CreateTestMatrix();
    auto report = GenerateQcReport(matrix, {}, "test_sample");

    EXPECT_EQ(report.sample_id, "test_sample");
    EXPECT_FALSE(report.timestamp.empty());
    EXPECT_EQ(report.cells.size(), 10);
    EXPECT_EQ(report.paralogs.size(), 3);
    EXPECT_EQ(report.global.total_cells, 10);
}

TEST_F(ScParalogQcTest, GenerateQcReport_WithThresholds) {
    auto matrix = CreateTestMatrix();
    QcThresholds thresh;
    thresh.min_reads_per_cell = 5;
    auto report = GenerateQcReport(matrix, thresh, "filtered_sample");

    EXPECT_GE(report.global.cells_passing_qc, 0);
    EXPECT_LE(report.global.cells_passing_qc, 10);
}

TEST_F(ScParalogQcTest, GetCellsPassingQc_AllPass) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    QcThresholds thresh;
    thresh.min_reads_per_cell = 1;
    thresh.min_confidence = 0.0f;
    thresh.min_assignment_rate = 0.0f;
    thresh.max_entropy = 100.0f;

    auto passing = GetCellsPassingQc(metrics, thresh);
    EXPECT_EQ(passing.size(), 10);
}

TEST_F(ScParalogQcTest, GetCellsPassingQc_SomePass) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    QcThresholds thresh;
    thresh.min_reads_per_cell = 15;
    thresh.min_confidence = 0.0f;
    thresh.min_assignment_rate = 0.0f;
    thresh.max_entropy = 100.0f;

    auto passing = GetCellsPassingQc(metrics, thresh);
    EXPECT_LT(passing.size(), 10);
    EXPECT_GT(passing.size(), 0);
}

TEST_F(ScParalogQcTest, FilterMatrixByQc_Filters) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    QcThresholds thresh;
    thresh.min_reads_per_cell = 25;
    auto passing = GetCellsPassingQc(metrics, thresh);

    auto filtered = FilterMatrixByQc(matrix, passing);
    EXPECT_LT(filtered.GetCells().size(), 10);
}

TEST_F(ScParalogQcTest, ExportQcReportJson_CreatesFile) {
    auto matrix = CreateTestMatrix();
    auto report = GenerateQcReport(matrix, {}, "json_test");

    auto path = test_dir_ / "report.json";
    EXPECT_TRUE(ExportQcReportJson(report, path));
    EXPECT_TRUE(std::filesystem::exists(path));

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_TRUE(content.find("\"sample_id\": \"json_test\"") != std::string::npos);
    EXPECT_TRUE(content.find("\"total_cells\": 10") != std::string::npos);
}

TEST_F(ScParalogQcTest, ExportQcReportJson_ValidJson) {
    auto matrix = CreateTestMatrix();
    auto report = GenerateQcReport(matrix);

    auto path = test_dir_ / "report2.json";
    EXPECT_TRUE(ExportQcReportJson(report, path));

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_TRUE(content.front() == '{');
    EXPECT_TRUE(content.find("}\n") != std::string::npos);
}

TEST_F(ScParalogQcTest, ExportCellQcTsv_CreatesFile) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    auto path = test_dir_ / "cells.tsv";
    EXPECT_TRUE(ExportCellQcTsv(metrics, path));
    EXPECT_TRUE(std::filesystem::exists(path));

    std::ifstream in(path);
    std::string header;
    std::getline(in, header);
    EXPECT_TRUE(header.find("cell_barcode") != std::string::npos);
    EXPECT_TRUE(header.find("assignment_entropy") != std::string::npos);
}

TEST_F(ScParalogQcTest, ExportCellQcTsv_CorrectRows) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    auto path = test_dir_ / "cells2.tsv";
    EXPECT_TRUE(ExportCellQcTsv(metrics, path));

    std::ifstream in(path);
    int lines = 0;
    std::string line;
    while (std::getline(in, line)) lines++;
    EXPECT_EQ(lines, 11);
}

TEST_F(ScParalogQcTest, ExportParalogQcTsv_CreatesFile) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeParalogQcMetrics(matrix);

    auto path = test_dir_ / "paralogs.tsv";
    EXPECT_TRUE(ExportParalogQcTsv(metrics, path));
    EXPECT_TRUE(std::filesystem::exists(path));
}

TEST_F(ScParalogQcTest, ExportSummaryTsv_CreatesFile) {
    auto matrix = CreateTestMatrix();
    auto cell_metrics = ComputeCellQcMetrics(matrix);
    auto paralog_metrics = ComputeParalogQcMetrics(matrix);
    auto summary = ComputeGlobalQcSummary(matrix, cell_metrics, paralog_metrics);

    auto path = test_dir_ / "summary.tsv";
    EXPECT_TRUE(ExportSummaryTsv(summary, path));
    EXPECT_TRUE(std::filesystem::exists(path));

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_TRUE(content.find("total_cells\t10") != std::string::npos);
}

TEST_F(ScParalogQcTest, ExportQcReportTsv_CreatesAllFiles) {
    auto matrix = CreateTestMatrix();
    auto report = GenerateQcReport(matrix);

    auto output_dir = test_dir_ / "tsv_report";
    EXPECT_TRUE(ExportQcReportTsv(report, output_dir));
    EXPECT_TRUE(std::filesystem::exists(output_dir / "cells_qc.tsv"));
    EXPECT_TRUE(std::filesystem::exists(output_dir / "paralogs_qc.tsv"));
    EXPECT_TRUE(std::filesystem::exists(output_dir / "summary_qc.tsv"));
}

TEST_F(ScParalogQcTest, CellQcMetrics_Entropy) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    for (const auto& m : metrics) {
        EXPECT_GE(m.assignment_entropy, 0.0f);
    }
}

TEST_F(ScParalogQcTest, CellQcMetrics_Dominance) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeCellQcMetrics(matrix);

    for (const auto& m : metrics) {
        EXPECT_GE(m.dominance_score, 0.0f);
    }
}

TEST_F(ScParalogQcTest, ParalogQcMetrics_CV) {
    auto matrix = CreateTestMatrix();
    auto metrics = ComputeParalogQcMetrics(matrix);

    for (const auto& m : metrics) {
        EXPECT_GE(m.coefficient_of_variation, 0.0f);
    }
}

TEST_F(ScParalogQcTest, GlobalSummary_AssignmentRate) {
    auto matrix = CreateTestMatrix();
    auto cell_metrics = ComputeCellQcMetrics(matrix);
    auto paralog_metrics = ComputeParalogQcMetrics(matrix);
    auto summary = ComputeGlobalQcSummary(matrix, cell_metrics, paralog_metrics);

    EXPECT_GE(summary.global_assignment_rate, 0.0f);
    EXPECT_LE(summary.global_assignment_rate, 1.0f);
}

TEST_F(ScParalogQcTest, GlobalSummary_ConfidenceDistribution) {
    auto matrix = CreateTestMatrix();
    auto cell_metrics = ComputeCellQcMetrics(matrix);
    auto paralog_metrics = ComputeParalogQcMetrics(matrix);
    auto summary = ComputeGlobalQcSummary(matrix, cell_metrics, paralog_metrics);

    EXPECT_FALSE(summary.confidence_dist.histogram.empty());
    EXPECT_GE(summary.confidence_dist.mean, 0.0f);
    EXPECT_LE(summary.confidence_dist.mean, 1.0f);
}

TEST_F(ScParalogQcTest, EmptyMatrix_Metrics) {
    CellParalogMatrix matrix;
    CellParalogConfig config;
    matrix.Finalize(config);

    auto cell_metrics = ComputeCellQcMetrics(matrix);
    auto paralog_metrics = ComputeParalogQcMetrics(matrix);

    EXPECT_TRUE(cell_metrics.empty());
    EXPECT_TRUE(paralog_metrics.empty());
}

TEST_F(ScParalogQcTest, EmptyMatrix_GlobalSummary) {
    CellParalogMatrix matrix;
    CellParalogConfig config;
    matrix.Finalize(config);

    auto cell_metrics = ComputeCellQcMetrics(matrix);
    auto paralog_metrics = ComputeParalogQcMetrics(matrix);
    auto summary = ComputeGlobalQcSummary(matrix, cell_metrics, paralog_metrics);

    EXPECT_EQ(summary.total_cells, 0);
    EXPECT_EQ(summary.total_paralogs, 0);
    EXPECT_FLOAT_EQ(summary.global_assignment_rate, 0.0f);
}

TEST_F(ScParalogQcTest, QcReport_Timestamp) {
    auto matrix = CreateTestMatrix();
    auto report = GenerateQcReport(matrix);

    EXPECT_FALSE(report.timestamp.empty());
    EXPECT_TRUE(report.timestamp.find("T") != std::string::npos);
    EXPECT_TRUE(report.timestamp.find("Z") != std::string::npos);
}

}  // namespace llmap::singlecell::test
