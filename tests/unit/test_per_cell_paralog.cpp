// LLmap — Tests for per-cell paralog matrix.

#include <gtest/gtest.h>

#include "singlecell/per_cell_paralog.h"

#include <filesystem>
#include <fstream>

namespace llmap::singlecell {
namespace {

AlignmentRecord MakeRecordWithParalog(
    const std::string& read_id,
    const std::string& cell_barcode,
    const std::vector<std::pair<std::string, float>>& paralogs) {
    AlignmentRecord record;
    record.read_id = read_id;
    record.read_len = 100;
    record.status = AlignmentStatus::Mapped;
    record.cell_barcode = cell_barcode;

    ParalogCall call;
    call.inter_paralog = paralogs;
    record.paralog_assignment = call;

    return record;
}

class PerCellParalogTest : public ::testing::Test {
protected:
    void TearDown() override {
        for (const auto& path : temp_files_) {
            std::filesystem::remove(path);
        }
    }

    std::filesystem::path CreateTempFile(const std::string& suffix) {
        auto path = std::filesystem::temp_directory_path() /
                    ("test_per_cell_paralog_" + std::to_string(counter_++) + suffix);
        temp_files_.push_back(path);
        return path;
    }

private:
    std::vector<std::filesystem::path> temp_files_;
    int counter_{0};
};

TEST_F(PerCellParalogTest, EmptyMatrix) {
    CellParalogMatrix matrix;
    CellParalogConfig config;
    matrix.Finalize(config);

    EXPECT_TRUE(matrix.IsFinalized());
    EXPECT_TRUE(matrix.GetEntries().empty());
    EXPECT_TRUE(matrix.GetCells().empty());
    EXPECT_TRUE(matrix.GetParalogs().empty());

    auto stats = matrix.GetStats();
    EXPECT_EQ(stats.unique_cells, 0);
    EXPECT_EQ(stats.unique_paralogs, 0);
    EXPECT_EQ(stats.matrix_entries, 0);
}

TEST_F(PerCellParalogTest, SingleRecord) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.8f}}));

    CellParalogConfig config;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    ASSERT_EQ(entries.size(), 1);
    EXPECT_EQ(entries[0].cell_barcode, "ACGT");
    EXPECT_EQ(entries[0].paralog_id, "IGHG1");
    EXPECT_FLOAT_EQ(entries[0].probability, 0.8f);
    EXPECT_EQ(entries[0].read_count, 1);
}

TEST_F(PerCellParalogTest, MultipleReadsOneCell) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.6f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "ACGT", {{"IGHG1", 0.8f}}));
    matrix.AddRecord(MakeRecordWithParalog("read3", "ACGT", {{"IGHG1", 0.4f}}));

    CellParalogConfig config;
    config.aggregation = CellParalogConfig::AggregationMethod::Mean;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    ASSERT_EQ(entries.size(), 1);
    EXPECT_EQ(entries[0].read_count, 3);
    EXPECT_NEAR(entries[0].probability, 0.6f, 0.001f);
}

TEST_F(PerCellParalogTest, AggregationMax) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.6f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "ACGT", {{"IGHG1", 0.9f}}));
    matrix.AddRecord(MakeRecordWithParalog("read3", "ACGT", {{"IGHG1", 0.4f}}));

    CellParalogConfig config;
    config.aggregation = CellParalogConfig::AggregationMethod::Max;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    ASSERT_EQ(entries.size(), 1);
    EXPECT_FLOAT_EQ(entries[0].probability, 0.9f);
}

TEST_F(PerCellParalogTest, AggregationSum) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.3f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "ACGT", {{"IGHG1", 0.4f}}));

    CellParalogConfig config;
    config.aggregation = CellParalogConfig::AggregationMethod::Sum;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    ASSERT_EQ(entries.size(), 1);
    EXPECT_FLOAT_EQ(entries[0].probability, 0.7f);
}

TEST_F(PerCellParalogTest, MultipleParalogsOneCell) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT",
        {{"IGHG1", 0.5f}, {"IGHG2", 0.3f}, {"IGHG3", 0.2f}}));

    CellParalogConfig config;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    EXPECT_EQ(entries.size(), 3);

    auto paralogs = matrix.GetParalogs();
    EXPECT_EQ(paralogs.size(), 3);
}

TEST_F(PerCellParalogTest, MultipleCells) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "CELL1", {{"IGHG1", 0.9f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "CELL2", {{"IGHG2", 0.8f}}));
    matrix.AddRecord(MakeRecordWithParalog("read3", "CELL3", {{"IGHG1", 0.7f}}));

    CellParalogConfig config;
    matrix.Finalize(config);

    auto cells = matrix.GetCells();
    EXPECT_EQ(cells.size(), 3);

    auto stats = matrix.GetStats();
    EXPECT_EQ(stats.unique_cells, 3);
    EXPECT_EQ(stats.unique_paralogs, 2);
}

TEST_F(PerCellParalogTest, MinProbabilityFilter) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT",
        {{"IGHG1", 0.5f}, {"IGHG2", 0.005f}}));

    CellParalogConfig config;
    config.min_probability = 0.01f;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    ASSERT_EQ(entries.size(), 1);
    EXPECT_EQ(entries[0].paralog_id, "IGHG1");
}

TEST_F(PerCellParalogTest, MinReadsPerCellFilter) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "CELL1", {{"IGHG1", 0.9f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "CELL2", {{"IGHG1", 0.8f}}));
    matrix.AddRecord(MakeRecordWithParalog("read3", "CELL2", {{"IGHG1", 0.7f}}));
    matrix.AddRecord(MakeRecordWithParalog("read4", "CELL2", {{"IGHG1", 0.6f}}));

    CellParalogConfig config;
    config.min_reads_per_cell = 3;
    matrix.Finalize(config);

    auto cells = matrix.GetCells();
    ASSERT_EQ(cells.size(), 1);
    EXPECT_EQ(cells[0], "CELL2");
}

TEST_F(PerCellParalogTest, NormalizeRows) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT",
        {{"IGHG1", 0.6f}, {"IGHG2", 0.4f}}));

    CellParalogConfig config;
    config.normalize_rows = true;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto entries = matrix.GetEntries();
    ASSERT_EQ(entries.size(), 2);

    float sum = 0.0f;
    for (const auto& e : entries) {
        sum += e.probability;
    }
    EXPECT_NEAR(sum, 1.0f, 0.001f);
}

TEST_F(PerCellParalogTest, GetEntriesForCell) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "CELL1", {{"IGHG1", 0.9f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "CELL2", {{"IGHG2", 0.8f}}));

    CellParalogConfig config;
    matrix.Finalize(config);

    auto entries = matrix.GetEntriesForCell("CELL1");
    ASSERT_EQ(entries.size(), 1);
    EXPECT_EQ(entries[0].paralog_id, "IGHG1");

    auto none = matrix.GetEntriesForCell("NOTEXIST");
    EXPECT_TRUE(none.empty());
}

TEST_F(PerCellParalogTest, GetDenseMatrix) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "A", {{"P1", 0.1f}, {"P2", 0.2f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "B", {{"P1", 0.3f}}));

    CellParalogConfig config;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto [dense, cells, paralogs] = matrix.GetDenseMatrix();

    EXPECT_EQ(cells.size(), 2);
    EXPECT_EQ(paralogs.size(), 2);
    EXPECT_EQ(dense.size(), 2);
    EXPECT_EQ(dense[0].size(), 2);
}

TEST_F(PerCellParalogTest, BuilderPattern) {
    auto matrix = CellParalogMatrixBuilder()
        .AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.8f}}))
        .AddRecord(MakeRecordWithParalog("read2", "ACGT", {{"IGHG2", 0.7f}}))
        .Build();

    EXPECT_TRUE(matrix.IsFinalized());
    EXPECT_EQ(matrix.GetEntries().size(), 2);
}

TEST_F(PerCellParalogTest, RecordWithoutBarcode) {
    CellParalogMatrix matrix;

    AlignmentRecord record;
    record.read_id = "read1";
    record.status = AlignmentStatus::Mapped;
    ParalogCall call;
    call.inter_paralog = {{"IGHG1", 0.9f}};
    record.paralog_assignment = call;

    matrix.AddRecord(record);

    CellParalogConfig config;
    matrix.Finalize(config);

    EXPECT_TRUE(matrix.GetEntries().empty());
}

TEST_F(PerCellParalogTest, RecordWithoutParalog) {
    CellParalogMatrix matrix;

    AlignmentRecord record;
    record.read_id = "read1";
    record.status = AlignmentStatus::Mapped;
    record.cell_barcode = "ACGT";

    matrix.AddRecord(record);

    CellParalogConfig config;
    matrix.Finalize(config);

    EXPECT_TRUE(matrix.GetEntries().empty());
}

TEST_F(PerCellParalogTest, ExportToCSV) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.8f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "TGCA", {{"IGHG2", 0.7f}}));

    CellParalogConfig config;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto path = CreateTempFile(".csv");
    EXPECT_TRUE(ExportToCSV(matrix, path));

    std::ifstream file(path);
    EXPECT_TRUE(file.is_open());

    std::string line;
    std::getline(file, line);
    EXPECT_EQ(line, "cell_barcode,paralog_id,probability,read_count,mean_confidence");

    int data_lines = 0;
    while (std::getline(file, line)) {
        ++data_lines;
    }
    EXPECT_EQ(data_lines, 2);
}

TEST_F(PerCellParalogTest, ExportToTSV) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.8f}}));

    CellParalogConfig config;
    matrix.Finalize(config);

    auto path = CreateTempFile(".tsv");
    EXPECT_TRUE(ExportToTSV(matrix, path));

    std::ifstream file(path);
    std::string line;
    std::getline(file, line);
    EXPECT_NE(line.find('\t'), std::string::npos);
}

TEST_F(PerCellParalogTest, ExportToDenseCSV) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "A", {{"P1", 0.1f}, {"P2", 0.2f}}));
    matrix.AddRecord(MakeRecordWithParalog("read2", "B", {{"P1", 0.3f}}));

    CellParalogConfig config;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto path = CreateTempFile(".csv");
    EXPECT_TRUE(ExportToDenseCSV(matrix, path));

    std::ifstream file(path);
    std::string line;

    std::getline(file, line);
    EXPECT_EQ(line.substr(0, 12), "cell_barcode");

    int data_lines = 0;
    while (std::getline(file, line)) {
        ++data_lines;
    }
    EXPECT_EQ(data_lines, 2);
}

TEST_F(PerCellParalogTest, ImportFromCSV) {
    auto path = CreateTempFile(".csv");

    {
        std::ofstream file(path);
        file << "cell_barcode,paralog_id,probability,read_count,mean_confidence\n";
        file << "ACGT,IGHG1,0.8,5,0.9\n";
        file << "TGCA,IGHG2,0.7,3,0.85\n";
    }

    auto matrix = ImportFromCSV(path);
    ASSERT_TRUE(matrix.has_value());

    auto entries = matrix->GetEntries();
    EXPECT_EQ(entries.size(), 2);
}

TEST_F(PerCellParalogTest, ImportFromCSVMissingFile) {
    auto result = ImportFromCSV("/nonexistent/path/file.csv");
    EXPECT_FALSE(result.has_value());
}

TEST_F(PerCellParalogTest, ExportToH5ADNotSupported) {
    CellParalogMatrix matrix;
    CellParalogConfig config;
    matrix.Finalize(config);

    auto path = CreateTempFile(".h5ad");
    EXPECT_FALSE(ExportToH5AD(matrix, path));
}

TEST_F(PerCellParalogTest, SparsityCalculation) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("r1", "C1", {{"P1", 0.5f}}));
    matrix.AddRecord(MakeRecordWithParalog("r2", "C1", {{"P2", 0.5f}}));
    matrix.AddRecord(MakeRecordWithParalog("r3", "C2", {{"P1", 0.5f}}));

    CellParalogConfig config;
    config.min_probability = 0.0f;
    matrix.Finalize(config);

    auto stats = matrix.GetStats();
    EXPECT_EQ(stats.unique_cells, 2);
    EXPECT_EQ(stats.unique_paralogs, 2);
    EXPECT_EQ(stats.matrix_entries, 3);
    EXPECT_NEAR(stats.sparsity, 0.25f, 0.001f);
}

TEST_F(PerCellParalogTest, ClearMatrix) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT", {{"IGHG1", 0.8f}}));

    CellParalogConfig config;
    matrix.Finalize(config);
    EXPECT_TRUE(matrix.IsFinalized());

    matrix.Clear();
    EXPECT_FALSE(matrix.IsFinalized());
    EXPECT_TRUE(matrix.GetEntries().empty());
}

TEST_F(PerCellParalogTest, ParalogAccumulatorBasic) {
    ParalogAccumulator acc;
    acc.Add(0.5f, 0.9f);
    acc.Add(0.7f, 0.8f);

    EXPECT_EQ(acc.count, 2);
    EXPECT_FLOAT_EQ(acc.GetMeanConfidence(), 0.85f);
}

TEST_F(PerCellParalogTest, AddRecordsVector) {
    std::vector<AlignmentRecord> records;
    records.push_back(MakeRecordWithParalog("r1", "C1", {{"P1", 0.5f}}));
    records.push_back(MakeRecordWithParalog("r2", "C2", {{"P2", 0.6f}}));

    CellParalogMatrix matrix;
    matrix.AddRecords(records);

    CellParalogConfig config;
    matrix.Finalize(config);

    EXPECT_EQ(matrix.GetCells().size(), 2);
}

TEST_F(PerCellParalogTest, CSVWithSpecialCharacters) {
    CellParalogMatrix matrix;
    matrix.AddRecord(MakeRecordWithParalog("read1", "ACGT,special", {{"IGHG\"1", 0.8f}}));

    CellParalogConfig config;
    matrix.Finalize(config);

    auto path = CreateTempFile(".csv");
    EXPECT_TRUE(ExportToCSV(matrix, path));

    auto imported = ImportFromCSV(path);
    ASSERT_TRUE(imported.has_value());
    EXPECT_EQ(imported->GetEntries().size(), 1);
}

}  // namespace
}  // namespace llmap::singlecell
