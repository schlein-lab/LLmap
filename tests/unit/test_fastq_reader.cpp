// Unit tests for FastqReader

#include <gtest/gtest.h>

#include "io/fastq_reader.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

namespace llmap::io {
namespace {

class FastqReaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_fastq_test";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path CreateTestFastq(const std::string& name,
                                           const std::string& content) {
        auto path = test_dir_ / name;
        std::ofstream out(path);
        out << content;
        return path;
    }

    std::filesystem::path test_dir_;
};

// ========== Basic Functionality ==========

TEST_F(FastqReaderTest, OpenNonexistentFile) {
    auto reader = FastqReader::Open("/nonexistent/file.fastq");
    EXPECT_EQ(reader, nullptr);
}

TEST_F(FastqReaderTest, ReadEmptyFile) {
    auto path = CreateTestFastq("empty.fastq", "");
    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    EXPECT_FALSE(record.has_value());
}

TEST_F(FastqReaderTest, ReadSingleRecord) {
    auto path = CreateTestFastq("single.fastq",
        "@read1\n"
        "ACGTACGT\n"
        "+\n"
        "IIIIIIII\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->id, "read1");
    EXPECT_EQ(record->sequence, "ACGTACGT");
    EXPECT_EQ(record->quality, "IIIIIIII");
    EXPECT_EQ(record->Length(), 8);

    auto next = reader->Next();
    EXPECT_FALSE(next.has_value());
}

TEST_F(FastqReaderTest, ReadMultipleRecords) {
    auto path = CreateTestFastq("multi.fastq",
        "@read1\n"
        "ACGT\n"
        "+\n"
        "IIII\n"
        "@read2\n"
        "TGCA\n"
        "+\n"
        "HHHH\n"
        "@read3\n"
        "GGCC\n"
        "+\n"
        "FFFF\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto records = reader->ReadAll();
    ASSERT_EQ(records.size(), 3);

    EXPECT_EQ(records[0].id, "read1");
    EXPECT_EQ(records[0].sequence, "ACGT");

    EXPECT_EQ(records[1].id, "read2");
    EXPECT_EQ(records[1].sequence, "TGCA");

    EXPECT_EQ(records[2].id, "read3");
    EXPECT_EQ(records[2].sequence, "GGCC");
}

TEST_F(FastqReaderTest, NextBatch) {
    auto path = CreateTestFastq("batch.fastq",
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nTGCA\n+\nHHHH\n"
        "@read3\nGGCC\n+\nFFFF\n"
        "@read4\nCCAA\n+\nEEEE\n"
        "@read5\nTTTT\n+\nDDDD\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto batch1 = reader->NextBatch(2);
    ASSERT_EQ(batch1.size(), 2);
    EXPECT_EQ(batch1[0].id, "read1");
    EXPECT_EQ(batch1[1].id, "read2");

    auto batch2 = reader->NextBatch(2);
    ASSERT_EQ(batch2.size(), 2);
    EXPECT_EQ(batch2[0].id, "read3");
    EXPECT_EQ(batch2[1].id, "read4");

    auto batch3 = reader->NextBatch(2);
    ASSERT_EQ(batch3.size(), 1);
    EXPECT_EQ(batch3[0].id, "read5");

    auto batch4 = reader->NextBatch(2);
    EXPECT_TRUE(batch4.empty());
}

// ========== Sequence Handling ==========

TEST_F(FastqReaderTest, UppercaseSequence) {
    auto path = CreateTestFastq("lower.fastq",
        "@read1\n"
        "acgtACGT\n"
        "+\n"
        "IIIIIIII\n"
    );

    FastqReaderConfig config;
    config.uppercase_sequence = true;

    auto reader = FastqReader::Open(path, config);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->sequence, "ACGTACGT");
}

TEST_F(FastqReaderTest, PreserveLowercase) {
    auto path = CreateTestFastq("lower.fastq",
        "@read1\n"
        "acgtACGT\n"
        "+\n"
        "IIIIIIII\n"
    );

    FastqReaderConfig config;
    config.uppercase_sequence = false;

    auto reader = FastqReader::Open(path, config);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->sequence, "acgtACGT");
}

TEST_F(FastqReaderTest, HandleNBases) {
    auto path = CreateTestFastq("n_bases.fastq",
        "@read1\n"
        "ACNTNACG\n"
        "+\n"
        "IIIIIIII\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->sequence, "ACNTNACG");
}

// ========== Configuration Options ==========

TEST_F(FastqReaderTest, MaxRecordsLimit) {
    auto path = CreateTestFastq("many.fastq",
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nTGCA\n+\nHHHH\n"
        "@read3\nGGCC\n+\nFFFF\n"
    );

    FastqReaderConfig config;
    config.max_records = 2;

    auto reader = FastqReader::Open(path, config);
    ASSERT_NE(reader, nullptr);

    auto records = reader->ReadAll();
    EXPECT_EQ(records.size(), 2);
}

TEST_F(FastqReaderTest, MinLengthFilter) {
    auto path = CreateTestFastq("varied.fastq",
        "@read1\nAC\n+\nII\n"
        "@read2\nACGT\n+\nIIII\n"
        "@read3\nACGTACGT\n+\nIIIIIIII\n"
    );

    FastqReaderConfig config;
    config.min_length = 4;

    auto reader = FastqReader::Open(path, config);
    ASSERT_NE(reader, nullptr);

    auto records = reader->ReadAll();
    ASSERT_EQ(records.size(), 2);
    EXPECT_EQ(records[0].id, "read2");
    EXPECT_EQ(records[1].id, "read3");
}

TEST_F(FastqReaderTest, MaxLengthFilter) {
    auto path = CreateTestFastq("varied.fastq",
        "@read1\nAC\n+\nII\n"
        "@read2\nACGT\n+\nIIII\n"
        "@read3\nACGTACGT\n+\nIIIIIIII\n"
    );

    FastqReaderConfig config;
    config.max_length = 4;

    auto reader = FastqReader::Open(path, config);
    ASSERT_NE(reader, nullptr);

    auto records = reader->ReadAll();
    ASSERT_EQ(records.size(), 2);
    EXPECT_EQ(records[0].id, "read1");
    EXPECT_EQ(records[1].id, "read2");
}

// ========== Statistics ==========

TEST_F(FastqReaderTest, StatsAfterRead) {
    auto path = CreateTestFastq("stats.fastq",
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nTGCAACGT\n+\nHHHHHHHH\n"
        "@read3\nGG\n+\nFF\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    reader->ReadAll();

    auto stats = reader->GetStats();
    EXPECT_EQ(stats.total_records, 3);
    EXPECT_EQ(stats.total_bases, 14);  // 4 + 8 + 2
    EXPECT_EQ(stats.min_length, 2);
    EXPECT_EQ(stats.max_length, 8);
    EXPECT_NEAR(stats.avg_length, 14.0f / 3.0f, 0.01f);
}

// ========== Edge Cases ==========

TEST_F(FastqReaderTest, WhitespaceHandling) {
    auto path = CreateTestFastq("whitespace.fastq",
        "@read1  \n"
        "ACGT  \n"
        "+\n"
        "IIII  \n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->sequence, "ACGT");
    EXPECT_EQ(record->quality, "IIII");
}

TEST_F(FastqReaderTest, ReadIdWithComment) {
    auto path = CreateTestFastq("comment.fastq",
        "@read1 some comment here\n"
        "ACGT\n"
        "+\n"
        "IIII\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    auto record = reader->Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->id, "read1 some comment here");
}

TEST_F(FastqReaderTest, HasMore) {
    auto path = CreateTestFastq("two.fastq",
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nTGCA\n+\nHHHH\n"
    );

    auto reader = FastqReader::Open(path);
    ASSERT_NE(reader, nullptr);

    EXPECT_TRUE(reader->HasMore());
    reader->Next();
    EXPECT_TRUE(reader->HasMore());
    reader->Next();
    EXPECT_FALSE(reader->HasMore());
}

// ========== Utility Functions ==========

TEST_F(FastqReaderTest, ValidateFastqRecord) {
    FastqRecord valid{"read1", "ACGT", "IIII"};
    EXPECT_TRUE(ValidateFastqRecord(valid));

    FastqRecord empty_id{"", "ACGT", "IIII"};
    EXPECT_FALSE(ValidateFastqRecord(empty_id));

    FastqRecord empty_seq{"read1", "", ""};
    EXPECT_FALSE(ValidateFastqRecord(empty_seq));

    FastqRecord length_mismatch{"read1", "ACGT", "III"};
    EXPECT_FALSE(ValidateFastqRecord(length_mismatch));

    FastqRecord invalid_char{"read1", "ACGX", "IIII"};
    EXPECT_FALSE(ValidateFastqRecord(invalid_char));
}

TEST_F(FastqReaderTest, IsGzipFileDetection) {
    auto plain = CreateTestFastq("plain.fastq", "@read1\nACGT\n+\nIIII\n");
    EXPECT_FALSE(IsGzipFile(plain));
}

TEST_F(FastqReaderTest, ReadFastqConvenience) {
    auto path = CreateTestFastq("conv.fastq",
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nTGCA\n+\nHHHH\n"
    );

    auto records = ReadFastq(path);
    EXPECT_EQ(records.size(), 2);

    auto limited = ReadFastq(path, 1);
    EXPECT_EQ(limited.size(), 1);
}

TEST_F(FastqReaderTest, CountFastqRecords) {
    auto path = CreateTestFastq("count.fastq",
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nTGCA\n+\nHHHH\n"
        "@read3\nGGCC\n+\nFFFF\n"
    );

    EXPECT_EQ(CountFastqRecords(path), 3);
}

}  // namespace
}  // namespace llmap::io
