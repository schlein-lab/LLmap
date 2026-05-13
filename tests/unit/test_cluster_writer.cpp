// Unit tests for ClusterWriter

#include <gtest/gtest.h>

#include "output/cluster_writer.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace llmap::output {
namespace {

class ClusterWriterTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_writer_test";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::string ReadFile(const std::filesystem::path& path) {
        std::ifstream in(path);
        std::stringstream buffer;
        buffer << in.rdbuf();
        return buffer.str();
    }

    std::vector<std::string> ReadLines(const std::filesystem::path& path) {
        std::vector<std::string> lines;
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line)) {
            lines.push_back(line);
        }
        return lines;
    }

    std::filesystem::path test_dir_;
};

// ========== Basic Functionality ==========

TEST_F(ClusterWriterTest, CreateWriter) {
    auto path = test_dir_ / "out.tsv";
    auto writer = ClusterWriter::Create(path);
    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(std::filesystem::exists(path));
}

TEST_F(ClusterWriterTest, WriteHeader) {
    auto path = test_dir_ / "header.tsv";

    ClusterWriterConfig config;
    config.include_header = true;

    auto writer = ClusterWriter::Create(path, config);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    auto content = ReadFile(path);
    EXPECT_TRUE(content.find("read_id") != std::string::npos);
    EXPECT_TRUE(content.find("cluster_id") != std::string::npos);
    EXPECT_TRUE(content.find("confidence") != std::string::npos);
}

TEST_F(ClusterWriterTest, WriteNoHeader) {
    auto path = test_dir_ / "no_header.tsv";

    ClusterWriterConfig config;
    config.include_header = false;

    auto writer = ClusterWriter::Create(path, config);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    auto content = ReadFile(path);
    EXPECT_TRUE(content.empty());
}

TEST_F(ClusterWriterTest, WriteSingleAssignment) {
    auto path = test_dir_ / "single.tsv";
    auto writer = ClusterWriter::Create(path);
    ASSERT_NE(writer, nullptr);

    ClusterAssignment assignment;
    assignment.read_id = "read1";
    assignment.cluster_id = 5;
    assignment.confidence = 0.95f;
    assignment.is_representative = true;
    assignment.cluster_size = 10;
    assignment.read_length = 1000;

    EXPECT_TRUE(writer->Write(assignment));
    writer->Close();

    auto lines = ReadLines(path);
    ASSERT_GE(lines.size(), 2);  // header + 1 record

    EXPECT_TRUE(lines[1].find("read1") != std::string::npos);
    EXPECT_TRUE(lines[1].find("5") != std::string::npos);
    EXPECT_TRUE(lines[1].find("true") != std::string::npos);
}

TEST_F(ClusterWriterTest, WriteMultipleAssignments) {
    auto path = test_dir_ / "multi.tsv";
    auto writer = ClusterWriter::Create(path);
    ASSERT_NE(writer, nullptr);

    for (int i = 0; i < 10; ++i) {
        ClusterAssignment assignment;
        assignment.read_id = "read" + std::to_string(i);
        assignment.cluster_id = i % 3;
        assignment.confidence = 0.8f + i * 0.01f;
        assignment.is_representative = (i % 3 == 0);
        assignment.cluster_size = 100;
        assignment.read_length = 500 + i * 10;

        EXPECT_TRUE(writer->Write(assignment));
    }

    writer->Close();

    auto lines = ReadLines(path);
    EXPECT_EQ(lines.size(), 11);  // header + 10 records
}

TEST_F(ClusterWriterTest, WriteBatch) {
    auto path = test_dir_ / "batch.tsv";
    auto writer = ClusterWriter::Create(path);
    ASSERT_NE(writer, nullptr);

    std::vector<ClusterAssignment> assignments;
    for (int i = 0; i < 5; ++i) {
        ClusterAssignment a;
        a.read_id = "read" + std::to_string(i);
        a.cluster_id = i;
        a.confidence = 0.9f;
        a.is_representative = false;
        a.cluster_size = 50;
        a.read_length = 200;
        assignments.push_back(a);
    }

    EXPECT_TRUE(writer->WriteBatch(assignments));
    writer->Close();

    auto lines = ReadLines(path);
    EXPECT_EQ(lines.size(), 6);  // header + 5 records
}

// ========== Statistics ==========

TEST_F(ClusterWriterTest, StatsTracking) {
    auto path = test_dir_ / "stats.tsv";
    auto writer = ClusterWriter::Create(path);
    ASSERT_NE(writer, nullptr);

    // Write 5 records across 3 clusters, 2 representatives
    std::vector<std::tuple<uint32_t, bool>> data = {
        {0, true}, {0, false}, {1, true}, {1, false}, {2, false}
    };

    for (size_t i = 0; i < data.size(); ++i) {
        ClusterAssignment a;
        a.read_id = "read" + std::to_string(i);
        a.cluster_id = std::get<0>(data[i]);
        a.confidence = 0.9f;
        a.is_representative = std::get<1>(data[i]);
        a.cluster_size = 100;
        a.read_length = 500;
        writer->Write(a);
    }

    writer->Close();

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.records_written, 5);
    EXPECT_EQ(stats.clusters_written, 3);
    EXPECT_EQ(stats.representatives_written, 2);
}

// ========== Configuration ==========

TEST_F(ClusterWriterTest, CustomDelimiter) {
    auto path = test_dir_ / "comma.csv";

    ClusterWriterConfig config;
    config.delimiter = ",";

    auto writer = ClusterWriter::Create(path, config);
    ASSERT_NE(writer, nullptr);

    ClusterAssignment a;
    a.read_id = "read1";
    a.cluster_id = 0;
    a.confidence = 0.9f;
    a.is_representative = false;
    a.cluster_size = 10;
    a.read_length = 100;

    writer->Write(a);
    writer->Close();

    auto content = ReadFile(path);
    EXPECT_TRUE(content.find(",") != std::string::npos);
    EXPECT_TRUE(content.find("\t") == std::string::npos);
}

TEST_F(ClusterWriterTest, ExcludeReadLength) {
    auto path = test_dir_ / "no_len.tsv";

    ClusterWriterConfig config;
    config.include_sequence_length = false;

    auto writer = ClusterWriter::Create(path, config);
    ASSERT_NE(writer, nullptr);

    ClusterAssignment a;
    a.read_id = "read1";
    a.cluster_id = 0;
    a.confidence = 0.9f;
    a.is_representative = false;
    a.cluster_size = 10;
    a.read_length = 12345;

    writer->Write(a);
    writer->Close();

    auto content = ReadFile(path);
    EXPECT_TRUE(content.find("12345") == std::string::npos);
}

// ========== Error Handling ==========

TEST_F(ClusterWriterTest, ParquetNotSupported) {
    auto path = test_dir_ / "out.parquet";

    ClusterWriterConfig config;
    config.format = ClusterOutputFormat::Parquet;

    auto writer = ClusterWriter::Create(path, config);
    EXPECT_EQ(writer, nullptr);
}

// ========== Convenience Function ==========

TEST_F(ClusterWriterTest, WriteClusterAssignmentsConvenience) {
    auto path = test_dir_ / "conv.tsv";

    std::vector<ClusterAssignment> assignments;
    for (int i = 0; i < 3; ++i) {
        ClusterAssignment a;
        a.read_id = "read" + std::to_string(i);
        a.cluster_id = i;
        a.confidence = 0.95f;
        a.is_representative = (i == 0);
        a.cluster_size = 20;
        a.read_length = 300;
        assignments.push_back(a);
    }

    EXPECT_TRUE(WriteClusterAssignments(path, assignments));

    auto lines = ReadLines(path);
    EXPECT_EQ(lines.size(), 4);  // header + 3 records
}

}  // namespace
}  // namespace llmap::output
