// LLmap — Unit tests for Parquet probabilistic output writer.

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/alignment_record.h"
#include "output/parquet_writer.h"

namespace llmap::output {
namespace {

// Helper to read CSV file contents
std::string ReadFile(const std::filesystem::path& path) {
    std::ifstream f(path);
    std::stringstream buf;
    buf << f.rdbuf();
    return buf.str();
}

// Helper to count lines in a string (excluding header)
std::size_t CountDataLines(const std::string& content) {
    std::size_t count = 0;
    std::istringstream iss(content);
    std::string line;
    bool header_skipped = false;
    while (std::getline(iss, line)) {
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }
        if (!line.empty()) {
            ++count;
        }
    }
    return count;
}

// ========== ParquetWriter creation tests ==========

TEST(ParquetWriterTest, CreateWriter) {
    auto tmp = std::filesystem::temp_directory_path() / "test_create.parquet";
    auto writer = ParquetWriter::Create(tmp);

    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->Close());

    // Clean up (may be .csv if Arrow not available)
    std::filesystem::remove(tmp);
    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    std::filesystem::remove(csv_path);
}

TEST(ParquetWriterTest, CreateWithCustomConfig) {
    auto tmp = std::filesystem::temp_directory_path() / "test_custom.parquet";

    ParquetWriterConfig config;
    config.sample_name = "testsample";
    config.compression = "gzip";
    config.min_probability = 0.1f;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    std::filesystem::remove(tmp);
    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    std::filesystem::remove(csv_path);
}

TEST(ParquetWriterTest, ArrowAvailable) {
#ifdef LLMAP_HAS_ARROW
    EXPECT_TRUE(ParquetWriter::ArrowAvailable());
#else
    EXPECT_FALSE(ParquetWriter::ArrowAvailable());
#endif
}

// ========== CSV fallback tests ==========

TEST(ParquetWriterTest, WritesToCsvWhenArrowUnavailable) {
    auto tmp = std::filesystem::temp_directory_path() / "test_csv_fallback.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;  // Force CSV

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;
    primary.cigar.ops = "1000M";
    primary.score = 950;
    primary.nm = 5;

    auto rec = make_mapped("read001", 1000, primary);
    rec.collapsed_at_level = 2;
    rec.collapsed_at_iteration = 15;
    rec.confidence_scores = {0.95f};

    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");

    auto content = ReadFile(csv_path);

    // Check header
    EXPECT_NE(content.find("read_id,bucket_id,probability"), std::string::npos);

    // Check data
    EXPECT_NE(content.find("read001"), std::string::npos);
    EXPECT_NE(content.find("chr14"), std::string::npos);

    std::filesystem::remove(csv_path);
}

// ========== Record to entries conversion tests ==========

TEST(ParquetWriterTest, RecordToEntriesMapped) {
    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;

    auto rec = make_mapped("read001", 1000, primary);
    rec.collapsed_at_level = 2;
    rec.collapsed_at_iteration = 15;
    rec.confidence_scores = {0.95f};

    auto entries = RecordToEntries(rec);

    ASSERT_EQ(entries.size(), 1u);
    EXPECT_EQ(entries[0].read_id, "read001");
    EXPECT_EQ(entries[0].bucket_id, "chr14");
    EXPECT_FLOAT_EQ(entries[0].probability, 1.0f);
    EXPECT_FLOAT_EQ(entries[0].confidence, 0.95f);
    EXPECT_EQ(entries[0].level, 2);
    EXPECT_EQ(entries[0].iteration, 15u);
    EXPECT_TRUE(entries[0].is_collapsed);
}

TEST(ParquetWriterTest, RecordToEntriesMappedWithAlternatives) {
    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;

    AlignmentHit alt1;
    alt1.target_id = "chr22";
    alt1.start = 200000;
    alt1.end = 201000;

    AlignmentHit alt2;
    alt2.target_id = "chr1";
    alt2.start = 300000;
    alt2.end = 301000;

    auto rec = make_mapped("read002", 1000, primary, {alt1, alt2});
    rec.confidence_scores = {0.8f};

    auto entries = RecordToEntries(rec);

    ASSERT_EQ(entries.size(), 3u);

    // Primary
    EXPECT_EQ(entries[0].bucket_id, "chr14");
    EXPECT_FLOAT_EQ(entries[0].probability, 1.0f);
    EXPECT_TRUE(entries[0].is_collapsed);

    // Alternatives (with decreasing probability)
    EXPECT_EQ(entries[1].bucket_id, "chr22");
    EXPECT_LT(entries[1].probability, 1.0f);
    EXPECT_FALSE(entries[1].is_collapsed);

    EXPECT_EQ(entries[2].bucket_id, "chr1");
    EXPECT_LT(entries[2].probability, entries[1].probability);
    EXPECT_FALSE(entries[2].is_collapsed);
}

TEST(ParquetWriterTest, RecordToEntriesTentative) {
    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.45f},
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.35f}
    };

    auto rec = make_tentative("read003", 1000, targets,
                              RejectionReason::DidNotConverge);
    rec.collapsed_at_level = 1;
    rec.collapsed_at_iteration = 10;

    auto entries = RecordToEntries(rec);

    ASSERT_EQ(entries.size(), 2u);

    EXPECT_EQ(entries[0].read_id, "read003");
    EXPECT_EQ(entries[0].bucket_id, "chr14");
    EXPECT_FLOAT_EQ(entries[0].probability, 0.45f);
    EXPECT_FLOAT_EQ(entries[0].confidence, 0.85f);
    EXPECT_FALSE(entries[0].is_collapsed);

    EXPECT_EQ(entries[1].bucket_id, "chr22");
    EXPECT_FLOAT_EQ(entries[1].probability, 0.35f);
}

TEST(ParquetWriterTest, RecordToEntriesUnmapped) {
    auto rec = make_unmapped("read004", 500, RejectionReason::NoSeeds);

    auto entries = RecordToEntries(rec);

    ASSERT_EQ(entries.size(), 1u);
    EXPECT_EQ(entries[0].read_id, "read004");
    EXPECT_EQ(entries[0].bucket_id, "*");
    EXPECT_FLOAT_EQ(entries[0].probability, 0.0f);
    EXPECT_FALSE(entries[0].is_collapsed);
}

TEST(ParquetWriterTest, RecordToEntriesFilterByMinProbability) {
    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.45f},
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.05f}  // Low probability
    };

    auto rec = make_tentative("read005", 1000, targets,
                              RejectionReason::DidNotConverge);

    auto entries = RecordToEntries(rec, 0.1f);

    ASSERT_EQ(entries.size(), 1u);  // Only chr14 passes threshold
    EXPECT_EQ(entries[0].bucket_id, "chr14");
}

// ========== Write record tests ==========

TEST(ParquetWriterTest, WriteMappedRecord) {
    auto tmp = std::filesystem::temp_directory_path() / "test_mapped.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;
    primary.cigar.ops = "1000M";
    primary.score = 950;

    auto rec = make_mapped("mapped001", 1000, primary);
    rec.collapsed_at_level = 2;
    rec.collapsed_at_iteration = 15;
    rec.confidence_scores = {0.95f};

    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_NE(content.find("mapped001"), std::string::npos);
    EXPECT_NE(content.find("chr14"), std::string::npos);
    EXPECT_NE(content.find("1.000000"), std::string::npos);  // probability
    EXPECT_NE(content.find(",2,"), std::string::npos);       // level
    EXPECT_NE(content.find(",15,"), std::string::npos);      // iteration
    EXPECT_NE(content.find(",1"), std::string::npos);        // is_collapsed

    std::filesystem::remove(csv_path);
}

TEST(ParquetWriterTest, WriteUnmappedRecord) {
    auto tmp = std::filesystem::temp_directory_path() / "test_unmapped.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    auto rec = make_unmapped("unmapped001", 500, RejectionReason::NoSeeds);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_NE(content.find("unmapped001"), std::string::npos);
    EXPECT_NE(content.find("*"), std::string::npos);  // No bucket
    EXPECT_NE(content.find("0.000000"), std::string::npos);  // probability 0

    std::filesystem::remove(csv_path);
}

TEST(ParquetWriterTest, WriteTentativeRecord) {
    auto tmp = std::filesystem::temp_directory_path() / "test_tentative.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.45f},
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.35f}
    };

    auto rec = make_tentative("tentative001", 1000, targets,
                              RejectionReason::DidNotConverge);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_NE(content.find("tentative001"), std::string::npos);
    EXPECT_NE(content.find("chr14"), std::string::npos);
    EXPECT_NE(content.find("chr22"), std::string::npos);
    EXPECT_NE(content.find("0.450000"), std::string::npos);
    EXPECT_NE(content.find("0.350000"), std::string::npos);

    EXPECT_EQ(CountDataLines(content), 2u);

    std::filesystem::remove(csv_path);
}

// ========== Filter tests ==========

TEST(ParquetWriterTest, SkipUnmappedWhenDisabled) {
    auto tmp = std::filesystem::temp_directory_path() / "test_skip_unmapped.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;
    config.include_unmapped = false;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    auto rec = make_unmapped("skip001", 500, RejectionReason::NoSeeds);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_EQ(content.find("skip001"), std::string::npos);

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.unmapped_records, 1u);
    EXPECT_EQ(stats.entries_written, 0u);

    std::filesystem::remove(csv_path);
}

TEST(ParquetWriterTest, SkipTentativeWhenDisabled) {
    auto tmp = std::filesystem::temp_directory_path() / "test_skip_tentative.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;
    config.include_tentative = false;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    auto rec = make_tentative("skip002", 1000, {},
                              RejectionReason::DidNotConverge);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_EQ(content.find("skip002"), std::string::npos);

    std::filesystem::remove(csv_path);
}

TEST(ParquetWriterTest, FilterByMinProbability) {
    auto tmp = std::filesystem::temp_directory_path() / "test_min_prob.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;
    config.min_probability = 0.4f;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.50f},  // Above threshold
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.30f}   // Below threshold
    };

    auto rec = make_tentative("filter001", 1000, targets,
                              RejectionReason::DidNotConverge);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_NE(content.find("chr14"), std::string::npos);
    EXPECT_EQ(content.find("chr22"), std::string::npos);  // Filtered out

    std::filesystem::remove(csv_path);
}

// ========== Batch write tests ==========

TEST(ParquetWriterTest, WriteBatch) {
    auto tmp = std::filesystem::temp_directory_path() / "test_batch.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    std::vector<AlignmentRecord> records;

    // Mapped
    AlignmentHit hit1;
    hit1.target_id = "chr14";
    hit1.start = 100000;
    hit1.end = 101000;
    hit1.cigar.ops = "1000M";
    hit1.score = 950;
    auto rec1 = make_mapped("batch001", 1000, hit1);
    rec1.confidence_scores = {0.9f};
    records.push_back(rec1);

    // Unmapped
    records.push_back(make_unmapped("batch002", 500, RejectionReason::NoSeeds));

    // Tentative with 2 targets
    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.50f},
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.30f}
    };
    records.push_back(make_tentative("batch003", 800, targets,
                                     RejectionReason::DidNotConverge));

    EXPECT_TRUE(writer->WriteBatch(records));
    writer->Close();

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.records_written, 3u);
    EXPECT_EQ(stats.mapped_records, 1u);
    EXPECT_EQ(stats.unmapped_records, 1u);
    EXPECT_EQ(stats.tentative_records, 1u);
    EXPECT_EQ(stats.entries_written, 4u);  // 1 + 1 + 2

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    std::filesystem::remove(csv_path);
}

// ========== Direct entry write tests ==========

TEST(ParquetWriterTest, WriteEntries) {
    auto tmp = std::filesystem::temp_directory_path() / "test_entries.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    std::vector<ProbabilityEntry> entries = {
        {"read001", "chr14", 0.8f, 0.9f, 2, 10, true},
        {"read001", "chr22", 0.2f, 0.5f, 2, 10, false},
        {"read002", "chr1", 1.0f, 0.95f, 3, 20, true}
    };

    EXPECT_TRUE(writer->WriteEntries(entries));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_EQ(CountDataLines(content), 3u);

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.entries_written, 3u);

    std::filesystem::remove(csv_path);
}

// ========== Statistics tests ==========

TEST(ParquetWriterTest, TrackStatistics) {
    auto tmp = std::filesystem::temp_directory_path() / "test_stats.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    hit.cigar.ops = "1000M";
    hit.score = 950;

    AlignmentHit alt;
    alt.target_id = "chr22";
    alt.start = 200000;
    alt.end = 201000;

    auto rec = make_mapped("stats001", 1000, hit, {alt});
    rec.confidence_scores = {0.9f};
    writer->Write(rec);

    auto rec2 = make_unmapped("stats002", 500, RejectionReason::NoSeeds);
    writer->Write(rec2);

    writer->Close();

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.records_written, 2u);
    EXPECT_EQ(stats.mapped_records, 1u);
    EXPECT_EQ(stats.unmapped_records, 1u);
    EXPECT_EQ(stats.entries_written, 3u);  // 1 primary + 1 alt + 1 unmapped
    EXPECT_GT(stats.write_time_ms, 0.0f);

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    std::filesystem::remove(csv_path);
}

// ========== Convenience function tests ==========

TEST(ParquetWriterTest, WriteParquetConvenience) {
    auto tmp = std::filesystem::temp_directory_path() / "test_convenience.parquet";

    std::vector<AlignmentRecord> records;

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    hit.cigar.ops = "1000M";
    hit.score = 950;
    records.push_back(make_mapped("conv001", 1000, hit));
    records.push_back(make_unmapped("conv002", 500, RejectionReason::NoSeeds));

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    EXPECT_TRUE(WriteParquet(tmp, records, config));

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    EXPECT_NE(content.find("conv001"), std::string::npos);
    EXPECT_NE(content.find("conv002"), std::string::npos);

    std::filesystem::remove(csv_path);
}

// ========== Error handling tests ==========

TEST(ParquetWriterTest, FailsOnInvalidPath) {
    auto writer = ParquetWriter::Create("/nonexistent/path/test.parquet");
    EXPECT_EQ(writer, nullptr);
}

TEST(ParquetWriterTest, ReturnsErrorAfterClose) {
    auto tmp = std::filesystem::temp_directory_path() / "test_closed.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    auto rec = make_mapped("fail001", 1000, hit);

    // After close, CSV file handle is invalid
    EXPECT_FALSE(writer->Write(rec));

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    std::filesystem::remove(csv_path);
}

// ========== Lossless invariant tests ==========

TEST(ParquetWriterTest, LosslessInvariant_NoSilentDrops) {
    auto tmp = std::filesystem::temp_directory_path() / "test_lossless.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;
    config.include_unmapped = true;
    config.include_tentative = true;
    config.min_probability = 0.0f;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    const std::size_t N = 100;
    std::vector<AlignmentRecord> records;
    records.reserve(N);

    for (std::size_t i = 0; i < N; ++i) {
        if (i % 3 == 0) {
            AlignmentHit hit;
            hit.target_id = "chr14";
            hit.start = 100000 + i * 1000;
            hit.end = hit.start + 1000;
            records.push_back(make_mapped("read_" + std::to_string(i), 1000, hit));
        } else if (i % 3 == 1) {
            records.push_back(make_unmapped("read_" + std::to_string(i), 500,
                                            RejectionReason::NoSeeds));
        } else {
            std::vector<TentativeTarget> targets = {
                {"chr14", 100000, 101000, 50, 800, 0.85f, 0.6f}
            };
            records.push_back(make_tentative("read_" + std::to_string(i), 800,
                                             targets, RejectionReason::DidNotConverge));
        }
    }

    EXPECT_TRUE(writer->WriteBatch(records));
    writer->Close();

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.records_written, N);
    EXPECT_EQ(stats.entries_written, N);  // Each record produces at least 1 entry

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);
    EXPECT_EQ(CountDataLines(content), N);

    std::filesystem::remove(csv_path);
}

// ========== Round-trip tests (schema validation) ==========

TEST(ParquetWriterTest, CSVSchemaCorrect) {
    auto tmp = std::filesystem::temp_directory_path() / "test_schema.parquet";

    ParquetWriterConfig config;
    config.format = ParquetOutputFormat::CSV;

    auto writer = ParquetWriter::Create(tmp, config);
    ASSERT_NE(writer, nullptr);

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    auto rec = make_mapped("schema001", 1000, hit);
    rec.collapsed_at_level = 2;
    rec.collapsed_at_iteration = 15;
    rec.confidence_scores = {0.95f};

    writer->Write(rec);
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto content = ReadFile(csv_path);

    // Verify header
    std::istringstream iss(content);
    std::string header;
    std::getline(iss, header);
    EXPECT_EQ(header, "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed");

    // Verify data line has 7 fields
    std::string data;
    std::getline(iss, data);
    std::size_t comma_count = 0;
    for (char c : data) {
        if (c == ',') ++comma_count;
    }
    EXPECT_EQ(comma_count, 6u);  // 7 fields = 6 commas

    std::filesystem::remove(csv_path);
}

}  // namespace
}  // namespace llmap::output
