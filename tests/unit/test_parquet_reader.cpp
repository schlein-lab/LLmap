// LLmap — Unit tests for Parquet probabilistic reader and round-trip validation.

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "core/alignment_record.h"
#include "output/parquet_reader.h"
#include "output/parquet_writer.h"

namespace llmap::output {
namespace {

// ========== ParseCSVLine tests ==========

TEST(ParquetReaderTest, ParseCSVLineValid) {
    auto entry = ParseCSVLine("read001,chr14,0.95,0.85,2,15,1");
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->read_id, "read001");
    EXPECT_EQ(entry->bucket_id, "chr14");
    EXPECT_FLOAT_EQ(entry->probability, 0.95f);
    EXPECT_FLOAT_EQ(entry->confidence, 0.85f);
    EXPECT_EQ(entry->level, 2);
    EXPECT_EQ(entry->iteration, 15u);
    EXPECT_TRUE(entry->is_collapsed);
}

TEST(ParquetReaderTest, ParseCSVLineNotCollapsed) {
    auto entry = ParseCSVLine("read002,chr22,0.5,0.4,1,10,0");
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->read_id, "read002");
    EXPECT_EQ(entry->bucket_id, "chr22");
    EXPECT_FALSE(entry->is_collapsed);
}

TEST(ParquetReaderTest, ParseCSVLineEmpty) {
    auto entry = ParseCSVLine("");
    EXPECT_FALSE(entry.has_value());
}

TEST(ParquetReaderTest, ParseCSVLineTooFewFields) {
    auto entry = ParseCSVLine("read001,chr14,0.95");
    EXPECT_FALSE(entry.has_value());
}

TEST(ParquetReaderTest, ParseCSVLineWithWhitespace) {
    auto entry = ParseCSVLine("  read001,chr14,0.95,0.85,2,15,1  ");
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->read_id, "read001");
}

TEST(ParquetReaderTest, ParseCSVLineUnmapped) {
    auto entry = ParseCSVLine("unmapped001,*,0.0,0.0,0,0,0");
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->bucket_id, "*");
    EXPECT_FLOAT_EQ(entry->probability, 0.0f);
}

// ========== ParquetReader open tests ==========

TEST(ParquetReaderTest, OpenNonexistent) {
    auto reader = ParquetReader::Open("/nonexistent/file.parquet");
    EXPECT_EQ(reader, nullptr);
}

TEST(ParquetReaderTest, OpenEmptyCSV) {
    auto tmp = std::filesystem::temp_directory_path() / "test_empty.csv";
    std::ofstream f(tmp);
    f.close();

    auto reader = ParquetReader::Open(tmp);
    EXPECT_EQ(reader, nullptr);

    std::filesystem::remove(tmp);
}

TEST(ParquetReaderTest, OpenValidCSV) {
    auto tmp = std::filesystem::temp_directory_path() / "test_valid.csv";
    std::ofstream f(tmp);
    f << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
    f << "read001,chr14,0.95,0.85,2,15,1\n";
    f.close();

    auto reader = ParquetReader::Open(tmp);
    ASSERT_NE(reader, nullptr);
    EXPECT_EQ(reader->DetectedFormat(), ParquetOutputFormat::CSV);

    std::filesystem::remove(tmp);
}

// ========== ReadAll tests ==========

TEST(ParquetReaderTest, ReadAllFromCSV) {
    auto tmp = std::filesystem::temp_directory_path() / "test_readall.csv";
    std::ofstream f(tmp);
    f << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
    f << "read001,chr14,0.95,0.85,2,15,1\n";
    f << "read002,chr22,0.80,0.70,1,10,0\n";
    f << "read003,*,0.0,0.0,0,0,0\n";
    f.close();

    auto reader = ParquetReader::Open(tmp);
    ASSERT_NE(reader, nullptr);

    auto entries = reader->ReadAll();
    ASSERT_EQ(entries.size(), 3u);

    EXPECT_EQ(entries[0].read_id, "read001");
    EXPECT_EQ(entries[0].bucket_id, "chr14");
    EXPECT_TRUE(entries[0].is_collapsed);

    EXPECT_EQ(entries[1].read_id, "read002");
    EXPECT_FALSE(entries[1].is_collapsed);

    EXPECT_EQ(entries[2].read_id, "read003");
    EXPECT_EQ(entries[2].bucket_id, "*");

    auto stats = reader->GetStats();
    EXPECT_EQ(stats.entries_read, 3u);
    EXPECT_EQ(stats.unique_reads, 3u);
    EXPECT_EQ(stats.collapsed_entries, 1u);
    EXPECT_EQ(stats.unmapped_entries, 1u);

    std::filesystem::remove(tmp);
}

TEST(ParquetReaderTest, ReadAllWithMinProbabilityFilter) {
    auto tmp = std::filesystem::temp_directory_path() / "test_minprob.csv";
    std::ofstream f(tmp);
    f << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
    f << "read001,chr14,0.95,0.85,2,15,1\n";
    f << "read002,chr22,0.10,0.70,1,10,0\n";
    f.close();

    ParquetReaderConfig config;
    config.min_probability = 0.5f;

    auto reader = ParquetReader::Open(tmp, config);
    ASSERT_NE(reader, nullptr);

    auto entries = reader->ReadAll();
    ASSERT_EQ(entries.size(), 1u);
    EXPECT_EQ(entries[0].read_id, "read001");

    std::filesystem::remove(tmp);
}

TEST(ParquetReaderTest, ReadAllSkipUnmapped) {
    auto tmp = std::filesystem::temp_directory_path() / "test_skipunmapped.csv";
    std::ofstream f(tmp);
    f << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
    f << "read001,chr14,0.95,0.85,2,15,1\n";
    f << "read002,*,0.0,0.0,0,0,0\n";
    f.close();

    ParquetReaderConfig config;
    config.skip_unmapped = true;

    auto reader = ParquetReader::Open(tmp, config);
    ASSERT_NE(reader, nullptr);

    auto entries = reader->ReadAll();
    ASSERT_EQ(entries.size(), 1u);
    EXPECT_EQ(entries[0].read_id, "read001");

    std::filesystem::remove(tmp);
}

// ========== ReadBatch tests ==========

TEST(ParquetReaderTest, ReadBatch) {
    auto tmp = std::filesystem::temp_directory_path() / "test_batch.csv";
    std::ofstream f(tmp);
    f << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
    for (int i = 0; i < 10; ++i) {
        f << "read" << i << ",chr14,0.95,0.85,2,15,1\n";
    }
    f.close();

    auto reader = ParquetReader::Open(tmp);
    ASSERT_NE(reader, nullptr);

    auto batch1 = reader->ReadBatch(3);
    EXPECT_EQ(batch1.size(), 3u);
    EXPECT_TRUE(reader->HasMore());

    auto batch2 = reader->ReadBatch(3);
    EXPECT_EQ(batch2.size(), 3u);
    EXPECT_TRUE(reader->HasMore());

    auto batch3 = reader->ReadBatch(10);
    EXPECT_EQ(batch3.size(), 4u);
    EXPECT_FALSE(reader->HasMore());

    std::filesystem::remove(tmp);
}

// ========== GroupByReadId tests ==========

TEST(ParquetReaderTest, GroupByReadId) {
    std::vector<ProbabilityEntry> entries = {
        {"read001", "chr14", 0.8f, 0.9f, 2, 15, true},
        {"read001", "chr22", 0.2f, 0.5f, 2, 15, false},
        {"read002", "chr1", 1.0f, 0.95f, 3, 20, true},
        {"read001", "chr1", 0.1f, 0.3f, 2, 15, false},
    };

    auto groups = GroupByReadId(entries);
    EXPECT_EQ(groups.size(), 2u);

    std::size_t read001_count = 0;
    std::size_t read002_count = 0;
    for (const auto& group : groups) {
        if (!group.empty() && group[0].read_id == "read001") {
            read001_count = group.size();
        }
        if (!group.empty() && group[0].read_id == "read002") {
            read002_count = group.size();
        }
    }
    EXPECT_EQ(read001_count, 3u);
    EXPECT_EQ(read002_count, 1u);
}

// ========== ValidateRoundTrip tests ==========

TEST(ParquetReaderTest, ValidateRoundTripSuccess) {
    std::vector<ProbabilityEntry> entries = {
        {"read001", "chr14", 0.95f, 0.85f, 2, 15, true},
        {"read002", "chr22", 0.80f, 0.70f, 1, 10, false},
    };

    auto result = ValidateRoundTrip(entries, entries);
    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.mismatches, 0u);
    EXPECT_EQ(result.entries_original, 2u);
    EXPECT_EQ(result.entries_reread, 2u);
}

TEST(ParquetReaderTest, ValidateRoundTripCountMismatch) {
    std::vector<ProbabilityEntry> original = {
        {"read001", "chr14", 0.95f, 0.85f, 2, 15, true},
    };
    std::vector<ProbabilityEntry> reread = {
        {"read001", "chr14", 0.95f, 0.85f, 2, 15, true},
        {"read002", "chr22", 0.80f, 0.70f, 1, 10, false},
    };

    auto result = ValidateRoundTrip(original, reread);
    EXPECT_FALSE(result.success);
    EXPECT_NE(result.error.find("mismatch"), std::string::npos);
}

TEST(ParquetReaderTest, ValidateRoundTripValueMismatch) {
    std::vector<ProbabilityEntry> original = {
        {"read001", "chr14", 0.95f, 0.85f, 2, 15, true},
    };
    std::vector<ProbabilityEntry> reread = {
        {"read001", "chr14", 0.50f, 0.85f, 2, 15, true},  // Different probability
    };

    auto result = ValidateRoundTrip(original, reread);
    EXPECT_FALSE(result.success);
    EXPECT_EQ(result.mismatches, 1u);
}

TEST(ParquetReaderTest, ValidateRoundTripTolerance) {
    std::vector<ProbabilityEntry> original = {
        {"read001", "chr14", 0.950000f, 0.85f, 2, 15, true},
    };
    std::vector<ProbabilityEntry> reread = {
        {"read001", "chr14", 0.950001f, 0.85f, 2, 15, true},  // Within tolerance
    };

    auto result = ValidateRoundTrip(original, reread, 1e-4f);
    EXPECT_TRUE(result.success);
}

// ========== Full round-trip tests ==========

TEST(ParquetReaderTest, FullRoundTripCSV) {
    auto tmp = std::filesystem::temp_directory_path() / "test_roundtrip.parquet";

    // Generate some alignment records
    std::vector<AlignmentRecord> records;

    AlignmentHit hit1;
    hit1.target_id = "chr14";
    hit1.start = 100000;
    hit1.end = 101000;
    auto rec1 = make_mapped("roundtrip001", 1000, hit1);
    rec1.collapsed_at_level = 2;
    rec1.collapsed_at_iteration = 15;
    rec1.confidence_scores = {0.95f};
    records.push_back(rec1);

    AlignmentHit hit2;
    hit2.target_id = "chr22";
    hit2.start = 200000;
    hit2.end = 201000;
    auto rec2 = make_mapped("roundtrip002", 1000, hit2);
    rec2.collapsed_at_level = 3;
    rec2.collapsed_at_iteration = 20;
    rec2.confidence_scores = {0.85f};
    records.push_back(rec2);

    records.push_back(make_unmapped("roundtrip003", 500, RejectionReason::NoSeeds));

    // Write
    ParquetWriterConfig write_config;
    write_config.format = ParquetOutputFormat::CSV;
    auto writer = ParquetWriter::Create(tmp, write_config);
    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->WriteBatch(records));
    writer->Close();

    // Read back
    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto reader = ParquetReader::Open(csv_path);
    ASSERT_NE(reader, nullptr);
    auto reread_entries = reader->ReadAll();

    // Also convert original records to entries for comparison
    std::vector<ProbabilityEntry> original_entries;
    for (const auto& rec : records) {
        auto entries = RecordToEntries(rec);
        for (auto& e : entries) {
            original_entries.push_back(std::move(e));
        }
    }

    // Validate round-trip
    auto result = ValidateRoundTrip(original_entries, reread_entries);
    EXPECT_TRUE(result.success) << result.error;
    EXPECT_EQ(result.entries_original, result.entries_reread);

    std::filesystem::remove(csv_path);
}

TEST(ParquetReaderTest, FullRoundTripWithAlternatives) {
    auto tmp = std::filesystem::temp_directory_path() / "test_rt_alts.parquet";

    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;

    AlignmentHit alt1;
    alt1.target_id = "chr22";
    alt1.start = 200000;
    alt1.end = 201000;

    auto rec = make_mapped("rtalt001", 1000, primary, {alt1});
    rec.collapsed_at_level = 2;
    rec.collapsed_at_iteration = 15;
    rec.confidence_scores = {0.90f};

    ParquetWriterConfig write_config;
    write_config.format = ParquetOutputFormat::CSV;
    auto writer = ParquetWriter::Create(tmp, write_config);
    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto reader = ParquetReader::Open(csv_path);
    ASSERT_NE(reader, nullptr);
    auto reread = reader->ReadAll();

    EXPECT_EQ(reread.size(), 2u);  // Primary + 1 alternative

    auto stats = reader->GetStats();
    EXPECT_EQ(stats.entries_read, 2u);
    EXPECT_EQ(stats.collapsed_entries, 1u);

    std::filesystem::remove(csv_path);
}

TEST(ParquetReaderTest, FullRoundTripTentative) {
    auto tmp = std::filesystem::temp_directory_path() / "test_rt_tent.parquet";

    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.45f},
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.35f},
    };

    auto rec = make_tentative("rttent001", 1000, targets,
                              RejectionReason::DidNotConverge);
    rec.collapsed_at_level = 1;
    rec.collapsed_at_iteration = 10;

    ParquetWriterConfig write_config;
    write_config.format = ParquetOutputFormat::CSV;
    auto writer = ParquetWriter::Create(tmp, write_config);
    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto reader = ParquetReader::Open(csv_path);
    ASSERT_NE(reader, nullptr);
    auto reread = reader->ReadAll();

    EXPECT_EQ(reread.size(), 2u);

    bool found_chr14 = false;
    bool found_chr22 = false;
    for (const auto& e : reread) {
        if (e.bucket_id == "chr14") {
            found_chr14 = true;
            EXPECT_FLOAT_EQ(e.probability, 0.45f);
        }
        if (e.bucket_id == "chr22") {
            found_chr22 = true;
            EXPECT_FLOAT_EQ(e.probability, 0.35f);
        }
    }
    EXPECT_TRUE(found_chr14);
    EXPECT_TRUE(found_chr22);

    std::filesystem::remove(csv_path);
}

// ========== ReadParquet convenience function ==========

TEST(ParquetReaderTest, ReadParquetConvenience) {
    auto tmp = std::filesystem::temp_directory_path() / "test_convenience.csv";
    std::ofstream f(tmp);
    f << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
    f << "conv001,chr14,0.95,0.85,2,15,1\n";
    f << "conv002,chr22,0.80,0.70,1,10,0\n";
    f.close();

    auto entries = ReadParquet(tmp);
    EXPECT_EQ(entries.size(), 2u);
    EXPECT_EQ(entries[0].read_id, "conv001");
    EXPECT_EQ(entries[1].read_id, "conv002");

    std::filesystem::remove(tmp);
}

// ========== Lossless invariant test ==========

TEST(ParquetReaderTest, LosslessInvariant_RoundTrip) {
    auto tmp = std::filesystem::temp_directory_path() / "test_lossless_rt.parquet";

    const std::size_t N = 100;
    std::vector<AlignmentRecord> records;
    records.reserve(N);

    for (std::size_t i = 0; i < N; ++i) {
        if (i % 3 == 0) {
            AlignmentHit hit;
            hit.target_id = "chr14";
            hit.start = 100000 + static_cast<std::uint32_t>(i * 1000);
            hit.end = hit.start + 1000;
            auto rec = make_mapped("lossless_" + std::to_string(i), 1000, hit);
            rec.collapsed_at_level = 2;
            rec.collapsed_at_iteration = 15;
            rec.confidence_scores = {0.95f};
            records.push_back(rec);
        } else if (i % 3 == 1) {
            records.push_back(make_unmapped("lossless_" + std::to_string(i), 500,
                                            RejectionReason::NoSeeds));
        } else {
            std::vector<TentativeTarget> targets = {
                {"chr14", 100000, 101000, 50, 800, 0.85f, 0.6f}
            };
            auto rec = make_tentative("lossless_" + std::to_string(i), 800,
                                      targets, RejectionReason::DidNotConverge);
            rec.collapsed_at_level = 1;
            rec.collapsed_at_iteration = 10;
            records.push_back(rec);
        }
    }

    // Write
    ParquetWriterConfig write_config;
    write_config.format = ParquetOutputFormat::CSV;
    write_config.include_unmapped = true;
    write_config.include_tentative = true;
    auto writer = ParquetWriter::Create(tmp, write_config);
    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->WriteBatch(records));

    auto write_stats = writer->GetStats();
    EXPECT_EQ(write_stats.records_written, N);

    writer->Close();

    // Read back
    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto reader = ParquetReader::Open(csv_path);
    ASSERT_NE(reader, nullptr);
    auto reread = reader->ReadAll();

    // Convert original for comparison
    std::vector<ProbabilityEntry> original;
    for (const auto& rec : records) {
        auto entries = RecordToEntries(rec);
        for (auto& e : entries) {
            original.push_back(std::move(e));
        }
    }

    // Validate: all N records produce exactly N entries (one per record min)
    EXPECT_EQ(reread.size(), original.size());

    auto result = ValidateRoundTrip(original, reread);
    EXPECT_TRUE(result.success) << result.error;

    std::filesystem::remove(csv_path);
}

// ========== Large-scale lossless invariant ==========

TEST(ParquetReaderTest, LosslessInvariant_10K) {
    auto tmp = std::filesystem::temp_directory_path() / "test_lossless_10k.parquet";

    const std::size_t N = 10000;
    std::vector<ProbabilityEntry> original;
    original.reserve(N);

    for (std::size_t i = 0; i < N; ++i) {
        ProbabilityEntry entry{
            .read_id = "read_" + std::to_string(i),
            .bucket_id = (i % 10 == 0) ? "*" : "chr" + std::to_string(i % 22 + 1),
            .probability = static_cast<float>(i % 100) / 100.0f,
            .confidence = 0.9f,
            .level = static_cast<std::uint8_t>(i % 4),
            .iteration = static_cast<std::uint32_t>(i % 50),
            .is_collapsed = (i % 3 == 0),
        };
        original.push_back(entry);
    }

    // Write directly as entries
    ParquetWriterConfig write_config;
    write_config.format = ParquetOutputFormat::CSV;
    auto writer = ParquetWriter::Create(tmp, write_config);
    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->WriteEntries(original));
    writer->Close();

    // Read back
    auto csv_path = tmp;
    csv_path.replace_extension(".csv");
    auto entries = ReadParquet(csv_path);

    EXPECT_EQ(entries.size(), N);

    auto result = ValidateRoundTrip(original, entries);
    EXPECT_TRUE(result.success) << result.error;

    std::filesystem::remove(csv_path);
}

}  // namespace
}  // namespace llmap::output
