// LLmap — Unit tests for BAM/SAM writer.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/alignment_record.h"
#include "output/bam_writer.h"

namespace llmap::output {
namespace {

// Helper to read file contents
std::string ReadFile(const std::filesystem::path& path) {
    std::ifstream f(path);
    std::stringstream buf;
    buf << f.rdbuf();
    return buf.str();
}

// Helper to create test references
std::vector<ReferenceSequence> TestRefs() {
    return {
        {"chr14", 107043718},
        {"chr22", 50818468},
        {"chrUn_GL000220v1", 161802}
    };
}

// ========== BamWriter creation tests ==========

TEST(BamWriterTest, CreateWithSamFormat) {
    auto tmp = std::filesystem::temp_directory_path() / "test_create.sam";
    auto refs = TestRefs();
    auto writer = BamWriter::Create(tmp, refs);

    ASSERT_NE(writer, nullptr);
    EXPECT_TRUE(writer->Close());

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, CreateWithCustomConfig) {
    auto tmp = std::filesystem::temp_directory_path() / "test_custom.sam";
    auto refs = TestRefs();

    BamWriterConfig config;
    config.sample_name = "testsample";
    config.read_group_id = "RG1";
    config.platform = "ILLUMINA";

    auto writer = BamWriter::Create(tmp, refs, config);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_NE(content.find("SM:testsample"), std::string::npos);
    EXPECT_NE(content.find("ID:RG1"), std::string::npos);
    EXPECT_NE(content.find("PL:ILLUMINA"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, HtslibAvailable) {
#ifdef LLMAP_HAS_HTSLIB
    EXPECT_TRUE(BamWriter::HtslibAvailable());
#else
    EXPECT_FALSE(BamWriter::HtslibAvailable());
#endif
}

// ========== SAM header tests ==========

TEST(BamWriterTest, WritesValidHeader) {
    auto tmp = std::filesystem::temp_directory_path() / "test_header.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    auto content = ReadFile(tmp);

    // Check header lines
    EXPECT_NE(content.find("@HD\tVN:1.6"), std::string::npos);
    EXPECT_NE(content.find("@SQ\tSN:chr14\tLN:107043718"), std::string::npos);
    EXPECT_NE(content.find("@SQ\tSN:chr22\tLN:50818468"), std::string::npos);
    EXPECT_NE(content.find("@RG\tID:1"), std::string::npos);
    EXPECT_NE(content.find("@PG\tID:llmap"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, NoHeaderWhenDisabled) {
    auto tmp = std::filesystem::temp_directory_path() / "test_noheader.sam";
    auto refs = TestRefs();

    BamWriterConfig config;
    config.write_header = false;

    auto writer = BamWriter::Create(tmp, refs, config);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_EQ(content.find("@HD"), std::string::npos);
    EXPECT_EQ(content.find("@SQ"), std::string::npos);

    std::filesystem::remove(tmp);
}

// ========== Mapped record tests ==========

TEST(BamWriterTest, WriteMappedRecord) {
    auto tmp = std::filesystem::temp_directory_path() / "test_mapped.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);

    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 105000000;
    primary.end = 105001000;
    primary.cigar.ops = "1000M";
    primary.score = 1000;
    primary.nm = 5;

    auto rec = make_mapped("read001", 1000, primary);
    rec.cluster_id = 42;
    rec.collapsed_at_level = 2;
    rec.collapsed_at_iteration = 15;

    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto content = ReadFile(tmp);

    // Check alignment line
    EXPECT_NE(content.find("read001"), std::string::npos);
    EXPECT_NE(content.find("chr14"), std::string::npos);
    EXPECT_NE(content.find("105000001"), std::string::npos);  // 1-based
    EXPECT_NE(content.find("1000M"), std::string::npos);
    EXPECT_NE(content.find("AS:i:1000"), std::string::npos);
    EXPECT_NE(content.find("NM:i:5"), std::string::npos);
    EXPECT_NE(content.find("XC:i:42"), std::string::npos);
    EXPECT_NE(content.find("XL:i:2"), std::string::npos);
    EXPECT_NE(content.find("XI:i:15"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, WriteMappedWithAlternatives) {
    auto tmp = std::filesystem::temp_directory_path() / "test_alts.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);

    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;
    primary.cigar.ops = "1000M";
    primary.score = 950;
    primary.nm = 10;

    AlignmentHit alt1;
    alt1.target_id = "chr22";
    alt1.start = 200000;
    alt1.end = 201000;
    alt1.cigar.ops = "1000M";
    alt1.score = 900;
    alt1.nm = 15;

    auto rec = make_mapped("read002", 1000, primary, {alt1});
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_NE(content.find("XA:Z:chr22"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, MappedWithParalogTags) {
    auto tmp = std::filesystem::temp_directory_path() / "test_paralog.sam";
    auto refs = TestRefs();

    BamWriterConfig config;
    config.include_paralog_tags = true;

    auto writer = BamWriter::Create(tmp, refs, config);
    ASSERT_NE(writer, nullptr);

    AlignmentHit primary;
    primary.target_id = "chr14";
    primary.start = 100000;
    primary.end = 101000;
    primary.cigar.ops = "1000M";
    primary.score = 950;
    primary.nm = 10;

    auto rec = make_mapped("read003", 1000, primary);
    rec.paralog_assignment = ParalogCall{};
    rec.paralog_assignment->n_discriminating_psvs = 5;
    rec.paralog_assignment->p_canonical = 0.85f;
    rec.paralog_assignment->inter_paralog = {{"IGHV1-2", 0.85f}, {"IGHV1-3", 0.15f}};

    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_NE(content.find("PD:i:5"), std::string::npos);
    EXPECT_NE(content.find("PC:f:0.85"), std::string::npos);
    EXPECT_NE(content.find("PP:Z:IGHV1-2"), std::string::npos);

    std::filesystem::remove(tmp);
}

// ========== Unmapped record tests ==========

TEST(BamWriterTest, WriteUnmappedRecord) {
    auto tmp = std::filesystem::temp_directory_path() / "test_unmapped.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);

    auto rec = make_unmapped("read004", 500, RejectionReason::NoSeeds);
    rec.cluster_id = 99;

    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_NE(content.find("read004"), std::string::npos);
    EXPECT_NE(content.find("\t4\t"), std::string::npos);  // FLAG=4 (unmapped)
    EXPECT_NE(content.find("XR:Z:NO_SEEDS"), std::string::npos);
    EXPECT_NE(content.find("XC:i:99"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, SkipUnmappedWhenDisabled) {
    auto tmp = std::filesystem::temp_directory_path() / "test_skip_unmapped.sam";
    auto refs = TestRefs();

    BamWriterConfig config;
    config.write_unmapped = false;

    auto writer = BamWriter::Create(tmp, refs, config);
    ASSERT_NE(writer, nullptr);

    auto rec = make_unmapped("read005", 500, RejectionReason::LowComplexity);
    EXPECT_TRUE(writer->Write(rec));  // Returns true but skips
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_EQ(content.find("read005"), std::string::npos);

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.unmapped_written, 0u);

    std::filesystem::remove(tmp);
}

// ========== Tentative record tests ==========

TEST(BamWriterTest, WriteTentativeRecord) {
    auto tmp = std::filesystem::temp_directory_path() / "test_tentative.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);

    std::vector<TentativeTarget> targets = {
        {"chr14", 100000, 101000, 50, 800, 0.85f, 0.45f},
        {"chr22", 200000, 201000, 40, 700, 0.80f, 0.35f}
    };

    auto rec = make_tentative("read006", 1000, targets,
                              RejectionReason::DidNotConverge);

    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_NE(content.find("read006"), std::string::npos);
    EXPECT_NE(content.find("XR:Z:NO_CONVERGE"), std::string::npos);
    EXPECT_NE(content.find("XT:Z:"), std::string::npos);
    EXPECT_NE(content.find("chr14:100000-101000"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(BamWriterTest, SkipTentativeWhenDisabled) {
    auto tmp = std::filesystem::temp_directory_path() / "test_skip_tentative.sam";
    auto refs = TestRefs();

    BamWriterConfig config;
    config.include_tentative = false;

    auto writer = BamWriter::Create(tmp, refs, config);
    ASSERT_NE(writer, nullptr);

    auto rec = make_tentative("read007", 1000, {},
                              RejectionReason::AmbiguousNoAnchor);
    EXPECT_TRUE(writer->Write(rec));
    writer->Close();

    auto content = ReadFile(tmp);
    EXPECT_EQ(content.find("read007"), std::string::npos);

    std::filesystem::remove(tmp);
}

// ========== Batch write tests ==========

TEST(BamWriterTest, WriteBatch) {
    auto tmp = std::filesystem::temp_directory_path() / "test_batch.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);

    std::vector<AlignmentRecord> records;

    // Mapped
    AlignmentHit hit1;
    hit1.target_id = "chr14";
    hit1.start = 100000;
    hit1.end = 101000;
    hit1.cigar.ops = "1000M";
    hit1.score = 950;
    hit1.nm = 5;
    records.push_back(make_mapped("batch001", 1000, hit1));

    // Unmapped
    records.push_back(make_unmapped("batch002", 500, RejectionReason::NoSeeds));

    // Tentative
    records.push_back(make_tentative("batch003", 800, {},
                                     RejectionReason::DidNotConverge));

    EXPECT_TRUE(writer->WriteBatch(records));
    writer->Close();

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.records_written, 3u);
    EXPECT_EQ(stats.mapped_written, 1u);
    EXPECT_EQ(stats.unmapped_written, 1u);
    EXPECT_EQ(stats.tentative_written, 1u);

    std::filesystem::remove(tmp);
}

// ========== Statistics tests ==========

TEST(BamWriterTest, TrackStatistics) {
    auto tmp = std::filesystem::temp_directory_path() / "test_stats.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    hit.cigar.ops = "1000M";
    hit.score = 950;
    hit.nm = 5;

    AlignmentHit alt;
    alt.target_id = "chr22";
    alt.start = 200000;
    alt.end = 201000;
    alt.cigar.ops = "1000M";
    alt.score = 900;
    alt.nm = 10;

    auto rec = make_mapped("stats001", 1000, hit, {alt});
    writer->Write(rec);

    auto rec2 = make_unmapped("stats002", 500, RejectionReason::NoSeeds);
    writer->Write(rec2);

    writer->Close();

    auto stats = writer->GetStats();
    EXPECT_EQ(stats.records_written, 2u);
    EXPECT_EQ(stats.mapped_written, 1u);
    EXPECT_EQ(stats.unmapped_written, 1u);
    EXPECT_EQ(stats.alternatives_written, 1u);
    EXPECT_GT(stats.write_time_ms, 0.0f);

    std::filesystem::remove(tmp);
}

// ========== CIGAR utility tests ==========

TEST(CigarTest, GenerateSimpleCigarExactMatch) {
    auto cig = cigar::GenerateSimpleCigar(1000, 0, 1000, 0);
    EXPECT_EQ(cig, "1000M");
}

TEST(CigarTest, GenerateSimpleCigarWithMismatches) {
    auto cig = cigar::GenerateSimpleCigar(1000, 0, 1000, 5);
    EXPECT_EQ(cig, "1000M");  // Same length, all differences are substitutions
}

TEST(CigarTest, GenerateSimpleCigarWithDeletion) {
    auto cig = cigar::GenerateSimpleCigar(1000, 0, 1010, 10);
    // Should contain a deletion
    EXPECT_NE(cig.find("D"), std::string::npos);
}

TEST(CigarTest, GenerateSimpleCigarWithInsertion) {
    auto cig = cigar::GenerateSimpleCigar(1010, 0, 1000, 10);
    // Should contain an insertion
    EXPECT_NE(cig.find("I"), std::string::npos);
}

TEST(CigarTest, CigarStatsMatch) {
    auto [qc, rc] = cigar::CigarStats("100M");
    EXPECT_EQ(qc, 100u);
    EXPECT_EQ(rc, 100u);
}

TEST(CigarTest, CigarStatsWithIndels) {
    auto [qc, rc] = cigar::CigarStats("50M10I50M10D50M");
    EXPECT_EQ(qc, 160u);  // 50+10+50+50
    EXPECT_EQ(rc, 160u);  // 50+50+10+50
}

TEST(CigarTest, CigarStatsWithSoftClip) {
    auto [qc, rc] = cigar::CigarStats("10S90M10S");
    EXPECT_EQ(qc, 110u);  // 10+90+10
    EXPECT_EQ(rc, 90u);   // only M consumes reference
}

TEST(CigarTest, ValidateCigarValid) {
    EXPECT_TRUE(cigar::ValidateCigar("100M"));
    EXPECT_TRUE(cigar::ValidateCigar("50M10I40M"));
    EXPECT_TRUE(cigar::ValidateCigar("10S90M"));
    EXPECT_TRUE(cigar::ValidateCigar("*"));
    EXPECT_TRUE(cigar::ValidateCigar(""));
}

TEST(CigarTest, ValidateCigarInvalid) {
    EXPECT_FALSE(cigar::ValidateCigar("M100"));  // Op before number
    EXPECT_FALSE(cigar::ValidateCigar("100"));   // No op
    EXPECT_FALSE(cigar::ValidateCigar("100Z"));  // Invalid op
}

// ========== Convenience function tests ==========

TEST(BamWriterTest, WriteAlignmentsConvenience) {
    auto tmp = std::filesystem::temp_directory_path() / "test_convenience.sam";
    auto refs = TestRefs();

    std::vector<AlignmentRecord> records;

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    hit.cigar.ops = "1000M";
    hit.score = 950;
    hit.nm = 5;
    records.push_back(make_mapped("conv001", 1000, hit));
    records.push_back(make_unmapped("conv002", 500, RejectionReason::NoSeeds));

    EXPECT_TRUE(WriteAlignments(tmp, records, refs));

    auto content = ReadFile(tmp);
    EXPECT_NE(content.find("conv001"), std::string::npos);
    EXPECT_NE(content.find("conv002"), std::string::npos);

    std::filesystem::remove(tmp);
}

// ========== Error handling tests ==========

TEST(BamWriterTest, FailsOnInvalidPath) {
    auto refs = TestRefs();
    auto writer = BamWriter::Create("/nonexistent/path/test.sam", refs);
    EXPECT_EQ(writer, nullptr);
}

TEST(BamWriterTest, ReturnsErrorAfterClose) {
    auto tmp = std::filesystem::temp_directory_path() / "test_closed.sam";
    auto refs = TestRefs();

    auto writer = BamWriter::Create(tmp, refs);
    ASSERT_NE(writer, nullptr);
    writer->Close();

    AlignmentHit hit;
    hit.target_id = "chr14";
    hit.start = 100000;
    hit.end = 101000;
    hit.cigar.ops = "1000M";
    auto rec = make_mapped("fail001", 1000, hit);

    EXPECT_FALSE(writer->Write(rec));
    EXPECT_FALSE(writer->LastError().empty());

    std::filesystem::remove(tmp);
}

}  // namespace
}  // namespace llmap::output
