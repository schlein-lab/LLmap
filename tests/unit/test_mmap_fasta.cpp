// Unit tests for MmapFastaReader

#include <gtest/gtest.h>

#include "io/mmap_fasta.h"

#include <filesystem>
#include <fstream>
#include <string>

namespace llmap::io {
namespace {

class MmapFastaTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_mmap_fasta_test";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path CreateTestFasta(const std::string& name,
                                          const std::string& content) {
        auto path = test_dir_ / name;
        std::ofstream out(path);
        out << content;
        return path;
    }

    std::filesystem::path test_dir_;
};

// ========== Basic Functionality ==========

TEST_F(MmapFastaTest, OpenNonexistentFile) {
    MmapFastaReader reader("/nonexistent/file.fasta");
    EXPECT_FALSE(reader.IsValid());
    EXPECT_FALSE(reader.LastError().empty());
}

TEST_F(MmapFastaTest, OpenEmptyFile) {
    auto path = CreateTestFasta("empty.fa", "");
    MmapFastaReader reader(path);
    EXPECT_FALSE(reader.IsValid());
}

TEST_F(MmapFastaTest, ReadSingleSequence) {
    auto path = CreateTestFasta("single.fa",
        ">chr1\n"
        "ACGTACGT\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());
    EXPECT_EQ(reader.NumSequences(), 1);

    auto seq = reader.GetSequence(0);
    ASSERT_TRUE(seq.has_value());
    EXPECT_EQ(seq->name, "chr1");
    EXPECT_EQ(seq->length, 8);

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, "ACGTACGT");
}

TEST_F(MmapFastaTest, ReadMultipleSequences) {
    auto path = CreateTestFasta("multi.fa",
        ">chr1\n"
        "ACGT\n"
        ">chr2\n"
        "TGCATGCA\n"
        ">chr3\n"
        "GGCC\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());
    EXPECT_EQ(reader.NumSequences(), 3);

    auto seq1 = reader.GetSequence(0);
    ASSERT_TRUE(seq1.has_value());
    EXPECT_EQ(seq1->name, "chr1");

    auto seq2 = reader.GetSequence(1);
    ASSERT_TRUE(seq2.has_value());
    EXPECT_EQ(seq2->name, "chr2");
    EXPECT_EQ(seq2->length, 8);

    auto seq3 = reader.GetSequence(2);
    ASSERT_TRUE(seq3.has_value());
    EXPECT_EQ(seq3->name, "chr3");
}

TEST_F(MmapFastaTest, GetSequenceByName) {
    auto path = CreateTestFasta("named.fa",
        ">chr1\n"
        "AAAA\n"
        ">chrX\n"
        "CCCC\n"
        ">chrM\n"
        "GGGG\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto seq = reader.GetSequence("chrX");
    ASSERT_TRUE(seq.has_value());
    EXPECT_EQ(seq->name, "chrX");
    EXPECT_EQ(seq->length, 4);

    std::string data = reader.GetSequenceData("chrX");
    EXPECT_EQ(data, "CCCC");

    auto missing = reader.GetSequence("chrY");
    EXPECT_FALSE(missing.has_value());
}

TEST_F(MmapFastaTest, SequenceNames) {
    auto path = CreateTestFasta("names.fa",
        ">seq_a\nACGT\n"
        ">seq_b\nTGCA\n"
        ">seq_c\nGGCC\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto names = reader.SequenceNames();
    ASSERT_EQ(names.size(), 3);
    EXPECT_EQ(names[0], "seq_a");
    EXPECT_EQ(names[1], "seq_b");
    EXPECT_EQ(names[2], "seq_c");
}

// ========== Multi-line Sequences ==========

TEST_F(MmapFastaTest, MultiLineSequence) {
    auto path = CreateTestFasta("multiline.fa",
        ">chr1\n"
        "ACGTACGT\n"
        "TGCATGCA\n"
        "GGCCGGCC\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());
    EXPECT_EQ(reader.NumSequences(), 1);

    auto seq = reader.GetSequence(0);
    ASSERT_TRUE(seq.has_value());
    EXPECT_EQ(seq->length, 24);  // 8 + 8 + 8

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, "ACGTACGTTGCATGCAGGCCGGCC");
}

TEST_F(MmapFastaTest, WrappedReference) {
    // Typical reference format with 80-char lines
    std::string seq(240, 'A');  // 240 bases
    std::string wrapped;
    for (size_t i = 0; i < seq.size(); i += 80) {
        wrapped += seq.substr(i, 80) + "\n";
    }

    auto path = CreateTestFasta("wrapped.fa", ">chr1\n" + wrapped);

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto info = reader.GetSequence(0);
    ASSERT_TRUE(info.has_value());
    EXPECT_EQ(info->length, 240);

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, seq);
}

// ========== Case Handling ==========

TEST_F(MmapFastaTest, LowercaseToUppercase) {
    auto path = CreateTestFasta("lower.fa",
        ">chr1\n"
        "acgtACGT\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, "ACGTACGT");  // All uppercase
}

TEST_F(MmapFastaTest, NBasesPreserved) {
    auto path = CreateTestFasta("nbases.fa",
        ">chr1\n"
        "ACGTnnnTGCA\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, "ACGTNNNTGCA");  // N bases preserved but uppercased
}

// ========== Subsequence Extraction ==========

TEST_F(MmapFastaTest, GetSubsequence) {
    auto path = CreateTestFasta("subseq.fa",
        ">chr1\n"
        "ABCDEFGHIJKLMNOP\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    // Get middle portion
    std::string sub = reader.GetSubsequence(0, 4, 4);
    EXPECT_EQ(sub, "EFGH");

    // Get from start
    sub = reader.GetSubsequence(0, 0, 3);
    EXPECT_EQ(sub, "ABC");

    // Get from end
    sub = reader.GetSubsequence(0, 14, 2);
    EXPECT_EQ(sub, "OP");

    // Get beyond end (returns what's available)
    sub = reader.GetSubsequence(0, 14, 10);
    EXPECT_EQ(sub, "OP");
}

TEST_F(MmapFastaTest, GetSubsequenceByName) {
    auto path = CreateTestFasta("subseq_name.fa",
        ">chr1\nAAAA\n"
        ">chr2\nBCDEFGHI\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    std::string sub = reader.GetSubsequence("chr2", 2, 3);
    EXPECT_EQ(sub, "DEF");
}

TEST_F(MmapFastaTest, SubsequenceAcrossLines) {
    auto path = CreateTestFasta("subseq_multiline.fa",
        ">chr1\n"
        "AAAA\n"
        "BBBB\n"
        "CCCC\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    // Crosses line boundary
    std::string sub = reader.GetSubsequence(0, 2, 6);
    EXPECT_EQ(sub, "AABBBB");
}

// ========== Raw Access ==========

TEST_F(MmapFastaTest, GetSequenceRaw) {
    auto path = CreateTestFasta("raw.fa",
        ">chr1\n"
        "ACGT\n"
        "TGCA\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto raw = reader.GetSequenceRaw(0);
    // Raw includes newlines
    EXPECT_TRUE(raw.find('\n') != std::string_view::npos);
    EXPECT_EQ(raw.size(), 10);  // "ACGT\nTGCA\n"
}

// ========== Statistics ==========

TEST_F(MmapFastaTest, GetStats) {
    auto path = CreateTestFasta("stats.fa",
        ">chr1\nACGT\n"
        ">chr2\nTGCATGCA\n"
        ">chr3\nGG\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto stats = reader.GetStats();
    EXPECT_GT(stats.file_size, 0);
    EXPECT_EQ(stats.num_sequences, 3);
    EXPECT_EQ(stats.total_bases, 14);  // 4 + 8 + 2
}

// ========== Move Semantics ==========

TEST_F(MmapFastaTest, MoveConstruct) {
    auto path = CreateTestFasta("move.fa", ">chr1\nACGT\n");

    MmapFastaReader reader1(path);
    ASSERT_TRUE(reader1.IsValid());

    MmapFastaReader reader2(std::move(reader1));
    EXPECT_TRUE(reader2.IsValid());
    EXPECT_EQ(reader2.NumSequences(), 1);
}

TEST_F(MmapFastaTest, MoveAssign) {
    auto path = CreateTestFasta("move.fa", ">chr1\nACGT\n");

    MmapFastaReader reader1(path);
    MmapFastaReader reader2("/nonexistent");

    reader2 = std::move(reader1);
    EXPECT_TRUE(reader2.IsValid());
    EXPECT_EQ(reader2.NumSequences(), 1);
}

// ========== Edge Cases ==========

TEST_F(MmapFastaTest, SequenceWithDescription) {
    auto path = CreateTestFasta("desc.fa",
        ">chr1 some description here\n"
        "ACGT\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto seq = reader.GetSequence(0);
    ASSERT_TRUE(seq.has_value());
    EXPECT_EQ(seq->name, "chr1");  // Name is just the ID, not the description
}

TEST_F(MmapFastaTest, EmptyLinesInFile) {
    auto path = CreateTestFasta("empty_lines.fa",
        "\n"
        ">chr1\n"
        "\n"
        "ACGT\n"
        "\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());
    EXPECT_EQ(reader.NumSequences(), 1);

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, "ACGT");
}

TEST_F(MmapFastaTest, WindowsLineEndings) {
    auto path = CreateTestFasta("windows.fa",
        ">chr1\r\n"
        "ACGT\r\n"
    );

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());
    EXPECT_EQ(reader.NumSequences(), 1);

    std::string data = reader.GetSequenceData(0);
    EXPECT_EQ(data, "ACGT");
}

TEST_F(MmapFastaTest, OutOfBoundsIndex) {
    auto path = CreateTestFasta("bounds.fa", ">chr1\nACGT\n");

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto seq = reader.GetSequence(999);
    EXPECT_FALSE(seq.has_value());

    auto raw = reader.GetSequenceRaw(999);
    EXPECT_TRUE(raw.empty());

    auto data = reader.GetSequenceData(999);
    EXPECT_TRUE(data.empty());
}

TEST_F(MmapFastaTest, InvalidSequenceName) {
    auto path = CreateTestFasta("invalid_name.fa", ">chr1\nACGT\n");

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    auto seq = reader.GetSequence("nonexistent");
    EXPECT_FALSE(seq.has_value());
}

// ========== Configuration Options ==========

TEST_F(MmapFastaTest, BuildIndexOnOpen) {
    auto path = CreateTestFasta("index.fa",
        ">chr1\nACGT\n"
        ">chr2\nTGCA\n"
    );

    MmapFastaConfig config;
    config.build_index_on_open = true;

    MmapFastaReader reader(path, config);
    ASSERT_TRUE(reader.IsValid());
    EXPECT_EQ(reader.NumSequences(), 2);
}

// ========== Memory Advice ==========

TEST_F(MmapFastaTest, AdviseSequential) {
    auto path = CreateTestFasta("advise.fa", ">chr1\nACGT\n");

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    // Should not crash
    reader.AdviseSequential();
}

TEST_F(MmapFastaTest, AdviseRandom) {
    auto path = CreateTestFasta("advise.fa", ">chr1\nACGT\n");

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    // Should not crash
    reader.AdviseRandom();
}

TEST_F(MmapFastaTest, AdviseWillNeed) {
    auto path = CreateTestFasta("advise.fa", ">chr1\nACGT\n>chr2\nTGCA\n");

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    // Should not crash
    reader.AdviseWillNeed(0);
    reader.AdviseWillNeed(1);
    reader.AdviseWillNeed(999);  // Out of bounds - should be no-op
}

TEST_F(MmapFastaTest, AdviseDontNeed) {
    auto path = CreateTestFasta("advise.fa", ">chr1\nACGT\n>chr2\nTGCA\n");

    MmapFastaReader reader(path);
    ASSERT_TRUE(reader.IsValid());

    // Should not crash
    reader.AdviseDontNeed(0);
    reader.AdviseDontNeed(1);
}

// ========== Utility Functions ==========

TEST_F(MmapFastaTest, IsFastaFileByExtension) {
    auto fa = CreateTestFasta("test.fa", ">chr1\nACGT\n");
    auto fasta = CreateTestFasta("test.fasta", ">chr1\nACGT\n");
    auto fna = CreateTestFasta("test.fna", ">chr1\nACGT\n");
    auto txt = CreateTestFasta("test.txt", ">chr1\nACGT\n");

    EXPECT_TRUE(IsFastaFile(fa));
    EXPECT_TRUE(IsFastaFile(fasta));
    EXPECT_TRUE(IsFastaFile(fna));
    EXPECT_FALSE(IsFastaFile(txt));
}

TEST_F(MmapFastaTest, IsFastaFileByContent) {
    auto valid = CreateTestFasta("valid.fa", ">chr1\nACGT\n");
    auto invalid = CreateTestFasta("invalid.fa", "not a fasta file\n");

    EXPECT_TRUE(IsFastaFile(valid));
    EXPECT_FALSE(IsFastaFile(invalid));
}

// ========== Path Access ==========

TEST_F(MmapFastaTest, PathAccessor) {
    auto path = CreateTestFasta("path.fa", ">chr1\nACGT\n");

    MmapFastaReader reader(path);
    EXPECT_EQ(reader.Path(), path);
}

}  // namespace
}  // namespace llmap::io
