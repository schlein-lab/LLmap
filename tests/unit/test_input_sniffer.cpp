// LLmap — input_sniffer unit tests.
//
// Covers SniffFormat (magic bytes / first line), ComputeFastaStats
// (median / N50 / n_seqs over a sample), and ResolveMode (the full §3
// Transcript-Mode resolution heuristic + override + reads_vs_assembly).
//
// All inputs are tiny synthetic temp files — no genomic data, no NAS.

#include "io/input_sniffer.h"
#include "core/transcript_mode.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

using llmap::core::TranscriptMode;
using llmap::io::ComputeFastaStats;
using llmap::io::FileFormat;
using llmap::io::ResolveMode;
using llmap::io::SniffFormat;

// Write `content` to a uniquely-named temp file with the given basename,
// return its path. The basename matters for the FLNC/isoseq heuristics.
std::string WriteTemp(const std::string& basename, const std::string& content) {
    namespace fs = std::filesystem;
    static int counter = 0;
    fs::path dir = fs::path(::testing::TempDir()) /
                   ("sniff_" + std::to_string(counter++));
    fs::create_directories(dir);
    fs::path p = dir / basename;
    std::ofstream out(p, std::ios::binary);
    out << content;
    out.close();
    return p.string();
}

// Build a FASTA string with n records each of length `len`.
std::string MakeFasta(std::size_t n, std::size_t len) {
    std::string s;
    for (std::size_t i = 0; i < n; ++i) {
        s += ">seq" + std::to_string(i) + "\n";
        s += std::string(len, 'A');
        s += "\n";
    }
    return s;
}

// ---------------------------------------------------------------------------
// SniffFormat
// ---------------------------------------------------------------------------

TEST(SniffFormat, Fasta) {
    const std::string p = WriteTemp("ref.fa", ">chr1\nACGTACGT\n");
    EXPECT_EQ(SniffFormat(p), FileFormat::Fasta);
}

TEST(SniffFormat, FastqVsSamDisambiguation) {
    // FASTQ: '@' + read name (not two uppercase + tab).
    const std::string fq =
        WriteTemp("reads.fq", "@read1\nACGTACGT\n+\nIIIIIIII\n");
    EXPECT_EQ(SniffFormat(fq), FileFormat::Fastq);

    // SAM: '@' + two uppercase header tag + tab.
    const std::string sam =
        WriteTemp("aln.sam", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n");
    EXPECT_EQ(SniffFormat(sam), FileFormat::Sam);
}

TEST(SniffFormat, EmptyOrMissing) {
    EXPECT_EQ(SniffFormat("/nonexistent/path/xyz"), FileFormat::Unknown);
    const std::string empty = WriteTemp("empty.dat", "");
    EXPECT_EQ(SniffFormat(empty), FileFormat::Unknown);
}

// ---------------------------------------------------------------------------
// ComputeFastaStats — median / N50
// ---------------------------------------------------------------------------

TEST(FastaStats, MedianAndN50KnownLengths) {
    // Lengths: 100, 200, 300, 400, 500.
    std::string s;
    const std::size_t lens[] = {100, 200, 300, 400, 500};
    for (std::size_t i = 0; i < 5; ++i) {
        s += ">s" + std::to_string(i) + "\n" + std::string(lens[i], 'C') + "\n";
    }
    const std::string p = WriteTemp("multi.fa", s);
    auto st = ComputeFastaStats(p);
    EXPECT_EQ(st.n_seqs, 5u);
    EXPECT_EQ(st.median_len, 300u);  // middle of sorted [100..500]
    // N50: total=1500, half=750. Sorted desc 500,400,300 -> cumsum 500,900>=750
    // => N50 = 400.
    EXPECT_EQ(st.n50, 400u);
}

TEST(FastaStats, MultiLineSequence) {
    const std::string p = WriteTemp("wrap.fa", ">a\nACGT\nACGT\nAC\n");  // len 10
    auto st = ComputeFastaStats(p);
    EXPECT_EQ(st.n_seqs, 1u);
    EXPECT_EQ(st.median_len, 10u);
}

// ---------------------------------------------------------------------------
// ResolveMode — the §3 heuristic
// ---------------------------------------------------------------------------

TEST(ResolveMode, OverrideWins) {
    const std::string p = WriteTemp("ref.fa", MakeFasta(3, 200000));
    auto r = ResolveMode(p, TranscriptMode::Transcript,
                         /*has_reads=*/true, /*has_assembly=*/false);
    EXPECT_EQ(r.mode, TranscriptMode::Transcript);  // override beats sniff
}

TEST(ResolveMode, ReadsPlusAssembly) {
    const std::string p = WriteTemp("reads.fq", "@r\nACGT\n+\nIIII\n");
    auto r = ResolveMode(p, TranscriptMode::Auto,
                         /*has_reads=*/true, /*has_assembly=*/true);
    EXPECT_EQ(r.mode, TranscriptMode::ReadsVsAssembly);
}

TEST(ResolveMode, FastaAssembly) {
    // 3 contigs of 200 kb -> median>50k, N50>100k, n_seqs<5k.
    const std::string p = WriteTemp("hifiasm.hap1.fa", MakeFasta(3, 200000));
    auto r = ResolveMode(p, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::Assembly);
    ASSERT_TRUE(r.fasta_stats.has_value());
    EXPECT_EQ(r.fasta_stats->n_seqs, 3u);
}

TEST(ResolveMode, FastaShortReadDefaultsGenome) {
    // Short FASTA records -> not assembly, not FLNC -> GenomeReads.
    const std::string p = WriteTemp("contigs.fa", MakeFasta(10, 150));
    auto r = ResolveMode(p, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::GenomeReads);
}

TEST(ResolveMode, FastaFlncByBasename) {
    // Medium-length FASTA but few records -> length rule misses; basename
    // hint ("flnc") promotes it to Transcript.
    const std::string p = WriteTemp("sample.flnc.fa", MakeFasta(5, 1500));
    auto r = ResolveMode(p, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::Transcript);
}

TEST(ResolveMode, FastqFlncBasename) {
    const std::string p = WriteTemp("lib.isoseq.fastq", "@r\nACGT\n+\nIIII\n");
    auto r = ResolveMode(p, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::Transcript);
}

TEST(ResolveMode, FastqDnaBasename) {
    const std::string p = WriteTemp("wgs_reads.fastq", "@r\nACGT\n+\nIIII\n");
    auto r = ResolveMode(p, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::GenomeReads);
}

TEST(ResolveMode, SamSpliceAwarePg) {
    const std::string sam = WriteTemp(
        "aln.sam",
        "@HD\tVN:1.6\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        "@PG\tID:minimap2\tPN:minimap2\tCL:minimap2 -ax splice ref.fa r.fq\n");
    auto r = ResolveMode(sam, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::Transcript);
}

TEST(ResolveMode, SamDnaPgDefaultsGenome) {
    const std::string sam = WriteTemp(
        "aln.sam",
        "@HD\tVN:1.6\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        "@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref.fa r.fq\n");
    auto r = ResolveMode(sam, TranscriptMode::Auto, true, false);
    EXPECT_EQ(r.mode, TranscriptMode::GenomeReads);
}

}  // namespace
