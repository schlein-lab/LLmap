// LLmap — JunctionDb tests.
//
// Synthesise tiny inputs for each loader and verify the merge semantics:
// GENCODE + GTEx + ChessDB + circBase all contribute boolean flags on a
// single JunctionEvidence; HasAnyEvidence covers the OR; Lookup returns
// the merged struct.

#include "annot/junction_db.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

using namespace llmap::annot;

namespace {

std::filesystem::path WriteTmp(const std::string& name,
                                const std::string& body) {
    auto p = std::filesystem::temp_directory_path() / name;
    std::ofstream out(p);
    out << body;
    return p;
}

}  // namespace

TEST(JunctionDb, EmptyDbReturnsAllFlagsFalse) {
    JunctionDb d;
    auto ev = d.Lookup("chr14", 1000, 5000);
    EXPECT_FALSE(ev.in_gencode);
    EXPECT_FALSE(ev.in_gtex);
    EXPECT_FALSE(d.HasAnyEvidence("chr14", 1000, 5000));
    EXPECT_EQ(d.TotalJunctions(), 0u);
}

TEST(JunctionDb, LoadGencodeMinimalGff) {
    auto p = WriteTmp("llmap_junc_gencode.gff3",
        "##gff-version 3\n"
        "chr14\tGENCODE\texon\t100\t300\t.\t+\t.\ttranscript_id=ENST.1\n"
        "chr14\tGENCODE\texon\t500\t1000\t.\t+\t.\ttranscript_id=ENST.1\n");
    JunctionDb d;
    ASSERT_TRUE(d.LoadGencode(p));
    // expected single junction donor=300-1=299, acceptor=500-1=499
    auto ev = d.Lookup("chr14", 299, 499);
    EXPECT_TRUE(ev.in_gencode);
    EXPECT_EQ(d.TotalJunctions(), 1u);
    std::filesystem::remove(p);
}

TEST(JunctionDb, LoadGtexBedMergesIntoExistingRecord) {
    JunctionDb d;
    // pre-populate via GENCODE
    auto gff = WriteTmp("llmap_junc_gencode2.gff3",
        "chr14\tGENCODE\texon\t100\t300\t.\t+\t.\ttranscript_id=ENST.1\n"
        "chr14\tGENCODE\texon\t500\t1000\t.\t+\t.\ttranscript_id=ENST.1\n");
    d.LoadGencode(gff);

    // GTEx BED: chrom, donor, acceptor, name, sample_count, strand
    auto bed = WriteTmp("llmap_junc_gtex.bed",
        "chr14\t299\t499\tj1\t127\t+\n");
    d.LoadGtexJunctions(bed);

    auto ev = d.Lookup("chr14", 299, 499);
    EXPECT_TRUE(ev.in_gencode);
    EXPECT_TRUE(ev.in_gtex);
    EXPECT_EQ(ev.gtex_sample_count, 127u);
    EXPECT_TRUE(d.HasAnyEvidence("chr14", 299, 499));

    std::filesystem::remove(gff);
    std::filesystem::remove(bed);
}

TEST(JunctionDb, LoadCircRnaInsertsBackSpliceFormat) {
    JunctionDb d;
    // circRNA BED: chrom, back-splice acceptor, back-splice donor, name, ., strand
    auto bed = WriteTmp("llmap_junc_circ.bed",
        "chr14\t1000\t5000\tcirc1\t.\t+\n");
    d.LoadCircRnaDb(bed);

    auto ev = d.Lookup("chr14", /*donor=*/5000, /*acceptor=*/1000);
    EXPECT_TRUE(ev.in_circ_db);
    std::filesystem::remove(bed);
}

TEST(JunctionDb, JunctionsOnChromCounter) {
    JunctionDb d;
    auto gff = WriteTmp("llmap_junc_count.gff3",
        "chr14\tGENCODE\texon\t100\t300\t.\t+\t.\ttranscript_id=ENST.1\n"
        "chr14\tGENCODE\texon\t500\t1000\t.\t+\t.\ttranscript_id=ENST.1\n"
        "chr15\tGENCODE\texon\t10\t50\t.\t+\t.\ttranscript_id=ENST.2\n"
        "chr15\tGENCODE\texon\t200\t300\t.\t+\t.\ttranscript_id=ENST.2\n");
    d.LoadGencode(gff);
    EXPECT_EQ(d.JunctionsOnChrom("chr14"), 1u);
    EXPECT_EQ(d.JunctionsOnChrom("chr15"), 1u);
    EXPECT_EQ(d.JunctionsOnChrom("chr_missing"), 0u);
    EXPECT_EQ(d.TotalJunctions(), 2u);
    std::filesystem::remove(gff);
}

TEST(JunctionDb, MissingFileLoadFailsGracefully) {
    JunctionDb d;
    EXPECT_FALSE(d.LoadGencode("/tmp/does_not_exist_zzzz_42.gff3"));
    EXPECT_FALSE(d.LoadGtexJunctions("/tmp/does_not_exist_zzzz_42.bed"));
    EXPECT_FALSE(d.LoadCircRnaDb("/tmp/does_not_exist_zzzz_42.bed"));
    EXPECT_FALSE(d.LoadChessDbJunctions("/tmp/does_not_exist_zzzz_42.gtf"));
    EXPECT_EQ(d.TotalJunctions(), 0u);
}
