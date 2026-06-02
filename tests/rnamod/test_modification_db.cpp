// LLmap — ModificationDb tests.

#include "rnamod/modification_db.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

using namespace llmap::rnamod;

namespace {

std::filesystem::path WriteTmp(const std::string& name,
                                const std::string& body) {
    auto p = std::filesystem::temp_directory_path() / name;
    std::ofstream out(p);
    out << body;
    return p;
}

}  // namespace

TEST(ModificationDb, NameRoundTripCoversAllKnown) {
    using K = ModificationKind;
    constexpr K all[] = {
        K::Unknown, K::M6A, K::M6Am, K::M5C, K::Hm5C, K::M1A,
        K::Pseudouridine, K::Inosine, K::M7G, K::M3C, K::M2_O_Methyl,
        K::Ac4C, K::M6_2A, K::M1G, K::M22G, K::M5U, K::T6A, K::Mcm5U,
        K::Queosine, K::A2I_Editing, K::C2U_Editing,
        K::Cap0, K::Cap1, K::Cap2,
        K::Polya_Standard, K::Polya_Short, K::Polya_Modified,
        K::Oligo_U_Tail, K::NovelUnclassifiedMod,
    };
    for (auto k : all) {
        const char* n = ModificationKindName(k);
        ASSERT_NE(n, nullptr);
        auto parsed = ParseModificationKind(n);
        ASSERT_TRUE(parsed.has_value()) << "name='" << n << "'";
        EXPECT_EQ(*parsed, k);
    }
    EXPECT_FALSE(ParseModificationKind("xxx").has_value());
}

TEST(ModificationDb, M6APriorRecognisesDrachContext) {
    ModificationDb d;
    d.LoadDefaults();
    EXPECT_NEAR(d.PriorProbability(ModificationKind::M6A, "GGACT"), 0.30f, 0.001f)
        << "GGACT is DRACH (D=G, R=G, A, C, H=T)";
    EXPECT_NEAR(d.PriorProbability(ModificationKind::M6A, "AGACA"), 0.30f, 0.001f);
    EXPECT_NEAR(d.PriorProbability(ModificationKind::M6A, "CCCCC"), 0.005f, 1e-4f)
        << "CCCCC is non-DRACH (R must be A/G)";
}

TEST(ModificationDb, C2UPriorRecognisesAidWrcHotspot) {
    ModificationDb d;
    d.LoadDefaults();
    // Context: [pos-2][pos-1][C][pos+1][pos+2]
    // For hotspot we look at middle 4: R G Y W → e.g. AGCT
    // ctx[1]=A (R), ctx[2]=G, ctx[3]=C (Y), ctx[4]=A (W)
    EXPECT_NEAR(d.PriorProbability(ModificationKind::C2U_Editing, "AAGCA"),
                 0.10f, 0.001f);
    // Non-hotspot
    EXPECT_NEAR(d.PriorProbability(ModificationKind::C2U_Editing, "CCCCC"),
                 0.005f, 1e-4f);
}

TEST(ModificationDb, LoadREPICAttachesSiteToTranscript) {
    auto bed = WriteTmp("llmap_modtest_repic.bed",
        "chr14\t12345\t12346\tENST.1\t127\t+\n"
        "chr14\t99999\t100000\tENST.2\t8\t-\n");
    ModificationDb d;
    EXPECT_TRUE(d.LoadREPIC(bed));
    auto s = d.KnownAt("ENST.1", 12345);
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0], ModificationKind::M6A);
    EXPECT_TRUE(d.KnownAt("ENST.1", 0).empty());
    EXPECT_EQ(d.TotalSites(), 2u);
    std::filesystem::remove(bed);
}

TEST(ModificationDb, LoadREDIportalDistinguishesA2IFromC2U) {
    auto tsv = WriteTmp("llmap_modtest_redi.tsv",
        "chrom\tpos\tref\talt\ttranscript_id\n"
        "chr1\t100\tA\tI\tENST.A\n"
        "chr1\t200\tC\tU\tENST.B\n");
    ModificationDb d;
    EXPECT_TRUE(d.LoadREDIportal(tsv));
    auto s_a = d.KnownAt("ENST.A", 100);
    ASSERT_EQ(s_a.size(), 1u);
    EXPECT_EQ(s_a[0], ModificationKind::A2I_Editing);
    auto s_b = d.KnownAt("ENST.B", 200);
    ASSERT_EQ(s_b.size(), 1u);
    EXPECT_EQ(s_b[0], ModificationKind::C2U_Editing);
    std::filesystem::remove(tsv);
}

TEST(ModificationDb, MissingFileLoadFails) {
    ModificationDb d;
    EXPECT_FALSE(d.LoadREPIC("/tmp/does_not_exist_xyz_42.bed"));
    EXPECT_FALSE(d.LoadREDIportal("/tmp/does_not_exist_xyz_42.tsv"));
    EXPECT_FALSE(d.LoadModomicsExtended("/tmp/does_not_exist_xyz_42.tsv"));
    EXPECT_EQ(d.TotalSites(), 0u);
}
