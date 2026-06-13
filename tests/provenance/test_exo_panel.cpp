// LLmap — exo_panel tests.

#include "provenance/exo_panel.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

namespace {

using llmap::provenance::ExoCategory;
using llmap::provenance::ExoCategoryName;
using llmap::provenance::ExoPanel;

TEST(ExoPanel, BuiltinStarterLookup) {
    ExoPanel p;
    p.LoadBuiltinStarter();
    EXPECT_GE(p.Size(), 5u);

    const auto* ebv = p.Lookup("NC_007605");
    ASSERT_NE(ebv, nullptr);
    EXPECT_EQ(ebv->taxon, "EBV");
    EXPECT_EQ(ebv->pv_tag, "exo:ebv");
    EXPECT_EQ(ebv->category, ExoCategory::Virus);
    EXPECT_FALSE(ebv->verified);

    const auto* phix = p.Lookup("NC_001422");
    ASSERT_NE(phix, nullptr);
    EXPECT_EQ(phix->pv_tag, "spikein:phix");
    EXPECT_EQ(phix->category, ExoCategory::SpikeIn);

    EXPECT_EQ(p.Lookup("chr1"), nullptr);  // host ref → not exogenous
}

TEST(ExoPanel, LoadTsvAuthoritative) {
    namespace fs = std::filesystem;
    const fs::path f = fs::path(::testing::TempDir()) / "exo_panel.tsv";
    {
        std::ofstream out(f);
        out << "# taxon ref_id pv_tag category\n";
        out << "Bradyrhizobium\tNZ_KITOME1\texo:kitome\tbacteria\n";
        out << "Malassezia\tNW_FUNGUS1\texo:fungus\tfungus\n";
    }
    ExoPanel p;
    ASSERT_TRUE(p.LoadTsv(f.string()));
    EXPECT_EQ(p.Size(), 2u);

    const auto* k = p.Lookup("NZ_KITOME1");
    ASSERT_NE(k, nullptr);
    EXPECT_EQ(k->pv_tag, "exo:kitome");
    EXPECT_EQ(k->category, ExoCategory::Bacteria);
    EXPECT_TRUE(k->verified);

    const auto* fung = p.Lookup("NW_FUNGUS1");
    ASSERT_NE(fung, nullptr);
    EXPECT_EQ(fung->category, ExoCategory::Fungus);
    fs::remove(f);
}

TEST(ExoPanel, LoadTsvMissingFileFails) {
    ExoPanel p;
    EXPECT_FALSE(p.LoadTsv("/nonexistent/exo.tsv"));
}

TEST(ExoPanel, CategoryNames) {
    EXPECT_STREQ(ExoCategoryName(ExoCategory::Virus), "virus");
    EXPECT_STREQ(ExoCategoryName(ExoCategory::Bacteria), "bacteria");
    EXPECT_STREQ(ExoCategoryName(ExoCategory::SpikeIn), "spikein");
    EXPECT_STREQ(ExoCategoryName(ExoCategory::Fungus), "fungus");
    EXPECT_STREQ(ExoCategoryName(ExoCategory::Unknown), "unknown");
}

}  // namespace
