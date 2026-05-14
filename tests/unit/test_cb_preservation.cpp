// LLmap — Unit tests for single-cell barcode preservation.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

#include "singlecell/cb_preservation.h"

namespace llmap::singlecell {
namespace {

// --- TagType conversion tests ---

TEST(TagTypeTest, ToStringRoundTrip) {
    EXPECT_EQ("CB", TagTypeToString(TagType::CB));
    EXPECT_EQ("CR", TagTypeToString(TagType::CR));
    EXPECT_EQ("CY", TagTypeToString(TagType::CY));
    EXPECT_EQ("UB", TagTypeToString(TagType::UB));
    EXPECT_EQ("UR", TagTypeToString(TagType::UR));
    EXPECT_EQ("UY", TagTypeToString(TagType::UY));
    EXPECT_EQ("RG", TagTypeToString(TagType::RG));
    EXPECT_EQ("BC", TagTypeToString(TagType::BC));
    EXPECT_EQ("QT", TagTypeToString(TagType::QT));
}

TEST(TagTypeTest, StringToTagType) {
    EXPECT_EQ(TagType::CB, StringToTagType("CB"));
    EXPECT_EQ(TagType::UB, StringToTagType("UB"));
    EXPECT_EQ(TagType::RG, StringToTagType("RG"));
    EXPECT_EQ(std::nullopt, StringToTagType("XX"));
    EXPECT_EQ(std::nullopt, StringToTagType(""));
}

// --- TagValue tests ---

TEST(TagValueTest, FromString) {
    auto tv = TagValue::FromString("ACGT");
    EXPECT_TRUE(tv.is_string);
    EXPECT_EQ("ACGT", tv.string_value);
}

TEST(TagValueTest, FromInt) {
    auto tv = TagValue::FromInt(42);
    EXPECT_FALSE(tv.is_string);
    EXPECT_EQ(42, tv.int_value);
}

TEST(TagValueTest, Equality) {
    auto s1 = TagValue::FromString("ABC");
    auto s2 = TagValue::FromString("ABC");
    auto s3 = TagValue::FromString("DEF");
    auto i1 = TagValue::FromInt(10);
    auto i2 = TagValue::FromInt(10);
    auto i3 = TagValue::FromInt(20);

    EXPECT_EQ(s1, s2);
    EXPECT_NE(s1, s3);
    EXPECT_EQ(i1, i2);
    EXPECT_NE(i1, i3);
    EXPECT_NE(s1, i1);
}

// --- SingleCellTags tests ---

TEST(SingleCellTagsTest, DefaultConstruction) {
    SingleCellTags tags;
    EXPECT_FALSE(tags.HasCellBarcode());
    EXPECT_FALSE(tags.HasUmi());
    EXPECT_FALSE(tags.HasAnyTags());
}

TEST(SingleCellTagsTest, HasCellBarcode) {
    SingleCellTags tags;
    tags.cell_barcode = "ACGTACGT";
    EXPECT_TRUE(tags.HasCellBarcode());

    SingleCellTags tags2;
    tags2.cell_barcode_raw = "ACGTACGT";
    EXPECT_TRUE(tags2.HasCellBarcode());
}

TEST(SingleCellTagsTest, HasUmi) {
    SingleCellTags tags;
    tags.umi = "AAAA";
    EXPECT_TRUE(tags.HasUmi());

    SingleCellTags tags2;
    tags2.umi_raw = "TTTT";
    EXPECT_TRUE(tags2.HasUmi());
}

TEST(SingleCellTagsTest, GetEffectiveCellBarcode) {
    SingleCellTags tags;
    EXPECT_EQ(std::nullopt, tags.GetEffectiveCellBarcode());

    tags.cell_barcode_raw = "RAW_CB";
    EXPECT_EQ("RAW_CB", tags.GetEffectiveCellBarcode());

    tags.cell_barcode = "CORRECTED_CB";
    EXPECT_EQ("CORRECTED_CB", tags.GetEffectiveCellBarcode());
}

TEST(SingleCellTagsTest, GetEffectiveUmi) {
    SingleCellTags tags;
    EXPECT_EQ(std::nullopt, tags.GetEffectiveUmi());

    tags.umi_raw = "RAW_UMI";
    EXPECT_EQ("RAW_UMI", tags.GetEffectiveUmi());

    tags.umi = "CORRECTED_UMI";
    EXPECT_EQ("CORRECTED_UMI", tags.GetEffectiveUmi());
}

TEST(SingleCellTagsTest, MergeFrom) {
    SingleCellTags base;
    base.cell_barcode = "OLD_CB";
    base.umi = "OLD_UMI";

    SingleCellTags other;
    other.cell_barcode = "NEW_CB";
    other.read_group = "RG1";

    base.MergeFrom(other);
    EXPECT_EQ("NEW_CB", *base.cell_barcode);
    EXPECT_EQ("OLD_UMI", *base.umi);
    EXPECT_EQ("RG1", *base.read_group);
}

TEST(SingleCellTagsTest, Clear) {
    SingleCellTags tags;
    tags.cell_barcode = "CB";
    tags.umi = "UMI";
    tags.custom_tags["XX"] = TagValue::FromString("custom");

    EXPECT_TRUE(tags.HasAnyTags());
    tags.Clear();
    EXPECT_FALSE(tags.HasAnyTags());
}

TEST(SingleCellTagsTest, Equality) {
    SingleCellTags t1, t2;
    EXPECT_EQ(t1, t2);

    t1.cell_barcode = "CB";
    EXPECT_NE(t1, t2);

    t2.cell_barcode = "CB";
    EXPECT_EQ(t1, t2);
}

// --- Extraction tests ---

TEST(ExtractTagsTest, FromAuxStringBasic) {
    auto tags = ExtractTagsFromAux("CB:Z:ACGTACGT\tUB:Z:TTTTTTTT");
    EXPECT_EQ("ACGTACGT", *tags.cell_barcode);
    EXPECT_EQ("TTTTTTTT", *tags.umi);
}

TEST(ExtractTagsTest, FromAuxStringAllTags) {
    auto tags = ExtractTagsFromAux(
        "CB:Z:AAAA\tCR:Z:BBBB\tCY:Z:FFFF\t"
        "UB:Z:CCCC\tUR:Z:DDDD\tUY:Z:GGGG\t"
        "RG:Z:group1\tBC:Z:EEEE\tQT:Z:HHHH");

    EXPECT_EQ("AAAA", *tags.cell_barcode);
    EXPECT_EQ("BBBB", *tags.cell_barcode_raw);
    EXPECT_EQ("FFFF", *tags.cell_barcode_quality);
    EXPECT_EQ("CCCC", *tags.umi);
    EXPECT_EQ("DDDD", *tags.umi_raw);
    EXPECT_EQ("GGGG", *tags.umi_quality);
    EXPECT_EQ("group1", *tags.read_group);
    EXPECT_EQ("EEEE", *tags.sample_barcode);
    EXPECT_EQ("HHHH", *tags.sample_barcode_quality);
}

TEST(ExtractTagsTest, FromAuxStringCustomTags) {
    auto tags = ExtractTagsFromAux("CB:Z:ACGT\tXX:Z:custom\tNM:i:5");
    EXPECT_EQ("ACGT", *tags.cell_barcode);
    EXPECT_TRUE(tags.custom_tags.contains("XX"));
    EXPECT_EQ("custom", tags.custom_tags["XX"].string_value);
    EXPECT_TRUE(tags.custom_tags.contains("NM"));
    EXPECT_EQ(5, tags.custom_tags["NM"].int_value);
}

TEST(ExtractTagsTest, FromAuxStringEmpty) {
    auto tags = ExtractTagsFromAux("");
    EXPECT_FALSE(tags.HasAnyTags());
}

TEST(ExtractTagsTest, FromPairs) {
    std::vector<std::pair<std::string, std::string>> pairs = {
        {"CB", "AAAA"},
        {"UB", "TTTT"},
        {"RG", "group1"}
    };

    auto tags = ExtractTagsFromPairs(pairs);
    EXPECT_EQ("AAAA", *tags.cell_barcode);
    EXPECT_EQ("TTTT", *tags.umi);
    EXPECT_EQ("group1", *tags.read_group);
}

TEST(ExtractTagsTest, FromReadNameUnderscore) {
    BarcodeExtractionConfig config;
    config.expected_barcode_length = 0;
    config.expected_umi_length = 0;

    auto tags = ExtractTagsFromReadName("read_ACGTACGT_TTTTTTTT", config);
    EXPECT_EQ("ACGTACGT", *tags.cell_barcode_raw);
    EXPECT_EQ("TTTTTTTT", *tags.umi_raw);
}

TEST(ExtractTagsTest, FromReadNameColon) {
    BarcodeExtractionConfig config;
    config.expected_barcode_length = 0;
    config.expected_umi_length = 0;

    auto tags = ExtractTagsFromReadName("sample:ACGT:TTTT", config);
    EXPECT_EQ("ACGT", *tags.cell_barcode_raw);
    EXPECT_EQ("TTTT", *tags.umi_raw);
}

TEST(ExtractTagsTest, FromReadNameNoBarcode) {
    BarcodeExtractionConfig config;
    auto tags = ExtractTagsFromReadName("simple_read_name", config);
    EXPECT_FALSE(tags.HasCellBarcode());
}

// --- Format tests ---

TEST(FormatTagsTest, AsAuxBasic) {
    SingleCellTags tags;
    tags.cell_barcode = "ACGT";
    tags.umi = "TTTT";

    auto aux = FormatTagsAsAux(tags);
    EXPECT_NE(std::string::npos, aux.find("CB:Z:ACGT"));
    EXPECT_NE(std::string::npos, aux.find("UB:Z:TTTT"));
}

TEST(FormatTagsTest, AsAuxEmpty) {
    SingleCellTags tags;
    auto aux = FormatTagsAsAux(tags);
    EXPECT_TRUE(aux.empty());
}

TEST(FormatTagsTest, AsPairs) {
    SingleCellTags tags;
    tags.cell_barcode = "ACGT";
    tags.umi = "TTTT";
    tags.custom_tags["XX"] = TagValue::FromString("val");

    auto pairs = FormatTagsAsPairs(tags);
    EXPECT_EQ(3u, pairs.size());

    bool found_cb = false, found_ub = false, found_xx = false;
    for (const auto& [tag, val] : pairs) {
        if (tag == "CB" && val == "ACGT") found_cb = true;
        if (tag == "UB" && val == "TTTT") found_ub = true;
        if (tag == "XX" && val == "val") found_xx = true;
    }
    EXPECT_TRUE(found_cb);
    EXPECT_TRUE(found_ub);
    EXPECT_TRUE(found_xx);
}

TEST(FormatTagsTest, RoundTrip) {
    SingleCellTags original;
    original.cell_barcode = "AAAA";
    original.cell_barcode_raw = "BBBB";
    original.umi = "CCCC";
    original.read_group = "RG1";

    auto aux = FormatTagsAsAux(original);
    auto reconstructed = ExtractTagsFromAux(aux);

    EXPECT_EQ(*original.cell_barcode, *reconstructed.cell_barcode);
    EXPECT_EQ(*original.cell_barcode_raw, *reconstructed.cell_barcode_raw);
    EXPECT_EQ(*original.umi, *reconstructed.umi);
    EXPECT_EQ(*original.read_group, *reconstructed.read_group);
}

// --- Validation tests ---

TEST(ValidateBarcodeTest, ValidBarcodes) {
    BarcodeExtractionConfig config;
    config.expected_barcode_length = 16;

    EXPECT_TRUE(ValidateBarcode("ACGTACGTACGTACGT", config));
    EXPECT_TRUE(ValidateBarcode("acgtacgtacgtacgt", config));
}

TEST(ValidateBarcodeTest, WrongLength) {
    BarcodeExtractionConfig config;
    config.expected_barcode_length = 16;

    EXPECT_FALSE(ValidateBarcode("ACGT", config));

    config.allow_partial_barcode = true;
    EXPECT_TRUE(ValidateBarcode("ACGT", config));
}

TEST(ValidateBarcodeTest, WithN) {
    BarcodeExtractionConfig config;
    config.expected_barcode_length = 4;
    config.allow_n_in_barcode = false;

    EXPECT_FALSE(ValidateBarcode("ACNT", config));

    config.allow_n_in_barcode = true;
    EXPECT_TRUE(ValidateBarcode("ACNT", config));
}

TEST(ValidateBarcodeTest, InvalidCharacters) {
    BarcodeExtractionConfig config;
    config.expected_barcode_length = 4;

    EXPECT_FALSE(ValidateBarcode("ACGX", config));
    EXPECT_FALSE(ValidateBarcode("1234", config));
}

TEST(ValidateUmiTest, ValidUmis) {
    BarcodeExtractionConfig config;
    config.expected_umi_length = 12;

    EXPECT_TRUE(ValidateUmi("ACGTACGTACGT", config));
}

TEST(ValidateUmiTest, WrongLength) {
    BarcodeExtractionConfig config;
    config.expected_umi_length = 12;

    EXPECT_FALSE(ValidateUmi("ACGT", config));

    config.allow_partial_barcode = true;
    EXPECT_TRUE(ValidateUmi("ACGT", config));
}

// --- Whitelist tests ---

TEST(CellBarcodeWhitelistTest, AddAndContains) {
    CellBarcodeWhitelist wl;
    EXPECT_TRUE(wl.Empty());

    wl.Add("ACGT");
    wl.Add("TTTT");
    EXPECT_EQ(2u, wl.Size());
    EXPECT_FALSE(wl.Empty());

    EXPECT_TRUE(wl.Contains("ACGT"));
    EXPECT_TRUE(wl.Contains("TTTT"));
    EXPECT_FALSE(wl.Contains("GGGG"));
}

TEST(CellBarcodeWhitelistTest, FindNearestExact) {
    CellBarcodeWhitelist wl;
    wl.Add("ACGT");

    auto result = wl.FindNearest("ACGT", 1);
    EXPECT_EQ("ACGT", *result);
}

TEST(CellBarcodeWhitelistTest, FindNearestOneMismatch) {
    CellBarcodeWhitelist wl;
    wl.Add("ACGT");

    auto result = wl.FindNearest("ACGA", 1);  // One mismatch
    EXPECT_EQ("ACGT", *result);

    result = wl.FindNearest("AAAA", 1);  // Three mismatches
    EXPECT_EQ(std::nullopt, result);
}

TEST(CellBarcodeWhitelistTest, FindNearestZeroDistance) {
    CellBarcodeWhitelist wl;
    wl.Add("ACGT");

    auto result = wl.FindNearest("ACGA", 0);  // Exact only
    EXPECT_EQ(std::nullopt, result);

    result = wl.FindNearest("ACGT", 0);
    EXPECT_EQ("ACGT", *result);
}

TEST(CellBarcodeWhitelistTest, LoadFromFile) {
    std::filesystem::path temp_path =
        std::filesystem::temp_directory_path() / "test_whitelist.txt";

    {
        std::ofstream out(temp_path);
        out << "AAAA\n";
        out << "BBBB\n";
        out << "CCCC\tignored\n";  // Tab-separated, take first column
    }

    auto wl = CellBarcodeWhitelist::LoadFromFile(temp_path);
    EXPECT_TRUE(wl.has_value());
    EXPECT_EQ(3u, wl->Size());
    EXPECT_TRUE(wl->Contains("AAAA"));
    EXPECT_TRUE(wl->Contains("CCCC"));
    EXPECT_FALSE(wl->Contains("ignored"));

    std::filesystem::remove(temp_path);
}

TEST(CellBarcodeWhitelistTest, LoadFromNonexistentFile) {
    auto wl = CellBarcodeWhitelist::LoadFromFile("/nonexistent/path.txt");
    EXPECT_FALSE(wl.has_value());
}

}  // namespace
}  // namespace llmap::singlecell
