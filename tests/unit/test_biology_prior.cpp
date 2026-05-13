#include <gtest/gtest.h>

#include "claude_agent/biology_prior.h"

#include <filesystem>
#include <fstream>

using namespace llmap::claude_agent;

class BiologyPriorTest : public ::testing::Test {
protected:
    std::filesystem::path temp_dir_;

    void SetUp() override {
        temp_dir_ = std::filesystem::temp_directory_path() /
                    ("biology_prior_test_" + std::to_string(
                        std::hash<std::thread::id>{}(std::this_thread::get_id())));
        std::filesystem::create_directories(temp_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(temp_dir_);
    }
};

TEST_F(BiologyPriorTest, SerializeEmptyPrior) {
    BiologyPrior prior;
    prior.version = "1.0";
    prior.reference_sha256 = "abc123";

    std::string json = SerializeBiologyPrior(prior);

    EXPECT_NE(json.find("\"version\": \"1.0\""), std::string::npos);
    EXPECT_NE(json.find("\"reference_sha256\": \"abc123\""), std::string::npos);
    EXPECT_NE(json.find("\"buckets\": {"), std::string::npos);
    EXPECT_NE(json.find("\"regional_overrides\": {"), std::string::npos);
}

TEST_F(BiologyPriorTest, SerializePriorWithBucket) {
    BiologyPrior prior;
    prior.version = "1.0";
    prior.reference_sha256 = "abc123";

    BucketAnnotation ann;
    ann.bucket_id = 12345;
    ann.level = "L2";
    ann.annotation = "IGH-Constants C-region";
    ann.prior_weight = 1.2;
    ann.paralog_partner_bucket = 12346;
    ann.expected_coverage_multiplier = 2.0;
    ann.claude_rationale = "Known duplication in HPRC samples";
    prior.buckets[12345] = ann;

    std::string json = SerializeBiologyPrior(prior);

    EXPECT_NE(json.find("\"12345\""), std::string::npos);
    EXPECT_NE(json.find("\"level\": \"L2\""), std::string::npos);
    EXPECT_NE(json.find("\"annotation\": \"IGH-Constants C-region\""), std::string::npos);
    EXPECT_NE(json.find("\"prior_weight\": 1.2"), std::string::npos);
    EXPECT_NE(json.find("\"paralog_partner_bucket\": 12346"), std::string::npos);
    EXPECT_NE(json.find("\"expected_coverage_multiplier\": 2"), std::string::npos);
}

TEST_F(BiologyPriorTest, SerializePriorWithRegionalOverride) {
    BiologyPrior prior;
    prior.version = "1.0";

    RegionalOverride ovr;
    ovr.region = "chr14:105M-107M";
    ovr.sub_bucket_granularity_kb = 10;
    ovr.max_iter = 25;
    ovr.convergence_threshold = 0.95;
    prior.regional_overrides["chr14:105M-107M"] = ovr;

    std::string json = SerializeBiologyPrior(prior);

    EXPECT_NE(json.find("\"chr14:105M-107M\""), std::string::npos);
    EXPECT_NE(json.find("\"sub_bucket_granularity_kb\": 10"), std::string::npos);
    EXPECT_NE(json.find("\"max_iter\": 25"), std::string::npos);
    EXPECT_NE(json.find("\"convergence_threshold\": 0.95"), std::string::npos);
}

TEST_F(BiologyPriorTest, SerializeMultipleBuckets) {
    BiologyPrior prior;
    prior.version = "1.0";

    BucketAnnotation ann1;
    ann1.bucket_id = 100;
    ann1.level = "L1";
    ann1.annotation = "Region A";
    prior.buckets[100] = ann1;

    BucketAnnotation ann2;
    ann2.bucket_id = 200;
    ann2.level = "L2";
    ann2.annotation = "Region B";
    prior.buckets[200] = ann2;

    std::string json = SerializeBiologyPrior(prior);

    EXPECT_NE(json.find("\"100\""), std::string::npos);
    EXPECT_NE(json.find("\"200\""), std::string::npos);
    EXPECT_NE(json.find("\"Region A\""), std::string::npos);
    EXPECT_NE(json.find("\"Region B\""), std::string::npos);
}

TEST_F(BiologyPriorTest, WriteBiologyPriorToFile) {
    BiologyPrior prior;
    prior.version = "1.0";
    prior.reference_sha256 = "test_sha";

    auto output_path = temp_dir_ / "biology_prior.json";
    EXPECT_TRUE(WriteBiologyPrior(prior, output_path));
    EXPECT_TRUE(std::filesystem::exists(output_path));

    std::ifstream ifs(output_path);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"version\": \"1.0\""), std::string::npos);
}

TEST_F(BiologyPriorTest, WriteBiologyPriorToInvalidPath) {
    BiologyPrior prior;
    prior.version = "1.0";

    auto output_path = "/nonexistent/dir/biology_prior.json";
    EXPECT_FALSE(WriteBiologyPrior(prior, output_path));
}

TEST_F(BiologyPriorTest, ReadBiologyPriorNonexistent) {
    auto result = ReadBiologyPrior("/nonexistent/file.json");
    EXPECT_FALSE(result.has_value());
}

TEST_F(BiologyPriorTest, EscapeSpecialCharacters) {
    BiologyPrior prior;
    prior.version = "1.0";

    BucketAnnotation ann;
    ann.bucket_id = 1;
    ann.level = "L1";
    ann.annotation = "Test with \"quotes\" and \\backslash";
    ann.claude_rationale = "Line1\nLine2\tTabbed";
    prior.buckets[1] = ann;

    std::string json = SerializeBiologyPrior(prior);

    EXPECT_NE(json.find("\\\"quotes\\\""), std::string::npos);
    EXPECT_NE(json.find("\\\\backslash"), std::string::npos);
    EXPECT_NE(json.find("\\n"), std::string::npos);
    EXPECT_NE(json.find("\\t"), std::string::npos);
}

TEST_F(BiologyPriorTest, SerializeSampleParams) {
    SampleParams params;
    params.preset = "iso-seq";
    params.convergence_threshold = 0.95;
    params.max_iterations = 30;
    params.expected_coverage_profile = "B-cell";
    params.foundation_model = "evo-1.5b";
    params.region_adjustments["chr14"] = 1.5;
    params.region_adjustments["chr22"] = 0.8;

    std::string json = SerializeSampleParams(params);

    EXPECT_NE(json.find("\"preset\": \"iso-seq\""), std::string::npos);
    EXPECT_NE(json.find("\"convergence_threshold\": 0.95"), std::string::npos);
    EXPECT_NE(json.find("\"max_iterations\": 30"), std::string::npos);
    EXPECT_NE(json.find("\"expected_coverage_profile\": \"B-cell\""), std::string::npos);
    EXPECT_NE(json.find("\"foundation_model\": \"evo-1.5b\""), std::string::npos);
    EXPECT_NE(json.find("\"chr14\": 1.5"), std::string::npos);
}

TEST_F(BiologyPriorTest, WriteSampleParamsToFile) {
    SampleParams params;
    params.preset = "hifi";
    params.max_iterations = 20;

    auto output_path = temp_dir_ / "sample_params.json";
    EXPECT_TRUE(WriteSampleParams(params, output_path));
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST_F(BiologyPriorTest, WriteSampleParamsToInvalidPath) {
    SampleParams params;

    auto output_path = "/nonexistent/dir/sample_params.json";
    EXPECT_FALSE(WriteSampleParams(params, output_path));
}

TEST_F(BiologyPriorTest, ReadSampleParamsNonexistent) {
    auto result = ReadSampleParams("/nonexistent/file.json");
    EXPECT_FALSE(result.has_value());
}

TEST_F(BiologyPriorTest, DeserializeBiologyPriorStub) {
    // Current implementation returns empty prior - just verify it doesn't crash
    auto prior = DeserializeBiologyPrior("{}");
    EXPECT_TRUE(prior.empty());
}

TEST_F(BiologyPriorTest, DeserializeSampleParamsStub) {
    // Current implementation returns default params - verify it doesn't crash
    auto params = DeserializeSampleParams("{}");
    EXPECT_EQ("hifi", params.preset);  // default
}

TEST_F(BiologyPriorTest, RoundTripBiologyPriorFile) {
    BiologyPrior prior;
    prior.version = "1.0";
    prior.reference_sha256 = "roundtrip_test";

    BucketAnnotation ann;
    ann.bucket_id = 999;
    ann.level = "L3";
    ann.annotation = "Test annotation";
    prior.buckets[999] = ann;

    auto output_path = temp_dir_ / "roundtrip.json";
    EXPECT_TRUE(WriteBiologyPrior(prior, output_path));

    // Read back - currently returns empty, but file should exist
    EXPECT_TRUE(std::filesystem::exists(output_path));
    auto file_size = std::filesystem::file_size(output_path);
    EXPECT_GT(file_size, 0u);
}

TEST_F(BiologyPriorTest, BucketWithoutParalogPartner) {
    BiologyPrior prior;
    prior.version = "1.0";

    BucketAnnotation ann;
    ann.bucket_id = 1;
    ann.level = "L1";
    ann.annotation = "No paralog";
    // paralog_partner_bucket is std::nullopt
    prior.buckets[1] = ann;

    std::string json = SerializeBiologyPrior(prior);

    // Should not contain paralog_partner_bucket field
    EXPECT_EQ(json.find("paralog_partner_bucket"), std::string::npos);
}
