// Unit tests for configuration file support.

#include <gtest/gtest.h>

#include "core/config.h"

#include <cstdlib>
#include <fstream>
#include <filesystem>

namespace llmap {
namespace {

class ConfigTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() /
                    ("llmap_config_test_" + std::to_string(getpid()));
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::error_code ec;
        std::filesystem::remove_all(test_dir_, ec);
    }

    std::filesystem::path test_dir_;
};

TEST(ConfigValueTest, AsBoolTrue) {
    ConfigValue v;
    v.value = "true";
    EXPECT_TRUE(v.AsBool());

    v.value = "yes";
    EXPECT_TRUE(v.AsBool());

    v.value = "on";
    EXPECT_TRUE(v.AsBool());

    v.value = "1";
    EXPECT_TRUE(v.AsBool());
}

TEST(ConfigValueTest, AsBoolFalse) {
    ConfigValue v;
    v.value = "false";
    EXPECT_FALSE(v.AsBool());

    v.value = "no";
    EXPECT_FALSE(v.AsBool());

    v.value = "off";
    EXPECT_FALSE(v.AsBool());

    v.value = "0";
    EXPECT_FALSE(v.AsBool());
}

TEST(ConfigValueTest, AsBoolDefault) {
    ConfigValue v;
    v.value = "maybe";
    EXPECT_TRUE(v.AsBool(true));
    EXPECT_FALSE(v.AsBool(false));
}

TEST(ConfigValueTest, AsInt) {
    ConfigValue v;
    v.value = "42";
    EXPECT_EQ(v.AsInt(), 42);

    v.value = "-100";
    EXPECT_EQ(v.AsInt(), -100);

    v.value = "  123  ";
    EXPECT_EQ(v.AsInt(), 123);
}

TEST(ConfigValueTest, AsIntDefault) {
    ConfigValue v;
    v.value = "not_a_number";
    EXPECT_EQ(v.AsInt(99), 99);
}

TEST(ConfigValueTest, AsDouble) {
    ConfigValue v;
    v.value = "3.14";
    EXPECT_DOUBLE_EQ(v.AsDouble(), 3.14);

    v.value = "0.5";
    EXPECT_DOUBLE_EQ(v.AsDouble(), 0.5);
}

TEST(ConfigValueTest, AsDoubleDefault) {
    ConfigValue v;
    v.value = "invalid";
    EXPECT_DOUBLE_EQ(v.AsDouble(1.5), 1.5);
}

TEST(ConfigValueTest, AsString) {
    ConfigValue v;
    v.value = "hello world";
    EXPECT_EQ(v.AsString(), "hello world");
}

TEST(ConfigParserTest, ParseEmptyContent) {
    ConfigParser parser;
    auto result = parser.Parse("", "test.toml");
    ASSERT_TRUE(result.ok());
    auto config = result.value();
    EXPECT_EQ(config.align.kmer_size, 15);
}

TEST(ConfigParserTest, ParseComments) {
    ConfigParser parser;
    auto result = parser.Parse("# This is a comment\n# Another comment", "test.toml");
    ASSERT_TRUE(result.ok());
    EXPECT_EQ(parser.Values().size(), 0);
}

TEST(ConfigParserTest, ParseSimpleKeyValue) {
    ConfigParser parser;
    auto result = parser.Parse("reference = /path/to/ref.fa", "test.toml");
    ASSERT_TRUE(result.ok());

    auto val = parser.Get("", "reference");
    ASSERT_TRUE(val.has_value());
    EXPECT_EQ(val->value, "/path/to/ref.fa");
}

TEST(ConfigParserTest, ParseQuotedValue) {
    ConfigParser parser;
    auto result = parser.Parse("reference = \"/path/with spaces/ref.fa\"", "test.toml");
    ASSERT_TRUE(result.ok());

    auto val = parser.Get("", "reference");
    ASSERT_TRUE(val.has_value());
    EXPECT_EQ(val->value, "/path/with spaces/ref.fa");
}

TEST(ConfigParserTest, ParseSection) {
    ConfigParser parser;
    auto result = parser.Parse("[align]\nthreads = 8\nkmer_size = 17", "test.toml");
    ASSERT_TRUE(result.ok());

    auto threads = parser.Get("align", "threads");
    ASSERT_TRUE(threads.has_value());
    EXPECT_EQ(threads->AsInt(), 8);

    auto kmer = parser.Get("align", "kmer_size");
    ASSERT_TRUE(kmer.has_value());
    EXPECT_EQ(kmer->AsInt(), 17);
}

TEST(ConfigParserTest, ParseMultipleSections) {
    ConfigParser parser;
    std::string content = R"(
[align]
threads = 4
kmer_size = 19

[llm]
enabled = true
model = "claude-opus-4-20250514"
temperature = 0.5
)";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    EXPECT_EQ(parser.Get("align", "threads")->AsInt(), 4);
    EXPECT_EQ(parser.Get("align", "kmer_size")->AsInt(), 19);
    EXPECT_TRUE(parser.Get("llm", "enabled")->AsBool());
    EXPECT_EQ(parser.Get("llm", "model")->AsString(), "claude-opus-4-20250514");
    EXPECT_DOUBLE_EQ(parser.Get("llm", "temperature")->AsDouble(), 0.5);
}

TEST(ConfigParserTest, ParseUnclosedSection) {
    ConfigParser parser;
    auto result = parser.Parse("[align", "test.toml");
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code(), ErrorCode::kParseInvalidFormat);
}

TEST(ConfigParserTest, ParseMissingEquals) {
    ConfigParser parser;
    auto result = parser.Parse("invalid_line", "test.toml");
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code(), ErrorCode::kParseInvalidFormat);
}

TEST(ConfigParserTest, ParseEmptyKey) {
    ConfigParser parser;
    auto result = parser.Parse(" = value", "test.toml");
    ASSERT_FALSE(result.ok());
}

TEST(ConfigParserTest, PopulatesAlignConfig) {
    ConfigParser parser;
    std::string content = R"(
[align]
threads = 16
kmer_size = 21
window_size = 15
min_chain_score = 100
max_gaps_in_chain = 10
min_identity = 0.95
min_mapping_quality = 20.0
output_format = "sam"
include_unmapped = true
include_secondary = false
)";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    auto config = result.value();
    EXPECT_EQ(config.align.threads, 16);
    EXPECT_EQ(config.align.kmer_size, 21);
    EXPECT_EQ(config.align.window_size, 15);
    EXPECT_EQ(config.align.min_chain_score, 100);
    EXPECT_EQ(config.align.max_gaps_in_chain, 10);
    EXPECT_DOUBLE_EQ(config.align.min_identity, 0.95);
    EXPECT_DOUBLE_EQ(config.align.min_mapping_quality, 20.0);
    EXPECT_EQ(config.align.output_format, "sam");
    EXPECT_TRUE(config.align.include_unmapped);
    EXPECT_FALSE(config.align.include_secondary);
}

TEST(ConfigParserTest, PopulatesLlmConfig) {
    ConfigParser parser;
    std::string content = R"(
[llm]
api_key = "sk-test-key"
model = "claude-opus-4-20250514"
temperature = 0.3
max_tokens = 8192
stall_threshold = 0.05
work_dir = "/custom/llm/dir"
enabled = true
)";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    auto config = result.value();
    EXPECT_EQ(config.llm.api_key, "sk-test-key");
    EXPECT_EQ(config.llm.model, "claude-opus-4-20250514");
    EXPECT_DOUBLE_EQ(config.llm.temperature, 0.3);
    EXPECT_EQ(config.llm.max_tokens, 8192);
    EXPECT_DOUBLE_EQ(config.llm.stall_threshold, 0.05);
    EXPECT_EQ(config.llm.work_dir, "/custom/llm/dir");
    EXPECT_TRUE(config.llm.enabled);
}

TEST(ConfigParserTest, PopulatesSingleCellConfig) {
    ConfigParser parser;
    std::string content = R"TOML(
[singlecell]
cb_tag = "CR"
umi_tag = "UR"
cb_pattern = "^([ACGT]{16})"
cb_whitelist_file = "/path/to/whitelist.txt"
platform = "parse"
min_reads_per_cell = 500
min_confidence = 0.8
)TOML";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    auto config = result.value();
    EXPECT_EQ(config.singlecell.cb_tag, "CR");
    EXPECT_EQ(config.singlecell.umi_tag, "UR");
    EXPECT_EQ(config.singlecell.cb_pattern, "^([ACGT]{16})");
    EXPECT_EQ(config.singlecell.cb_whitelist_file, "/path/to/whitelist.txt");
    EXPECT_EQ(config.singlecell.platform, "parse");
    EXPECT_EQ(config.singlecell.min_reads_per_cell, 500);
    EXPECT_DOUBLE_EQ(config.singlecell.min_confidence, 0.8);
}

TEST(ConfigParserTest, PopulatesPsvConfig) {
    ConfigParser parser;
    std::string content = R"(
[psv]
catalog_file = "/path/to/psv.bed"
weight = 0.7
min_posterior = 0.95
enabled = true
)";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    auto config = result.value();
    EXPECT_EQ(config.psv.catalog_file, "/path/to/psv.bed");
    EXPECT_DOUBLE_EQ(config.psv.weight, 0.7);
    EXPECT_DOUBLE_EQ(config.psv.min_posterior, 0.95);
    EXPECT_TRUE(config.psv.enabled);
}

TEST(ConfigParserTest, PopulatesLoggingConfig) {
    ConfigParser parser;
    std::string content = R"(
[logging]
level = "debug"
format = "json"
file = "/var/log/llmap.log"
)";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    auto config = result.value();
    EXPECT_EQ(config.logging.level, LogLevel::kDebug);
    EXPECT_EQ(config.logging.format, LogFormat::kJson);
    EXPECT_EQ(config.logging.file, "/var/log/llmap.log");
}

TEST(ConfigParserTest, PopulatesGeneralConfig) {
    ConfigParser parser;
    std::string content = R"(
reference = "/path/to/reference.fa"
index = "/path/to/index.idx"
temp_dir = "/scratch/llmap"
)";
    auto result = parser.Parse(content, "test.toml");
    ASSERT_TRUE(result.ok());

    auto config = result.value();
    EXPECT_EQ(config.reference_path, "/path/to/reference.fa");
    EXPECT_EQ(config.index_path, "/path/to/index.idx");
    EXPECT_EQ(config.temp_dir, "/scratch/llmap");
}

TEST_F(ConfigTest, ParseFileNotFound) {
    ConfigParser parser;
    auto result = parser.ParseFile(test_dir_ / "nonexistent.toml");
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code(), ErrorCode::kIoFileNotFound);
}

TEST_F(ConfigTest, ParseFileValid) {
    auto config_path = test_dir_ / "test.toml";
    {
        std::ofstream f(config_path);
        f << "[align]\nthreads = 32\n";
    }

    ConfigParser parser;
    auto result = parser.ParseFile(config_path);
    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.value().align.threads, 32);
}

TEST(ConfigDefaultsTest, HasSensibleDefaults) {
    auto config = LLmapConfig::Defaults();

    EXPECT_EQ(config.align.threads, 0);
    EXPECT_EQ(config.align.kmer_size, 15);
    EXPECT_EQ(config.align.window_size, 10);
    EXPECT_EQ(config.align.output_format, "bam");
    EXPECT_FALSE(config.align.include_unmapped);
    EXPECT_TRUE(config.align.include_secondary);

    EXPECT_FALSE(config.llm.enabled);
    EXPECT_TRUE(config.llm.api_key.empty());

    EXPECT_EQ(config.singlecell.cb_tag, "CB");
    EXPECT_EQ(config.singlecell.platform, "10x-v3");

    EXPECT_FALSE(config.psv.enabled);

    EXPECT_EQ(config.logging.level, LogLevel::kInfo);
    EXPECT_EQ(config.logging.format, LogFormat::kText);
}

TEST(ConfigSearchPathsTest, ReturnsMultiplePaths) {
    auto paths = GetConfigSearchPaths();
    EXPECT_GE(paths.size(), 3);

    bool has_local = false;
    bool has_etc = false;
    for (const auto& p : paths) {
        if (p == "./llmap.toml" || p == "./.llmap.toml") has_local = true;
        if (p == "/etc/llmap/config.toml") has_etc = true;
    }
    EXPECT_TRUE(has_local);
    EXPECT_TRUE(has_etc);
}

TEST(ConfigOverridesTest, ApplyOverrides) {
    auto config = LLmapConfig::Defaults();

    std::vector<ConfigOverride> overrides = {
        {"threads", "64"},
        {"kmer_size", "21"},
        {"min_identity", "0.85"},
        {"reference", "/new/ref.fa"},
        {"llm.enabled", "true"},
        {"log_level", "debug"}
    };

    ApplyOverrides(config, overrides);

    EXPECT_EQ(config.align.threads, 64);
    EXPECT_EQ(config.align.kmer_size, 21);
    EXPECT_DOUBLE_EQ(config.align.min_identity, 0.85);
    EXPECT_EQ(config.reference_path, "/new/ref.fa");
    EXPECT_TRUE(config.llm.enabled);
    EXPECT_EQ(config.logging.level, LogLevel::kDebug);
}

TEST(ConfigValidationTest, ValidConfigPasses) {
    auto config = LLmapConfig::Defaults();
    auto result = ValidateConfig(config);
    EXPECT_TRUE(result.ok());
}

TEST(ConfigValidationTest, InvalidKmerSize) {
    auto config = LLmapConfig::Defaults();
    config.align.kmer_size = 3;
    auto result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());

    config.align.kmer_size = 35;
    result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());
}

TEST(ConfigValidationTest, InvalidMinIdentity) {
    auto config = LLmapConfig::Defaults();
    config.align.min_identity = -0.1;
    auto result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());

    config.align.min_identity = 1.5;
    result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());
}

TEST(ConfigValidationTest, InvalidLlmTemperature) {
    auto config = LLmapConfig::Defaults();
    config.llm.temperature = 3.0;
    auto result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());
}

TEST(ConfigValidationTest, LlmEnabledWithoutApiKey) {
    auto config = LLmapConfig::Defaults();
    config.llm.enabled = true;
    config.llm.api_key = "";
    auto result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());
    EXPECT_EQ(result.error().code(), ErrorCode::kConfigMissingRequired);
}

TEST(ConfigValidationTest, InvalidPsvWeight) {
    auto config = LLmapConfig::Defaults();
    config.psv.weight = 1.5;
    auto result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());
}

TEST(ConfigValidationTest, InvalidSingleCellConfidence) {
    auto config = LLmapConfig::Defaults();
    config.singlecell.min_confidence = -0.5;
    auto result = ValidateConfig(config);
    EXPECT_FALSE(result.ok());
}

TEST(ConfigToTomlTest, GeneratesValidToml) {
    auto config = LLmapConfig::Defaults();
    config.reference_path = "/path/to/ref.fa";
    config.align.threads = 8;
    config.align.kmer_size = 17;
    config.llm.enabled = true;
    config.llm.api_key = "test-key";

    std::string toml = ConfigToToml(config);

    EXPECT_NE(toml.find("reference = \"/path/to/ref.fa\""), std::string::npos);
    EXPECT_NE(toml.find("[align]"), std::string::npos);
    EXPECT_NE(toml.find("threads = 8"), std::string::npos);
    EXPECT_NE(toml.find("kmer_size = 17"), std::string::npos);
    EXPECT_NE(toml.find("[llm]"), std::string::npos);
    EXPECT_NE(toml.find("enabled = true"), std::string::npos);
}

TEST(ConfigToTomlTest, RoundTrip) {
    auto original = LLmapConfig::Defaults();
    original.reference_path = "/test/ref.fa";
    original.align.threads = 16;
    original.align.kmer_size = 19;
    original.align.min_identity = 0.92;
    original.llm.enabled = true;
    original.llm.api_key = "sk-test";
    original.llm.temperature = 0.5;
    original.singlecell.platform = "parse";
    original.psv.enabled = true;
    original.psv.weight = 0.6;
    original.logging.level = LogLevel::kDebug;

    std::string toml = ConfigToToml(original);

    ConfigParser parser;
    auto result = parser.Parse(toml, "generated.toml");
    ASSERT_TRUE(result.ok());

    auto parsed = result.value();
    EXPECT_EQ(parsed.reference_path, original.reference_path);
    EXPECT_EQ(parsed.align.threads, original.align.threads);
    EXPECT_EQ(parsed.align.kmer_size, original.align.kmer_size);
    EXPECT_DOUBLE_EQ(parsed.align.min_identity, original.align.min_identity);
    EXPECT_EQ(parsed.llm.enabled, original.llm.enabled);
    EXPECT_DOUBLE_EQ(parsed.llm.temperature, original.llm.temperature);
    EXPECT_EQ(parsed.singlecell.platform, original.singlecell.platform);
    EXPECT_EQ(parsed.psv.enabled, original.psv.enabled);
    EXPECT_DOUBLE_EQ(parsed.psv.weight, original.psv.weight);
    EXPECT_EQ(parsed.logging.level, original.logging.level);
}

TEST(LoadConfigTest, ReturnsDefaultsWhenNoFile) {
    auto result = LoadConfig();
    ASSERT_TRUE(result.ok());
    auto config = result.value();
    EXPECT_EQ(config.align.kmer_size, 15);
}

TEST_F(ConfigTest, LoadConfigFromFile) {
    auto config_path = test_dir_ / "custom.toml";
    {
        std::ofstream f(config_path);
        f << "[align]\nkmer_size = 23\n";
    }

    auto result = LoadConfigFromFile(config_path);
    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.value().align.kmer_size, 23);
    EXPECT_EQ(result.value().source_file, config_path);
}

}  // namespace
}  // namespace llmap
