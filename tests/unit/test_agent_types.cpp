#include <gtest/gtest.h>

#include "claude_agent/agent_types.h"

using namespace llmap::claude_agent;

class AgentTypesTest : public ::testing::Test {};

TEST_F(AgentTypesTest, ParseAgentModeOff) {
    EXPECT_EQ(AgentMode::Off, *ParseAgentMode("off"));
    EXPECT_EQ(AgentMode::Off, *ParseAgentMode("OFF"));
    EXPECT_EQ(AgentMode::Off, *ParseAgentMode("0"));
}

TEST_F(AgentTypesTest, ParseAgentModeIndexOnly) {
    EXPECT_EQ(AgentMode::IndexOnly, *ParseAgentMode("index-only"));
    EXPECT_EQ(AgentMode::IndexOnly, *ParseAgentMode("indexonly"));
    EXPECT_EQ(AgentMode::IndexOnly, *ParseAgentMode("INDEX-ONLY"));
    EXPECT_EQ(AgentMode::IndexOnly, *ParseAgentMode("1"));
}

TEST_F(AgentTypesTest, ParseAgentModeSampleAware) {
    EXPECT_EQ(AgentMode::SampleAware, *ParseAgentMode("sample-aware"));
    EXPECT_EQ(AgentMode::SampleAware, *ParseAgentMode("sampleaware"));
    EXPECT_EQ(AgentMode::SampleAware, *ParseAgentMode("SAMPLE-AWARE"));
    EXPECT_EQ(AgentMode::SampleAware, *ParseAgentMode("2"));
}

TEST_F(AgentTypesTest, ParseAgentModeSelfHealing) {
    EXPECT_EQ(AgentMode::SelfHealing, *ParseAgentMode("self-healing"));
    EXPECT_EQ(AgentMode::SelfHealing, *ParseAgentMode("selfhealing"));
    EXPECT_EQ(AgentMode::SelfHealing, *ParseAgentMode("SELF-HEALING"));
    EXPECT_EQ(AgentMode::SelfHealing, *ParseAgentMode("3"));
}

TEST_F(AgentTypesTest, ParseAgentModeResearch) {
    EXPECT_EQ(AgentMode::Research, *ParseAgentMode("research"));
    EXPECT_EQ(AgentMode::Research, *ParseAgentMode("RESEARCH"));
    EXPECT_EQ(AgentMode::Research, *ParseAgentMode("4"));
}

TEST_F(AgentTypesTest, ParseAgentModeInvalid) {
    EXPECT_FALSE(ParseAgentMode("invalid").has_value());
    EXPECT_FALSE(ParseAgentMode("").has_value());
    EXPECT_FALSE(ParseAgentMode("5").has_value());
    EXPECT_FALSE(ParseAgentMode("sample_aware").has_value());
}

TEST_F(AgentTypesTest, SessionTypeName) {
    EXPECT_STREQ("index-build", SessionTypeName(SessionType::IndexBuild));
    EXPECT_STREQ("sample-init", SessionTypeName(SessionType::SampleInit));
    EXPECT_STREQ("diagnostic", SessionTypeName(SessionType::Diagnostic));
    EXPECT_STREQ("reporter", SessionTypeName(SessionType::Reporter));
}

TEST_F(AgentTypesTest, AgentModeName) {
    EXPECT_STREQ("off", AgentModeName(AgentMode::Off));
    EXPECT_STREQ("index-only", AgentModeName(AgentMode::IndexOnly));
    EXPECT_STREQ("sample-aware", AgentModeName(AgentMode::SampleAware));
    EXPECT_STREQ("self-healing", AgentModeName(AgentMode::SelfHealing));
    EXPECT_STREQ("research", AgentModeName(AgentMode::Research));
}

TEST_F(AgentTypesTest, AgentConfigDefaults) {
    AgentConfig config;
    EXPECT_EQ(AgentMode::SampleAware, config.mode);
    EXPECT_TRUE(config.api_key.empty());
    EXPECT_EQ("claude-sonnet-4-20250514", config.model);
    EXPECT_EQ(8192u, config.max_tokens);
    EXPECT_DOUBLE_EQ(0.0, config.temperature);
    EXPECT_EQ(50u, config.rate_limit_rpm);
    EXPECT_EQ(3u, config.max_retries);
    EXPECT_FALSE(config.cache_dir.has_value());
}

TEST_F(AgentTypesTest, AgentResultDefaults) {
    AgentResult result;
    EXPECT_EQ(AgentStatus::Pending, result.status);
    EXPECT_TRUE(result.session_id.empty());
    EXPECT_EQ(0u, result.tool_calls);
    EXPECT_DOUBLE_EQ(0.0, result.cost_usd);
    EXPECT_FALSE(result.output_path.has_value());
    EXPECT_FALSE(result.error_message.has_value());
}

TEST_F(AgentTypesTest, BucketAnnotationDefaults) {
    BucketAnnotation ann;
    EXPECT_EQ(0u, ann.bucket_id);
    EXPECT_TRUE(ann.level.empty());
    EXPECT_TRUE(ann.annotation.empty());
    EXPECT_DOUBLE_EQ(1.0, ann.prior_weight);
    EXPECT_FALSE(ann.paralog_partner_bucket.has_value());
    EXPECT_DOUBLE_EQ(1.0, ann.expected_coverage_multiplier);
    EXPECT_TRUE(ann.claude_rationale.empty());
}

TEST_F(AgentTypesTest, BiologyPriorEmpty) {
    BiologyPrior prior;
    EXPECT_TRUE(prior.empty());
    EXPECT_EQ("1.0", prior.version);
    EXPECT_TRUE(prior.reference_sha256.empty());
}

TEST_F(AgentTypesTest, BiologyPriorNonEmpty) {
    BiologyPrior prior;
    BucketAnnotation ann;
    ann.bucket_id = 12345;
    ann.level = "L2";
    ann.annotation = "Test region";
    prior.buckets[12345] = ann;
    EXPECT_FALSE(prior.empty());
}

TEST_F(AgentTypesTest, RegionalOverrideDefaults) {
    RegionalOverride ovr;
    EXPECT_TRUE(ovr.region.empty());
    EXPECT_EQ(50u, ovr.sub_bucket_granularity_kb);
    EXPECT_EQ(20u, ovr.max_iter);
    EXPECT_DOUBLE_EQ(0.99, ovr.convergence_threshold);
}

TEST_F(AgentTypesTest, SampleParamsDefaults) {
    SampleParams params;
    EXPECT_EQ("hifi", params.preset);
    EXPECT_DOUBLE_EQ(0.99, params.convergence_threshold);
    EXPECT_EQ(20u, params.max_iterations);
    EXPECT_TRUE(params.expected_coverage_profile.empty());
    EXPECT_EQ("caduceus-ph", params.foundation_model);
}

TEST_F(AgentTypesTest, DiagnosticReportDefaults) {
    DiagnosticReport report;
    EXPECT_TRUE(report.stall_pattern.empty());
    EXPECT_TRUE(report.root_cause.empty());
    EXPECT_TRUE(report.resolution.empty());
    EXPECT_FALSE(report.custom_kernel_path.has_value());
    EXPECT_FALSE(report.kernel_hot_loaded);
}

TEST_F(AgentTypesTest, AnalysisReportDefaults) {
    AnalysisReport report;
    EXPECT_TRUE(report.sample_id.empty());
    EXPECT_TRUE(report.summary.empty());
    EXPECT_TRUE(report.region_notes.empty());
    EXPECT_TRUE(report.warnings.empty());
    EXPECT_TRUE(report.recommendations.empty());
}

TEST_F(AgentTypesTest, SessionOutputVariant) {
    SessionOutput output1 = BiologyPrior{};
    EXPECT_TRUE(std::holds_alternative<BiologyPrior>(output1));

    SessionOutput output2 = SampleParams{};
    EXPECT_TRUE(std::holds_alternative<SampleParams>(output2));

    SessionOutput output3 = DiagnosticReport{};
    EXPECT_TRUE(std::holds_alternative<DiagnosticReport>(output3));

    SessionOutput output4 = AnalysisReport{};
    EXPECT_TRUE(std::holds_alternative<AnalysisReport>(output4));
}
