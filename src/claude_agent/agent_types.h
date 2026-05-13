#pragma once

#include <chrono>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace llmap::claude_agent {

enum class AgentMode : uint8_t {
    Off = 0,
    IndexOnly = 1,
    SampleAware = 2,
    SelfHealing = 3,
    Research = 4
};

enum class SessionType : uint8_t {
    IndexBuild = 0,  // Session A
    SampleInit = 1,  // Session B
    Diagnostic = 2,  // Session C
    Reporter = 3     // Session D
};

constexpr const char* SessionTypeName(SessionType t) {
    switch (t) {
        case SessionType::IndexBuild: return "index-build";
        case SessionType::SampleInit: return "sample-init";
        case SessionType::Diagnostic: return "diagnostic";
        case SessionType::Reporter:   return "reporter";
    }
    return "unknown";
}

constexpr const char* AgentModeName(AgentMode m) {
    switch (m) {
        case AgentMode::Off:         return "off";
        case AgentMode::IndexOnly:   return "index-only";
        case AgentMode::SampleAware: return "sample-aware";
        case AgentMode::SelfHealing: return "self-healing";
        case AgentMode::Research:    return "research";
    }
    return "unknown";
}

std::optional<AgentMode> ParseAgentMode(std::string_view s);

struct AgentConfig {
    AgentMode mode{AgentMode::SampleAware};
    std::string api_key;
    std::string model{"claude-sonnet-4-20250514"};
    size_t max_tokens{8192};
    double temperature{0.0};
    size_t rate_limit_rpm{50};
    size_t max_retries{3};
    std::chrono::milliseconds retry_delay{1000};
    std::optional<std::string> cache_dir;
};

enum class AgentStatus : uint8_t {
    Pending = 0,
    Running = 1,
    Completed = 2,
    Failed = 3,
    TimedOut = 4
};

struct AgentResult {
    AgentStatus status{AgentStatus::Pending};
    std::string session_id;
    SessionType session_type{};
    size_t tool_calls{0};
    double cost_usd{0.0};
    std::chrono::milliseconds latency{0};
    std::optional<std::string> output_path;
    std::optional<std::string> error_message;
};

struct BucketAnnotation {
    uint64_t bucket_id{0};
    std::string level;
    std::string annotation;
    double prior_weight{1.0};
    std::optional<uint64_t> paralog_partner_bucket;
    double expected_coverage_multiplier{1.0};
    std::string claude_rationale;
};

struct RegionalOverride {
    std::string region;
    size_t sub_bucket_granularity_kb{50};
    size_t max_iter{20};
    double convergence_threshold{0.99};
};

struct BiologyPrior {
    std::string version{"1.0"};
    std::string reference_sha256;
    std::map<uint64_t, BucketAnnotation> buckets;
    std::map<std::string, RegionalOverride> regional_overrides;

    bool empty() const { return buckets.empty() && regional_overrides.empty(); }
};

struct SampleParams {
    std::string preset{"hifi"};
    double convergence_threshold{0.99};
    size_t max_iterations{20};
    std::string expected_coverage_profile;
    std::string foundation_model{"caduceus-ph"};
    std::map<std::string, double> region_adjustments;
};

struct DiagnosticReport {
    std::string stall_pattern;
    std::string root_cause;
    std::string resolution;
    std::optional<std::string> custom_kernel_path;
    bool kernel_hot_loaded{false};
};

struct AnalysisReport {
    std::string sample_id;
    std::string summary;
    std::map<std::string, std::string> region_notes;
    std::vector<std::string> warnings;
    std::vector<std::string> recommendations;
};

using SessionOutput = std::variant<
    BiologyPrior,
    SampleParams,
    DiagnosticReport,
    AnalysisReport
>;

}  // namespace llmap::claude_agent
