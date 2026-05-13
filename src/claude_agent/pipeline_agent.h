#pragma once

#include "agent_session.h"
#include "agent_types.h"
#include "cuda_sandbox.h"

#include <chrono>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace llmap::claude_agent {

enum class StallType : uint8_t {
    None = 0,
    NoProgress = 1,
    OscillatingProbabilities = 2,
    HighUncertainty = 3,
    ParalogConfusion = 4,
    CoverageAnomaly = 5
};

constexpr const char* StallTypeName(StallType t) {
    switch (t) {
        case StallType::None:                    return "none";
        case StallType::NoProgress:              return "no-progress";
        case StallType::OscillatingProbabilities:return "oscillating-probabilities";
        case StallType::HighUncertainty:         return "high-uncertainty";
        case StallType::ParalogConfusion:        return "paralog-confusion";
        case StallType::CoverageAnomaly:         return "coverage-anomaly";
    }
    return "unknown";
}

struct StallMetrics {
    StallType type{StallType::None};
    size_t iteration{0};
    size_t reads_affected{0};
    double entropy_delta{0.0};
    double max_oscillation{0.0};
    std::vector<std::string> affected_regions;
};

struct PipelineAgentConfig {
    AgentConfig agent;
    std::filesystem::path work_dir;
    bool enable_diagnostics{true};
    bool enable_reporter{true};
    bool enable_cuda_codegen{false};
    size_t stall_detection_window{5};
    double stall_entropy_threshold{0.01};
    double oscillation_threshold{0.1};
    std::chrono::seconds diagnostic_timeout{300};
};

struct DiagnosticContext {
    StallMetrics stall;
    std::filesystem::path wave_state_path;
    std::filesystem::path output_dir;
};

struct ReporterContext {
    std::filesystem::path alignment_path;
    std::filesystem::path wave_state_summary;
    std::string sample_id;
    std::filesystem::path output_dir;
};

struct ResolutionResult {
    bool resolved{false};
    std::optional<std::string> kernel_name;
    std::optional<LoadedKernel> loaded_kernel;
    DiagnosticReport report;
    std::chrono::milliseconds latency{0};
};

using IterationCallback = std::function<void(size_t iteration, double entropy)>;
using StallCallback = std::function<void(const StallMetrics& stall)>;

class StallDetector {
public:
    explicit StallDetector(size_t window_size = 5,
                           double entropy_threshold = 0.01,
                           double oscillation_threshold = 0.1);
    ~StallDetector();

    StallDetector(const StallDetector&) = delete;
    StallDetector& operator=(const StallDetector&) = delete;
    StallDetector(StallDetector&&) noexcept;
    StallDetector& operator=(StallDetector&&) noexcept;

    void RecordIteration(size_t iteration, double entropy,
                         size_t uncertain_reads = 0,
                         const std::vector<std::string>& regions = {});

    StallMetrics DetectStall() const;
    bool IsStalled() const;
    void Reset();

    size_t IterationCount() const;
    double CurrentEntropy() const;
    double EntropyDelta() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class PipelineAgent {
public:
    explicit PipelineAgent(PipelineAgentConfig config);
    ~PipelineAgent();

    PipelineAgent(const PipelineAgent&) = delete;
    PipelineAgent& operator=(const PipelineAgent&) = delete;
    PipelineAgent(PipelineAgent&&) noexcept;
    PipelineAgent& operator=(PipelineAgent&&) noexcept;

    void RecordIteration(size_t iteration, double entropy,
                         size_t uncertain_reads = 0,
                         const std::vector<std::string>& regions = {});

    bool IsStalled() const;
    StallMetrics GetStallMetrics() const;

    ResolutionResult DiagnoseAndResolve(const DiagnosticContext& ctx);

    AnalysisReport GenerateReport(const ReporterContext& ctx);

    std::optional<LoadedKernel> GetLoadedKernel(std::string_view name) const;

    void SetIterationCallback(IterationCallback cb);
    void SetStallCallback(StallCallback cb);

    const PipelineAgentConfig& GetConfig() const;
    std::string SessionId() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::filesystem::path WriteWaveStateJson(
    const std::filesystem::path& output_dir,
    size_t iteration,
    double entropy,
    const std::vector<std::pair<std::string, double>>& bucket_probabilities);

StallMetrics AnalyzeWaveStateForStall(const std::filesystem::path& wave_state_json);

}  // namespace llmap::claude_agent
