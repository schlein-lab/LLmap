// LLmap — Pipeline agent integration for alignment stall detection and resolution.

#include "claude_agent/pipeline_agent.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>

namespace llmap::claude_agent {

namespace {

std::atomic<uint64_t> g_session_counter{0};

std::string GenerateSessionId() {
    auto now = std::chrono::system_clock::now();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
        now.time_since_epoch()).count();
    auto counter = g_session_counter.fetch_add(1, std::memory_order_relaxed);
    std::ostringstream oss;
    oss << "pipeline-" << us << "-" << counter;
    return oss.str();
}

}  // namespace

// --- StallDetector ---

struct StallDetector::Impl {
    size_t window_size{5};
    double entropy_threshold{0.01};
    double oscillation_threshold{0.1};

    struct IterationData {
        size_t iteration{0};
        double entropy{0.0};
        size_t uncertain_reads{0};
        std::vector<std::string> regions;
    };

    std::vector<IterationData> history;
};

StallDetector::StallDetector(size_t window_size,
                             double entropy_threshold,
                             double oscillation_threshold)
    : impl_(std::make_unique<Impl>()) {
    impl_->window_size = window_size;
    impl_->entropy_threshold = entropy_threshold;
    impl_->oscillation_threshold = oscillation_threshold;
}

StallDetector::~StallDetector() = default;
StallDetector::StallDetector(StallDetector&&) noexcept = default;
StallDetector& StallDetector::operator=(StallDetector&&) noexcept = default;

void StallDetector::RecordIteration(size_t iteration, double entropy,
                                    size_t uncertain_reads,
                                    const std::vector<std::string>& regions) {
    impl_->history.push_back({iteration, entropy, uncertain_reads, regions});

    if (impl_->history.size() > impl_->window_size * 2) {
        impl_->history.erase(impl_->history.begin());
    }
}

StallMetrics StallDetector::DetectStall() const {
    StallMetrics metrics;

    if (impl_->history.size() < impl_->window_size) {
        return metrics;
    }

    const auto& recent = impl_->history;
    size_t n = recent.size();

    double first_entropy = recent[n - impl_->window_size].entropy;
    double last_entropy = recent.back().entropy;
    double delta = std::abs(last_entropy - first_entropy);

    // Check oscillation first - oscillating patterns are a specific stall type
    double max_osc = 0.0;
    size_t direction_changes = 0;
    for (size_t i = n - impl_->window_size + 1; i < n; ++i) {
        double diff = recent[i].entropy - recent[i - 1].entropy;
        double prev_diff = (i > n - impl_->window_size + 1)
            ? recent[i - 1].entropy - recent[i - 2].entropy
            : 0.0;

        if (i > n - impl_->window_size + 1) {
            if ((diff > 0 && prev_diff < 0) || (diff < 0 && prev_diff > 0)) {
                direction_changes++;
                max_osc = std::max(max_osc, std::abs(diff));
            }
        }
    }

    // Oscillation requires multiple direction changes and significant magnitude
    if (direction_changes >= 2 && max_osc > impl_->oscillation_threshold) {
        metrics.type = StallType::OscillatingProbabilities;
        metrics.iteration = recent.back().iteration;
        metrics.max_oscillation = max_osc;
        metrics.entropy_delta = delta;
        metrics.reads_affected = recent.back().uncertain_reads;
        metrics.affected_regions = recent.back().regions;
        return metrics;
    }

    if (delta < impl_->entropy_threshold) {
        metrics.type = StallType::NoProgress;
        metrics.iteration = recent.back().iteration;
        metrics.entropy_delta = delta;
        metrics.reads_affected = recent.back().uncertain_reads;
        metrics.affected_regions = recent.back().regions;
        return metrics;
    }

    if (recent.back().uncertain_reads > 0) {
        double uncertainty_ratio = static_cast<double>(recent.back().uncertain_reads) /
            std::max(1.0, static_cast<double>(recent[0].uncertain_reads));
        if (uncertainty_ratio > 0.8 && impl_->history.size() >= impl_->window_size) {
            metrics.type = StallType::HighUncertainty;
            metrics.iteration = recent.back().iteration;
            metrics.reads_affected = recent.back().uncertain_reads;
            metrics.affected_regions = recent.back().regions;
        }
    }

    return metrics;
}

bool StallDetector::IsStalled() const {
    return DetectStall().type != StallType::None;
}

void StallDetector::Reset() {
    impl_->history.clear();
}

size_t StallDetector::IterationCount() const {
    return impl_->history.size();
}

double StallDetector::CurrentEntropy() const {
    if (impl_->history.empty()) return 0.0;
    return impl_->history.back().entropy;
}

double StallDetector::EntropyDelta() const {
    if (impl_->history.size() < 2) return 0.0;
    return impl_->history.back().entropy - impl_->history[impl_->history.size() - 2].entropy;
}

// --- PipelineAgent ---

struct PipelineAgent::Impl {
    PipelineAgentConfig config;
    std::string session_id;
    StallDetector detector;
    std::unique_ptr<CudaSandbox> sandbox;
    IterationCallback iteration_cb;
    StallCallback stall_cb;

    explicit Impl(PipelineAgentConfig cfg)
        : config(std::move(cfg))
        , session_id(GenerateSessionId())
        , detector(cfg.stall_detection_window,
                   cfg.stall_entropy_threshold,
                   cfg.oscillation_threshold) {
        if (config.enable_cuda_codegen) {
            CudaSandbox::Config sandbox_cfg;
            sandbox_cfg.session_id = session_id;
            sandbox_cfg.compiler.sandbox_dir = config.work_dir / "sandbox";
            sandbox_cfg.compiler.output_dir = config.work_dir / "kernels";
            sandbox_cfg.audit_log_path = config.work_dir / "audit.log";
            sandbox = std::make_unique<CudaSandbox>(std::move(sandbox_cfg));
        }
    }
};

PipelineAgent::PipelineAgent(PipelineAgentConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {
}

PipelineAgent::~PipelineAgent() = default;
PipelineAgent::PipelineAgent(PipelineAgent&&) noexcept = default;
PipelineAgent& PipelineAgent::operator=(PipelineAgent&&) noexcept = default;

void PipelineAgent::RecordIteration(size_t iteration, double entropy,
                                    size_t uncertain_reads,
                                    const std::vector<std::string>& regions) {
    impl_->detector.RecordIteration(iteration, entropy, uncertain_reads, regions);

    if (impl_->iteration_cb) {
        impl_->iteration_cb(iteration, entropy);
    }

    if (impl_->detector.IsStalled() && impl_->stall_cb) {
        impl_->stall_cb(impl_->detector.DetectStall());
    }
}

bool PipelineAgent::IsStalled() const {
    return impl_->detector.IsStalled();
}

StallMetrics PipelineAgent::GetStallMetrics() const {
    return impl_->detector.DetectStall();
}

ResolutionResult PipelineAgent::DiagnoseAndResolve(const DiagnosticContext& ctx) {
    ResolutionResult result;
    auto start = std::chrono::steady_clock::now();

    if (!impl_->config.enable_diagnostics) {
        result.report.stall_pattern = StallTypeName(ctx.stall.type);
        result.report.root_cause = "diagnostics-disabled";
        result.report.resolution = "none";
        return result;
    }

    DiagnosticSession::Config diag_cfg;
    diag_cfg.wave_state_json = ctx.wave_state_path;
    diag_cfg.output_dir = ctx.output_dir;
    diag_cfg.enable_cuda_codegen = impl_->config.enable_cuda_codegen;
    diag_cfg.agent_config = impl_->config.agent;

    DiagnosticSession session(std::move(diag_cfg));
    auto agent_result = session.Run();

    if (agent_result.status != AgentStatus::Completed) {
        result.report.stall_pattern = StallTypeName(ctx.stall.type);
        result.report.root_cause = "agent-failed";
        result.report.resolution = agent_result.error_message.value_or("unknown error");
        return result;
    }

    result.report = session.GetDiagnosticReport();

    if (result.report.custom_kernel_path && impl_->sandbox) {
        std::ifstream kernel_file(*result.report.custom_kernel_path);
        if (kernel_file) {
            std::ostringstream oss;
            oss << kernel_file.rdbuf();
            std::string kernel_source = oss.str();

            auto exec_result = impl_->sandbox->CompileAndLoad(
                kernel_source, "diagnostic_kernel");

            if (exec_result.success) {
                auto kernel = impl_->sandbox->GetKernel("diagnostic_kernel");
                if (kernel) {
                    result.resolved = true;
                    result.kernel_name = "diagnostic_kernel";
                    result.loaded_kernel = *kernel;
                    result.report.kernel_hot_loaded = true;
                }
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    result.latency = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    return result;
}

AnalysisReport PipelineAgent::GenerateReport(const ReporterContext& ctx) {
    AnalysisReport report;
    report.sample_id = ctx.sample_id;

    if (!impl_->config.enable_reporter) {
        report.summary = "reporter-disabled";
        return report;
    }

    ReporterSession::Config rep_cfg;
    rep_cfg.alignment_bam = ctx.alignment_path;
    rep_cfg.wave_state_summary = ctx.wave_state_summary;
    rep_cfg.output_dir = ctx.output_dir;
    rep_cfg.sample_id = ctx.sample_id;
    rep_cfg.agent_config = impl_->config.agent;

    ReporterSession session(std::move(rep_cfg));
    auto agent_result = session.Run();

    if (agent_result.status == AgentStatus::Completed) {
        report = session.GetAnalysisReport();
    } else {
        report.summary = "agent-failed";
        report.warnings.push_back(
            agent_result.error_message.value_or("unknown error"));
    }

    return report;
}

std::optional<LoadedKernel> PipelineAgent::GetLoadedKernel(
    std::string_view name) const {
    if (!impl_->sandbox) return std::nullopt;
    return impl_->sandbox->GetKernel(name);
}

void PipelineAgent::SetIterationCallback(IterationCallback cb) {
    impl_->iteration_cb = std::move(cb);
}

void PipelineAgent::SetStallCallback(StallCallback cb) {
    impl_->stall_cb = std::move(cb);
}

const PipelineAgentConfig& PipelineAgent::GetConfig() const {
    return impl_->config;
}

std::string PipelineAgent::SessionId() const {
    return impl_->session_id;
}

// --- Utility functions ---

std::filesystem::path WriteWaveStateJson(
    const std::filesystem::path& output_dir,
    size_t iteration,
    double entropy,
    const std::vector<std::pair<std::string, double>>& bucket_probabilities) {

    std::filesystem::create_directories(output_dir);

    std::ostringstream filename;
    filename << "wave_state_iter_" << iteration << ".json";
    auto path = output_dir / filename.str();

    std::ofstream out(path);
    if (!out) return {};

    out << "{\n";
    out << "  \"iteration\": " << iteration << ",\n";
    out << "  \"entropy\": " << entropy << ",\n";
    out << "  \"bucket_probabilities\": [\n";

    for (size_t i = 0; i < bucket_probabilities.size(); ++i) {
        const auto& [bucket_id, prob] = bucket_probabilities[i];
        out << "    {\"bucket_id\": \"" << bucket_id << "\", \"probability\": " << prob << "}";
        if (i + 1 < bucket_probabilities.size()) out << ",";
        out << "\n";
    }

    out << "  ]\n";
    out << "}\n";

    return path;
}

StallMetrics AnalyzeWaveStateForStall(const std::filesystem::path& wave_state_json) {
    StallMetrics metrics;

    std::ifstream in(wave_state_json);
    if (!in) return metrics;

    std::string content((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());

    auto find_value = [&](const std::string& key) -> std::optional<double> {
        auto pos = content.find("\"" + key + "\"");
        if (pos == std::string::npos) return std::nullopt;
        pos = content.find(":", pos);
        if (pos == std::string::npos) return std::nullopt;
        pos++;
        while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\t')) {
            pos++;
        }
        size_t end = pos;
        while (end < content.size() &&
               (std::isdigit(content[end]) || content[end] == '.' || content[end] == '-')) {
            end++;
        }
        if (end == pos) return std::nullopt;
        return std::stod(content.substr(pos, end - pos));
    };

    if (auto iter = find_value("iteration")) {
        metrics.iteration = static_cast<size_t>(*iter);
    }
    if (auto entropy = find_value("entropy")) {
        metrics.entropy_delta = *entropy;
    }

    return metrics;
}

}  // namespace llmap::claude_agent
