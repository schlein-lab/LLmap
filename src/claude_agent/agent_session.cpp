#include "agent_session.h"

#include <chrono>
#include <sstream>
#include <random>

namespace llmap::claude_agent {

namespace {

std::string GenerateSessionId(SessionType type) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, 0xFFFFFFFF);

    std::ostringstream oss;
    oss << SessionTypeName(type) << "-"
        << std::hex << dist(gen) << dist(gen);
    return oss.str();
}

}  // namespace

// IndexBuildSession implementation

struct IndexBuildSession::Impl {
    Config config;
    AnthropicClient client;
    BiologyPrior prior;
    AgentResult result;

    explicit Impl(Config cfg)
        : config(std::move(cfg))
        , client(config.agent_config) {
    }
};

IndexBuildSession::IndexBuildSession(Config config)
    : AgentSession(GenerateSessionId(SessionType::IndexBuild))
    , impl_(std::make_unique<Impl>(std::move(config))) {
}

IndexBuildSession::~IndexBuildSession() = default;

std::future<AgentResult> IndexBuildSession::RunAsync() {
    return std::async(std::launch::async, [this]() {
        return Run();
    });
}

AgentResult IndexBuildSession::Run() {
    auto start = std::chrono::steady_clock::now();

    impl_->result.session_id = session_id_;
    impl_->result.session_type = SessionType::IndexBuild;
    impl_->result.status = AgentStatus::Running;

    std::ostringstream system_prompt;
    system_prompt << "You are an index-build agent for LLmap.\n"
                  << "Your goal is to analyze the reference genome and produce "
                  << "biology_prior.json with bucket-level annotations.\n"
                  << "Reference: " << impl_->config.reference_fasta << "\n";

    if (impl_->config.gff_file) {
        system_prompt << "GFF: " << *impl_->config.gff_file << "\n";
    }
    if (impl_->config.repeatmasker_bed) {
        system_prompt << "RepeatMasker: " << *impl_->config.repeatmasker_bed << "\n";
    }

    std::string user_msg = "Analyze the reference and identify regions that "
                           "require special handling for paralog disambiguation.";

    if (impl_->client.HasApiKey()) {
        auto turn = impl_->client.RunConversation(
            system_prompt.str(), user_msg, 50);

        impl_->result.tool_calls = impl_->client.TotalToolCalls();
        impl_->result.cost_usd = impl_->client.EstimateCostUsd();
        impl_->result.status = AgentStatus::Completed;
    } else {
        // No API key - generate empty prior
        impl_->prior.version = "1.0";
        impl_->result.status = AgentStatus::Completed;
    }

    auto end = std::chrono::steady_clock::now();
    impl_->result.latency = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start);

    auto output_path = impl_->config.output_dir / "biology_prior.json";
    impl_->result.output_path = output_path.string();

    return impl_->result;
}

BiologyPrior IndexBuildSession::GetBiologyPrior() const {
    return impl_->prior;
}

// SampleInitSession implementation

struct SampleInitSession::Impl {
    Config config;
    AnthropicClient client;
    SampleParams params;
    AgentResult result;

    explicit Impl(Config cfg)
        : config(std::move(cfg))
        , client(config.agent_config) {
    }
};

SampleInitSession::SampleInitSession(Config config)
    : AgentSession(GenerateSessionId(SessionType::SampleInit))
    , impl_(std::make_unique<Impl>(std::move(config))) {
}

SampleInitSession::~SampleInitSession() = default;

std::future<AgentResult> SampleInitSession::RunAsync() {
    return std::async(std::launch::async, [this]() {
        return Run();
    });
}

AgentResult SampleInitSession::Run() {
    auto start = std::chrono::steady_clock::now();

    impl_->result.session_id = session_id_;
    impl_->result.session_type = SessionType::SampleInit;
    impl_->result.status = AgentStatus::Running;

    std::ostringstream system_prompt;
    system_prompt << "You are a sample-init agent for LLmap.\n"
                  << "Your goal is to analyze the input FASTQ and choose "
                  << "optimal parameters for the alignment run.\n"
                  << "FASTQ: " << impl_->config.fastq_path << "\n";

    if (impl_->config.sample_sheet) {
        system_prompt << "Sample sheet: " << *impl_->config.sample_sheet << "\n";
    }

    std::string user_msg = "Analyze the sample and determine optimal parameters.";

    if (impl_->client.HasApiKey()) {
        auto turn = impl_->client.RunConversation(
            system_prompt.str(), user_msg, 20);

        impl_->result.tool_calls = impl_->client.TotalToolCalls();
        impl_->result.cost_usd = impl_->client.EstimateCostUsd();
        impl_->result.status = AgentStatus::Completed;
    } else {
        impl_->params.preset = "hifi";
        impl_->params.convergence_threshold = 0.99;
        impl_->params.max_iterations = 20;
        impl_->result.status = AgentStatus::Completed;
    }

    auto end = std::chrono::steady_clock::now();
    impl_->result.latency = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start);

    auto output_path = impl_->config.output_dir / "sample_params.json";
    impl_->result.output_path = output_path.string();

    return impl_->result;
}

SampleParams SampleInitSession::GetSampleParams() const {
    return impl_->params;
}

// DiagnosticSession implementation

struct DiagnosticSession::Impl {
    Config config;
    AnthropicClient client;
    DiagnosticReport report;
    AgentResult result;

    explicit Impl(Config cfg)
        : config(std::move(cfg))
        , client(config.agent_config) {
    }
};

DiagnosticSession::DiagnosticSession(Config config)
    : AgentSession(GenerateSessionId(SessionType::Diagnostic))
    , impl_(std::make_unique<Impl>(std::move(config))) {
}

DiagnosticSession::~DiagnosticSession() = default;

std::future<AgentResult> DiagnosticSession::RunAsync() {
    return std::async(std::launch::async, [this]() {
        return Run();
    });
}

AgentResult DiagnosticSession::Run() {
    auto start = std::chrono::steady_clock::now();

    impl_->result.session_id = session_id_;
    impl_->result.session_type = SessionType::Diagnostic;
    impl_->result.status = AgentStatus::Running;

    std::ostringstream system_prompt;
    system_prompt << "You are a diagnostic agent for LLmap.\n"
                  << "The EM iteration has stalled. Analyze the wave state "
                  << "and diagnose the issue.\n"
                  << "Wave state: " << impl_->config.wave_state_json << "\n";

    if (impl_->config.enable_cuda_codegen) {
        system_prompt << "CUDA codegen is ENABLED. You may write custom kernels.\n"
                      << "Sandbox script: " << impl_->config.sandbox_script << "\n";
    }

    std::string user_msg = "Diagnose why EM is not converging and propose a fix.";

    if (impl_->client.HasApiKey()) {
        auto turn = impl_->client.RunConversation(
            system_prompt.str(), user_msg, 100);

        impl_->result.tool_calls = impl_->client.TotalToolCalls();
        impl_->result.cost_usd = impl_->client.EstimateCostUsd();
        impl_->result.status = AgentStatus::Completed;
    } else {
        impl_->report.stall_pattern = "unknown";
        impl_->report.root_cause = "API key not configured";
        impl_->report.resolution = "none";
        impl_->result.status = AgentStatus::Completed;
    }

    auto end = std::chrono::steady_clock::now();
    impl_->result.latency = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start);

    return impl_->result;
}

DiagnosticReport DiagnosticSession::GetDiagnosticReport() const {
    return impl_->report;
}

// ReporterSession implementation

struct ReporterSession::Impl {
    Config config;
    AnthropicClient client;
    AnalysisReport report;
    AgentResult result;

    explicit Impl(Config cfg)
        : config(std::move(cfg))
        , client(config.agent_config) {
    }
};

ReporterSession::ReporterSession(Config config)
    : AgentSession(GenerateSessionId(SessionType::Reporter))
    , impl_(std::make_unique<Impl>(std::move(config))) {
}

ReporterSession::~ReporterSession() = default;

std::future<AgentResult> ReporterSession::RunAsync() {
    return std::async(std::launch::async, [this]() {
        return Run();
    });
}

AgentResult ReporterSession::Run() {
    auto start = std::chrono::steady_clock::now();

    impl_->result.session_id = session_id_;
    impl_->result.session_type = SessionType::Reporter;
    impl_->result.status = AgentStatus::Running;

    impl_->report.sample_id = impl_->config.sample_id;

    std::ostringstream system_prompt;
    system_prompt << "You are a reporter agent for LLmap.\n"
                  << "Your goal is to generate a diagnostic report for the "
                  << "alignment run.\n"
                  << "BAM: " << impl_->config.alignment_bam << "\n"
                  << "Wave state summary: " << impl_->config.wave_state_summary << "\n"
                  << "Sample ID: " << impl_->config.sample_id << "\n";

    std::string user_msg = "Generate a comprehensive analysis report.";

    if (impl_->client.HasApiKey()) {
        auto turn = impl_->client.RunConversation(
            system_prompt.str(), user_msg, 30);

        impl_->result.tool_calls = impl_->client.TotalToolCalls();
        impl_->result.cost_usd = impl_->client.EstimateCostUsd();
        impl_->result.status = AgentStatus::Completed;
    } else {
        impl_->report.summary = "Analysis completed (no agent API key)";
        impl_->result.status = AgentStatus::Completed;
    }

    auto end = std::chrono::steady_clock::now();
    impl_->result.latency = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start);

    auto output_path = impl_->config.output_dir / "report.md";
    impl_->result.output_path = output_path.string();

    return impl_->result;
}

AnalysisReport ReporterSession::GetAnalysisReport() const {
    return impl_->report;
}

}  // namespace llmap::claude_agent
