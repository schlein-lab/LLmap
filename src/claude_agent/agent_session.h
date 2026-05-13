#pragma once

#include "agent_types.h"
#include "anthropic_client.h"

#include <filesystem>
#include <future>
#include <memory>
#include <string>

namespace llmap::claude_agent {

class AgentSession {
public:
    virtual ~AgentSession() = default;

    virtual SessionType Type() const = 0;
    virtual std::future<AgentResult> RunAsync() = 0;
    virtual AgentResult Run() = 0;

    std::string SessionId() const { return session_id_; }

protected:
    explicit AgentSession(std::string session_id)
        : session_id_(std::move(session_id)) {
    }

    std::string session_id_;
};

class IndexBuildSession : public AgentSession {
public:
    struct Config {
        std::filesystem::path reference_fasta;
        std::filesystem::path output_dir;
        std::optional<std::filesystem::path> gff_file;
        std::optional<std::filesystem::path> repeatmasker_bed;
        AgentConfig agent_config;
    };

    explicit IndexBuildSession(Config config);
    ~IndexBuildSession() override;

    SessionType Type() const override { return SessionType::IndexBuild; }
    std::future<AgentResult> RunAsync() override;
    AgentResult Run() override;

    BiologyPrior GetBiologyPrior() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class SampleInitSession : public AgentSession {
public:
    struct Config {
        std::filesystem::path fastq_path;
        std::filesystem::path output_dir;
        std::optional<std::filesystem::path> sample_sheet;
        AgentConfig agent_config;
    };

    explicit SampleInitSession(Config config);
    ~SampleInitSession() override;

    SessionType Type() const override { return SessionType::SampleInit; }
    std::future<AgentResult> RunAsync() override;
    AgentResult Run() override;

    SampleParams GetSampleParams() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class DiagnosticSession : public AgentSession {
public:
    struct Config {
        std::filesystem::path wave_state_json;
        std::filesystem::path output_dir;
        bool enable_cuda_codegen{false};
        std::filesystem::path sandbox_script;
        AgentConfig agent_config;
    };

    explicit DiagnosticSession(Config config);
    ~DiagnosticSession() override;

    SessionType Type() const override { return SessionType::Diagnostic; }
    std::future<AgentResult> RunAsync() override;
    AgentResult Run() override;

    DiagnosticReport GetDiagnosticReport() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class ReporterSession : public AgentSession {
public:
    struct Config {
        std::filesystem::path alignment_bam;
        std::filesystem::path wave_state_summary;
        std::filesystem::path output_dir;
        std::string sample_id;
        AgentConfig agent_config;
    };

    explicit ReporterSession(Config config);
    ~ReporterSession() override;

    SessionType Type() const override { return SessionType::Reporter; }
    std::future<AgentResult> RunAsync() override;
    AgentResult Run() override;

    AnalysisReport GetAnalysisReport() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::unique_ptr<AgentSession> CreateSession(
    SessionType type,
    AgentConfig agent_config,
    const std::filesystem::path& work_dir
);

}  // namespace llmap::claude_agent
