#include <gtest/gtest.h>

#include "claude_agent/agent_session.h"

#include <filesystem>

using namespace llmap::claude_agent;

class AgentSessionTest : public ::testing::Test {
protected:
    std::filesystem::path temp_dir_;

    void SetUp() override {
        temp_dir_ = std::filesystem::temp_directory_path() /
                    ("agent_session_test_" + std::to_string(
                        std::hash<std::thread::id>{}(std::this_thread::get_id())));
        std::filesystem::create_directories(temp_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(temp_dir_);
    }
};

TEST_F(AgentSessionTest, IndexBuildSessionType) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));

    EXPECT_EQ(SessionType::IndexBuild, session.Type());
    EXPECT_FALSE(session.SessionId().empty());
    EXPECT_TRUE(session.SessionId().find("index-build") != std::string::npos);
}

TEST_F(AgentSessionTest, IndexBuildSessionRunWithoutApiKey) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
    EXPECT_EQ(SessionType::IndexBuild, result.session_type);
    EXPECT_EQ(0u, result.tool_calls);
    EXPECT_TRUE(result.output_path.has_value());
}

TEST_F(AgentSessionTest, IndexBuildSessionGetBiologyPrior) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));
    session.Run();

    auto prior = session.GetBiologyPrior();
    EXPECT_EQ("1.0", prior.version);
}

TEST_F(AgentSessionTest, SampleInitSessionType) {
    SampleInitSession::Config config;
    config.fastq_path = "/path/to/reads.fq";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    SampleInitSession session(std::move(config));

    EXPECT_EQ(SessionType::SampleInit, session.Type());
    EXPECT_TRUE(session.SessionId().find("sample-init") != std::string::npos);
}

TEST_F(AgentSessionTest, SampleInitSessionRunWithoutApiKey) {
    SampleInitSession::Config config;
    config.fastq_path = "/path/to/reads.fq";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    SampleInitSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
    EXPECT_EQ(SessionType::SampleInit, result.session_type);
}

TEST_F(AgentSessionTest, SampleInitSessionGetParams) {
    SampleInitSession::Config config;
    config.fastq_path = "/path/to/reads.fq";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    SampleInitSession session(std::move(config));
    session.Run();

    auto params = session.GetSampleParams();
    EXPECT_EQ("hifi", params.preset);
    EXPECT_DOUBLE_EQ(0.99, params.convergence_threshold);
}

TEST_F(AgentSessionTest, DiagnosticSessionType) {
    DiagnosticSession::Config config;
    config.wave_state_json = "/path/to/stalled.json";
    config.output_dir = temp_dir_;
    config.enable_cuda_codegen = false;
    config.agent_config.api_key = "";

    DiagnosticSession session(std::move(config));

    EXPECT_EQ(SessionType::Diagnostic, session.Type());
    EXPECT_TRUE(session.SessionId().find("diagnostic") != std::string::npos);
}

TEST_F(AgentSessionTest, DiagnosticSessionRunWithoutApiKey) {
    DiagnosticSession::Config config;
    config.wave_state_json = "/path/to/stalled.json";
    config.output_dir = temp_dir_;
    config.enable_cuda_codegen = false;
    config.agent_config.api_key = "";

    DiagnosticSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
    EXPECT_EQ(SessionType::Diagnostic, result.session_type);
}

TEST_F(AgentSessionTest, DiagnosticSessionGetReport) {
    DiagnosticSession::Config config;
    config.wave_state_json = "/path/to/stalled.json";
    config.output_dir = temp_dir_;
    config.enable_cuda_codegen = false;
    config.agent_config.api_key = "";

    DiagnosticSession session(std::move(config));
    session.Run();

    auto report = session.GetDiagnosticReport();
    EXPECT_EQ("API key not configured", report.root_cause);
}

TEST_F(AgentSessionTest, ReporterSessionType) {
    ReporterSession::Config config;
    config.alignment_bam = "/path/to/aligned.bam";
    config.wave_state_summary = "/path/to/summary.json";
    config.output_dir = temp_dir_;
    config.sample_id = "sample001";
    config.agent_config.api_key = "";

    ReporterSession session(std::move(config));

    EXPECT_EQ(SessionType::Reporter, session.Type());
    EXPECT_TRUE(session.SessionId().find("reporter") != std::string::npos);
}

TEST_F(AgentSessionTest, ReporterSessionRunWithoutApiKey) {
    ReporterSession::Config config;
    config.alignment_bam = "/path/to/aligned.bam";
    config.wave_state_summary = "/path/to/summary.json";
    config.output_dir = temp_dir_;
    config.sample_id = "sample001";
    config.agent_config.api_key = "";

    ReporterSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
    EXPECT_EQ(SessionType::Reporter, result.session_type);
    EXPECT_TRUE(result.output_path.has_value());
}

TEST_F(AgentSessionTest, ReporterSessionGetReport) {
    ReporterSession::Config config;
    config.alignment_bam = "/path/to/aligned.bam";
    config.wave_state_summary = "/path/to/summary.json";
    config.output_dir = temp_dir_;
    config.sample_id = "sample001";
    config.agent_config.api_key = "";

    ReporterSession session(std::move(config));
    session.Run();

    auto report = session.GetAnalysisReport();
    EXPECT_EQ("sample001", report.sample_id);
    EXPECT_FALSE(report.summary.empty());
}

TEST_F(AgentSessionTest, AsyncIndexBuildSession) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));
    auto future = session.RunAsync();

    auto result = future.get();
    EXPECT_EQ(AgentStatus::Completed, result.status);
}

TEST_F(AgentSessionTest, AsyncSampleInitSession) {
    SampleInitSession::Config config;
    config.fastq_path = "/path/to/reads.fq";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    SampleInitSession session(std::move(config));
    auto future = session.RunAsync();

    auto result = future.get();
    EXPECT_EQ(AgentStatus::Completed, result.status);
}

TEST_F(AgentSessionTest, SessionIdUnique) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    IndexBuildSession session1(config);
    IndexBuildSession session2(config);

    EXPECT_NE(session1.SessionId(), session2.SessionId());
}

TEST_F(AgentSessionTest, LatencyTracking) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));
    auto result = session.Run();

    // Latency is captured, may be 0ms for very fast operations
    EXPECT_GE(result.latency.count(), 0);
}

TEST_F(AgentSessionTest, DiagnosticWithCudaCodegen) {
    DiagnosticSession::Config config;
    config.wave_state_json = "/path/to/stalled.json";
    config.output_dir = temp_dir_;
    config.enable_cuda_codegen = true;
    config.sandbox_script = "/scripts/sandbox_compile.sh";
    config.agent_config.api_key = "";

    DiagnosticSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
}

TEST_F(AgentSessionTest, OptionalGffFile) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.gff_file = "/path/to/annotations.gff";
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
}

TEST_F(AgentSessionTest, OptionalRepeatMasker) {
    IndexBuildSession::Config config;
    config.reference_fasta = "/path/to/ref.fa";
    config.output_dir = temp_dir_;
    config.repeatmasker_bed = "/path/to/repeats.bed";
    config.agent_config.api_key = "";

    IndexBuildSession session(std::move(config));
    auto result = session.Run();

    EXPECT_EQ(AgentStatus::Completed, result.status);
}
