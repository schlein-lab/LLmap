// LLmap — Tests for pipeline agent integration.

#include <gtest/gtest.h>

#include "claude_agent/pipeline_agent.h"

#include <chrono>
#include <filesystem>
#include <fstream>

namespace llmap::claude_agent {
namespace {

class StallDetectorTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_stall_test";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path test_dir_;
};

TEST_F(StallDetectorTest, InitialStateNotStalled) {
    StallDetector detector;
    EXPECT_FALSE(detector.IsStalled());
    EXPECT_EQ(detector.IterationCount(), 0);
    EXPECT_EQ(detector.CurrentEntropy(), 0.0);
}

TEST_F(StallDetectorTest, RecordsIterations) {
    StallDetector detector;

    detector.RecordIteration(0, 2.5);
    EXPECT_EQ(detector.IterationCount(), 1);
    EXPECT_DOUBLE_EQ(detector.CurrentEntropy(), 2.5);

    detector.RecordIteration(1, 2.3);
    EXPECT_EQ(detector.IterationCount(), 2);
    EXPECT_DOUBLE_EQ(detector.CurrentEntropy(), 2.3);
}

TEST_F(StallDetectorTest, DetectsNoProgressStall) {
    StallDetector detector(5, 0.01, 0.1);

    for (int i = 0; i < 6; ++i) {
        detector.RecordIteration(i, 2.5);
    }

    EXPECT_TRUE(detector.IsStalled());
    auto metrics = detector.DetectStall();
    EXPECT_EQ(metrics.type, StallType::NoProgress);
    EXPECT_EQ(metrics.iteration, 5);
    EXPECT_LT(metrics.entropy_delta, 0.01);
}

TEST_F(StallDetectorTest, NoStallWithProgress) {
    StallDetector detector(5, 0.01, 0.1);

    for (int i = 0; i < 6; ++i) {
        detector.RecordIteration(i, 2.5 - i * 0.1);
    }

    EXPECT_FALSE(detector.IsStalled());
    auto metrics = detector.DetectStall();
    EXPECT_EQ(metrics.type, StallType::None);
}

TEST_F(StallDetectorTest, DetectsOscillation) {
    StallDetector detector(5, 0.001, 0.1);

    detector.RecordIteration(0, 2.0);
    detector.RecordIteration(1, 2.3);
    detector.RecordIteration(2, 2.0);
    detector.RecordIteration(3, 2.3);
    detector.RecordIteration(4, 2.0);
    detector.RecordIteration(5, 2.3);

    EXPECT_TRUE(detector.IsStalled());
    auto metrics = detector.DetectStall();
    EXPECT_EQ(metrics.type, StallType::OscillatingProbabilities);
    EXPECT_GT(metrics.max_oscillation, 0.1);
}

TEST_F(StallDetectorTest, ResetClearsHistory) {
    StallDetector detector;

    for (int i = 0; i < 10; ++i) {
        detector.RecordIteration(i, 2.5);
    }

    EXPECT_GT(detector.IterationCount(), 0);
    detector.Reset();
    EXPECT_EQ(detector.IterationCount(), 0);
    EXPECT_FALSE(detector.IsStalled());
}

TEST_F(StallDetectorTest, EntropyDeltaCalculation) {
    StallDetector detector;

    detector.RecordIteration(0, 2.5);
    EXPECT_NEAR(detector.EntropyDelta(), 0.0, 1e-10);

    detector.RecordIteration(1, 2.3);
    EXPECT_NEAR(detector.EntropyDelta(), -0.2, 1e-10);

    detector.RecordIteration(2, 2.4);
    EXPECT_NEAR(detector.EntropyDelta(), 0.1, 1e-10);
}

TEST_F(StallDetectorTest, TracksAffectedRegions) {
    StallDetector detector(3, 0.01, 0.1);

    std::vector<std::string> regions1 = {"chr1:1000-2000"};
    std::vector<std::string> regions2 = {"chr1:1000-2000", "chr2:5000-6000"};

    detector.RecordIteration(0, 2.5, 100, regions1);
    detector.RecordIteration(1, 2.5, 150, regions2);
    detector.RecordIteration(2, 2.5, 200, regions2);

    auto metrics = detector.DetectStall();
    EXPECT_EQ(metrics.type, StallType::NoProgress);
    EXPECT_EQ(metrics.reads_affected, 200);
    EXPECT_EQ(metrics.affected_regions.size(), 2);
}

// --- PipelineAgent tests ---

class PipelineAgentTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_pipeline_agent_test";
        std::filesystem::create_directories(test_dir_);

        config_.work_dir = test_dir_;
        config_.enable_diagnostics = true;
        config_.enable_reporter = true;
        config_.enable_cuda_codegen = false;
        config_.stall_detection_window = 3;
        config_.stall_entropy_threshold = 0.01;
        config_.oscillation_threshold = 0.1;
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path test_dir_;
    PipelineAgentConfig config_;
};

TEST_F(PipelineAgentTest, Construction) {
    PipelineAgent agent(config_);

    EXPECT_FALSE(agent.IsStalled());
    EXPECT_FALSE(agent.SessionId().empty());
    EXPECT_TRUE(agent.SessionId().starts_with("pipeline-"));
}

TEST_F(PipelineAgentTest, RecordsIterationsAndDetectsStall) {
    PipelineAgent agent(config_);

    for (int i = 0; i < 4; ++i) {
        agent.RecordIteration(i, 2.5);
    }

    EXPECT_TRUE(agent.IsStalled());
    auto metrics = agent.GetStallMetrics();
    EXPECT_EQ(metrics.type, StallType::NoProgress);
}

TEST_F(PipelineAgentTest, IterationCallbackInvoked) {
    PipelineAgent agent(config_);

    size_t callback_count = 0;
    double last_entropy = 0.0;

    agent.SetIterationCallback([&](size_t iter, double entropy) {
        callback_count++;
        last_entropy = entropy;
    });

    agent.RecordIteration(0, 2.5);
    agent.RecordIteration(1, 2.3);

    EXPECT_EQ(callback_count, 2);
    EXPECT_DOUBLE_EQ(last_entropy, 2.3);
}

TEST_F(PipelineAgentTest, StallCallbackInvoked) {
    PipelineAgent agent(config_);

    size_t stall_callback_count = 0;
    StallType detected_type = StallType::None;

    agent.SetStallCallback([&](const StallMetrics& stall) {
        stall_callback_count++;
        detected_type = stall.type;
    });

    for (int i = 0; i < 4; ++i) {
        agent.RecordIteration(i, 2.5);
    }

    EXPECT_GE(stall_callback_count, 1);
    EXPECT_EQ(detected_type, StallType::NoProgress);
}

TEST_F(PipelineAgentTest, ConfigAccessible) {
    PipelineAgent agent(config_);

    const auto& cfg = agent.GetConfig();
    EXPECT_EQ(cfg.work_dir, test_dir_);
    EXPECT_TRUE(cfg.enable_diagnostics);
    EXPECT_EQ(cfg.stall_detection_window, 3);
}

TEST_F(PipelineAgentTest, MoveConstruction) {
    PipelineAgent agent1(config_);
    agent1.RecordIteration(0, 2.5);

    PipelineAgent agent2 = std::move(agent1);
    EXPECT_FALSE(agent2.SessionId().empty());
}

// --- DiagnosticContext tests ---

TEST_F(PipelineAgentTest, DiagnoseWithDisabledDiagnostics) {
    config_.enable_diagnostics = false;
    PipelineAgent agent(config_);

    DiagnosticContext ctx;
    ctx.stall.type = StallType::NoProgress;
    ctx.stall.iteration = 10;
    ctx.wave_state_path = test_dir_ / "wave_state.json";
    ctx.output_dir = test_dir_ / "diagnostic_output";

    auto result = agent.DiagnoseAndResolve(ctx);

    EXPECT_FALSE(result.resolved);
    EXPECT_EQ(result.report.stall_pattern, "no-progress");
    EXPECT_EQ(result.report.root_cause, "diagnostics-disabled");
}

TEST_F(PipelineAgentTest, GenerateReportWithDisabledReporter) {
    config_.enable_reporter = false;
    PipelineAgent agent(config_);

    ReporterContext ctx;
    ctx.sample_id = "test_sample";
    ctx.alignment_path = test_dir_ / "alignments.sam";
    ctx.wave_state_summary = test_dir_ / "summary.json";
    ctx.output_dir = test_dir_ / "report_output";

    auto report = agent.GenerateReport(ctx);

    EXPECT_EQ(report.sample_id, "test_sample");
    EXPECT_EQ(report.summary, "reporter-disabled");
}

TEST_F(PipelineAgentTest, GetLoadedKernelWithoutSandbox) {
    PipelineAgent agent(config_);

    auto kernel = agent.GetLoadedKernel("test_kernel");
    EXPECT_FALSE(kernel.has_value());
}

// --- Utility function tests ---

TEST_F(PipelineAgentTest, WriteWaveStateJson) {
    std::vector<std::pair<std::string, double>> buckets = {
        {"bucket_0", 0.8},
        {"bucket_1", 0.15},
        {"bucket_2", 0.05}
    };

    auto path = WriteWaveStateJson(test_dir_, 5, 1.2, buckets);

    EXPECT_TRUE(std::filesystem::exists(path));
    EXPECT_EQ(path.filename().string(), "wave_state_iter_5.json");

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("\"iteration\": 5") != std::string::npos);
    EXPECT_TRUE(content.find("\"entropy\": 1.2") != std::string::npos);
    EXPECT_TRUE(content.find("bucket_0") != std::string::npos);
    EXPECT_TRUE(content.find("0.8") != std::string::npos);
}

TEST_F(PipelineAgentTest, WriteWaveStateJsonCreatesDirectory) {
    auto output_dir = test_dir_ / "nested" / "output" / "dir";
    std::vector<std::pair<std::string, double>> buckets;

    auto path = WriteWaveStateJson(output_dir, 0, 0.0, buckets);

    EXPECT_TRUE(std::filesystem::exists(output_dir));
    EXPECT_TRUE(std::filesystem::exists(path));
}

TEST_F(PipelineAgentTest, AnalyzeWaveStateForStall) {
    std::ofstream out(test_dir_ / "wave_state.json");
    out << R"({
        "iteration": 10,
        "entropy": 2.5,
        "bucket_probabilities": []
    })";
    out.close();

    auto metrics = AnalyzeWaveStateForStall(test_dir_ / "wave_state.json");

    EXPECT_EQ(metrics.iteration, 10);
    EXPECT_DOUBLE_EQ(metrics.entropy_delta, 2.5);
}

TEST_F(PipelineAgentTest, AnalyzeWaveStateForStallNonexistent) {
    auto metrics = AnalyzeWaveStateForStall(test_dir_ / "nonexistent.json");
    EXPECT_EQ(metrics.type, StallType::None);
    EXPECT_EQ(metrics.iteration, 0);
}

// --- StallType tests ---

TEST(StallTypeTest, StallTypeNames) {
    EXPECT_STREQ(StallTypeName(StallType::None), "none");
    EXPECT_STREQ(StallTypeName(StallType::NoProgress), "no-progress");
    EXPECT_STREQ(StallTypeName(StallType::OscillatingProbabilities), "oscillating-probabilities");
    EXPECT_STREQ(StallTypeName(StallType::HighUncertainty), "high-uncertainty");
    EXPECT_STREQ(StallTypeName(StallType::ParalogConfusion), "paralog-confusion");
    EXPECT_STREQ(StallTypeName(StallType::CoverageAnomaly), "coverage-anomaly");
}

// --- Integration-style tests ---

class PipelineAgentIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_pipeline_integration";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path test_dir_;
};

TEST_F(PipelineAgentIntegrationTest, SimulatedAlignmentWorkflow) {
    PipelineAgentConfig config;
    config.work_dir = test_dir_;
    config.enable_diagnostics = true;
    config.enable_reporter = true;
    config.enable_cuda_codegen = false;
    config.stall_detection_window = 3;
    config.stall_entropy_threshold = 0.05;

    PipelineAgent agent(config);

    bool stall_detected = false;
    agent.SetStallCallback([&](const StallMetrics&) {
        stall_detected = true;
    });

    agent.RecordIteration(0, 3.0, 1000, {"chr1:1000-50000"});
    agent.RecordIteration(1, 2.5, 800, {"chr1:1000-50000"});
    agent.RecordIteration(2, 2.0, 600, {"chr1:1000-50000"});

    EXPECT_FALSE(stall_detected);
    EXPECT_FALSE(agent.IsStalled());

    agent.RecordIteration(3, 2.0, 600, {"chr1:1000-50000"});
    agent.RecordIteration(4, 2.0, 600, {"chr1:1000-50000"});
    agent.RecordIteration(5, 2.0, 600, {"chr1:1000-50000"});

    EXPECT_TRUE(stall_detected);
    EXPECT_TRUE(agent.IsStalled());

    auto metrics = agent.GetStallMetrics();
    EXPECT_EQ(metrics.type, StallType::NoProgress);
    EXPECT_EQ(metrics.iteration, 5);
    EXPECT_EQ(metrics.reads_affected, 600);
}

TEST_F(PipelineAgentIntegrationTest, WriteAndAnalyzeWaveState) {
    std::vector<std::pair<std::string, double>> buckets = {
        {"IGHV1-2*01", 0.45},
        {"IGHV1-2*02", 0.35},
        {"IGHV1-3*01", 0.15},
        {"other", 0.05}
    };

    auto path = WriteWaveStateJson(test_dir_, 10, 1.8, buckets);
    ASSERT_TRUE(std::filesystem::exists(path));

    auto metrics = AnalyzeWaveStateForStall(path);
    EXPECT_EQ(metrics.iteration, 10);
    EXPECT_NEAR(metrics.entropy_delta, 1.8, 0.001);
}

TEST_F(PipelineAgentIntegrationTest, MultipleSessionsIndependent) {
    PipelineAgentConfig config1;
    config1.work_dir = test_dir_ / "session1";
    config1.stall_detection_window = 3;

    PipelineAgentConfig config2;
    config2.work_dir = test_dir_ / "session2";
    config2.stall_detection_window = 5;

    PipelineAgent agent1(config1);
    PipelineAgent agent2(config2);

    EXPECT_NE(agent1.SessionId(), agent2.SessionId());

    for (int i = 0; i < 4; ++i) {
        agent1.RecordIteration(i, 2.0);
        agent2.RecordIteration(i, 2.0);
    }

    EXPECT_TRUE(agent1.IsStalled());
    EXPECT_FALSE(agent2.IsStalled());

    for (int i = 4; i < 6; ++i) {
        agent2.RecordIteration(i, 2.0);
    }

    EXPECT_TRUE(agent2.IsStalled());
}

}  // namespace
}  // namespace llmap::claude_agent
