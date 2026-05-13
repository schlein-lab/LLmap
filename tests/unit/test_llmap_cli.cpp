// Unit tests for LLmap CLI (llmap_main.cpp)
//
// Tests CLI argument parsing, help output, and error handling.
// Full integration tests require FAISS and are in tests/integration/.

#include <gtest/gtest.h>

#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

namespace {

// Path to the llmap binary
std::string GetLlmapBinary() {
    // Try to find in build directory
    std::filesystem::path candidates[] = {
        "src/llmap",                           // In-tree build (from build/)
        "../build/src/llmap",                  // Out-of-tree (from tests/)
        "/home/christian/llmap-local/build/src/llmap"  // Absolute fallback
    };

    for (const auto& candidate : candidates) {
        if (std::filesystem::exists(candidate)) {
            return std::filesystem::absolute(candidate).string();
        }
    }
    return "llmap";  // Hope it's in PATH
}

// Execute a command and capture output
struct ExecResult {
    int exit_code;
    std::string output;
};

ExecResult Exec(const std::string& cmd) {
    std::array<char, 4096> buffer;
    std::string output;

    std::string full_cmd = cmd + " 2>&1";
    std::unique_ptr<FILE, decltype(&pclose)> pipe(
        popen(full_cmd.c_str(), "r"), pclose);

    if (!pipe) {
        return {-1, "popen failed"};
    }

    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        output += buffer.data();
    }

    int status = pclose(pipe.release());
    int exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;

    return {exit_code, output};
}

class LlmapCliTest : public ::testing::Test {
protected:
    void SetUp() override {
        llmap_bin_ = GetLlmapBinary();
        test_dir_ = std::filesystem::temp_directory_path() / "llmap_cli_test";
        std::filesystem::create_directories(test_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(test_dir_);
    }

    std::filesystem::path CreateTestFastq(const std::string& name, size_t num_reads = 10) {
        auto path = test_dir_ / name;
        std::ofstream out(path);

        for (size_t i = 0; i < num_reads; ++i) {
            out << "@read" << i << "\n";
            out << "ACGTACGTACGTACGTACGTACGTACGTACGT\n";
            out << "+\n";
            out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        }

        return path;
    }

    std::string llmap_bin_;
    std::filesystem::path test_dir_;
};

// ========== Version and Help ==========

TEST_F(LlmapCliTest, Version) {
    auto result = Exec(llmap_bin_ + " --version");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("llmap") != std::string::npos);
    EXPECT_TRUE(result.output.find("0.1.0") != std::string::npos);
}

TEST_F(LlmapCliTest, Help) {
    auto result = Exec(llmap_bin_ + " --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Usage:") != std::string::npos);
    EXPECT_TRUE(result.output.find("allpair") != std::string::npos);
    EXPECT_TRUE(result.output.find("index") != std::string::npos);
    EXPECT_TRUE(result.output.find("align") != std::string::npos);
}

TEST_F(LlmapCliTest, HelpShort) {
    auto result = Exec(llmap_bin_ + " -h");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Usage:") != std::string::npos);
}

TEST_F(LlmapCliTest, NoArgs) {
    auto result = Exec(llmap_bin_);

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("LLmap") != std::string::npos);
    EXPECT_TRUE(result.output.find("Usage:") != std::string::npos);
}

// ========== Unknown Command ==========

TEST_F(LlmapCliTest, UnknownCommand) {
    auto result = Exec(llmap_bin_ + " foobar");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not yet implemented") != std::string::npos ||
                result.output.find("foobar") != std::string::npos);
}

// ========== Allpair Subcommand ==========

TEST_F(LlmapCliTest, AllpairHelp) {
    auto result = Exec(llmap_bin_ + " allpair --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Stage 1 Self-Interference") != std::string::npos);
    EXPECT_TRUE(result.output.find("--reads") != std::string::npos);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos);
    EXPECT_TRUE(result.output.find("--knn") != std::string::npos);
    EXPECT_TRUE(result.output.find("--resolution") != std::string::npos);
}

TEST_F(LlmapCliTest, AllpairHelpShort) {
    auto result = Exec(llmap_bin_ + " allpair -h");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--reads") != std::string::npos);
}

TEST_F(LlmapCliTest, AllpairMissingReads) {
    auto result = Exec(llmap_bin_ + " allpair -o /tmp/out.tsv");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--reads") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, AllpairMissingOutput) {
    auto fastq = CreateTestFastq("test.fastq");
    auto result = Exec(llmap_bin_ + " allpair -r " + fastq.string());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, AllpairFileNotFound) {
    auto result = Exec(llmap_bin_ + " allpair -r /nonexistent/file.fastq -o /tmp/out.tsv");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos);
}

TEST_F(LlmapCliTest, AllpairUnknownOption) {
    auto result = Exec(llmap_bin_ + " allpair --unknown-option");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Unknown") != std::string::npos ||
                result.output.find("unknown") != std::string::npos);
}

// ========== Other Subcommands (Not Yet Implemented) ==========

TEST_F(LlmapCliTest, IndexNotImplemented) {
    auto result = Exec(llmap_bin_ + " index");

    // Should return non-zero with "not yet implemented" message
    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not yet implemented") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignNotImplemented) {
    auto result = Exec(llmap_bin_ + " align");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not yet implemented") != std::string::npos);
}

// ========== Regression: Banner Shows on Empty Args ==========

TEST_F(LlmapCliTest, BannerShowsCorrectTagline) {
    auto result = Exec(llmap_bin_);

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Lossless") != std::string::npos);
    EXPECT_TRUE(result.output.find("Wave-particle") != std::string::npos ||
                result.output.find("wave") != std::string::npos);
}

}  // namespace
