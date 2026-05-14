// Unit tests for LLmap version info (core/version.h, version_util.h)
//
// Tests version constants, format functions, and CLI --version output.

#include <gtest/gtest.h>

#include <array>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <regex>
#include <string>

#include "core/version.h"
#include "core/version_util.h"

namespace {

// Path to the llmap binary
std::string GetLlmapBinary() {
    std::filesystem::path candidates[] = {
        "src/llmap",
        "../build/src/llmap",
        "/home/christian/llmap-local/build/src/llmap"
    };

    for (const auto& candidate : candidates) {
        if (std::filesystem::exists(candidate)) {
            return std::filesystem::absolute(candidate).string();
        }
    }
    return "llmap";
}

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

}  // namespace

// ========== Version Constants ==========

TEST(VersionTest, VersionMajorIsNonNegative) {
    EXPECT_GE(llmap::kVersionMajor, 0);
}

TEST(VersionTest, VersionMinorIsNonNegative) {
    EXPECT_GE(llmap::kVersionMinor, 0);
}

TEST(VersionTest, VersionPatchIsNonNegative) {
    EXPECT_GE(llmap::kVersionPatch, 0);
}

TEST(VersionTest, VersionStringMatchesComponents) {
    std::string expected = std::to_string(llmap::kVersionMajor) + "." +
                           std::to_string(llmap::kVersionMinor) + "." +
                           std::to_string(llmap::kVersionPatch);
    EXPECT_EQ(std::string(llmap::kVersion), expected);
}

TEST(VersionTest, GitCommitIsNotEmpty) {
    EXPECT_FALSE(llmap::kGitCommit.empty());
}

TEST(VersionTest, GitCommitIsValidLength) {
    // Either "unknown" or 8-char hex
    if (llmap::kGitCommit != "unknown") {
        EXPECT_EQ(llmap::kGitCommit.size(), 8);
    }
}

TEST(VersionTest, BuildDateIsISO8601) {
    std::regex iso8601_pattern(R"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)");
    EXPECT_TRUE(std::regex_match(std::string(llmap::kBuildDate), iso8601_pattern));
}

TEST(VersionTest, BuildTypeIsNotEmpty) {
    EXPECT_FALSE(llmap::kBuildType.empty());
}

TEST(VersionTest, CompilerIdIsNotEmpty) {
    EXPECT_FALSE(llmap::kCompilerId.empty());
}

TEST(VersionTest, CompilerVersionIsNotEmpty) {
    EXPECT_FALSE(llmap::kCompilerVersion.empty());
}

TEST(VersionTest, FeaturesStringIsNotEmpty) {
    EXPECT_FALSE(llmap::kFeatures.empty());
}

TEST(VersionTest, FeaturesContainsCuda) {
    std::string features(llmap::kFeatures);
    EXPECT_TRUE(features.find("cuda") != std::string::npos);
}

TEST(VersionTest, FeaturesContainsFaiss) {
    std::string features(llmap::kFeatures);
    EXPECT_TRUE(features.find("faiss") != std::string::npos);
}

TEST(VersionTest, FeaturesContainsClaude) {
    std::string features(llmap::kFeatures);
    EXPECT_TRUE(features.find("claude") != std::string::npos);
}

TEST(VersionTest, HomepageUrlIsValid) {
    EXPECT_EQ(std::string(llmap::kHomepageUrl), "https://losslessmap.com");
}

// ========== Version Utility Functions ==========

TEST(VersionUtilTest, FormatVersionMatchesConstant) {
    EXPECT_EQ(llmap::FormatVersion(), std::string(llmap::kVersion));
}

TEST(VersionUtilTest, FormatVersionShortStartsWithLlmap) {
    std::string short_version = llmap::FormatVersionShort();
    EXPECT_EQ(short_version.substr(0, 6), "llmap ");
}

TEST(VersionUtilTest, FormatVersionShortContainsVersion) {
    std::string short_version = llmap::FormatVersionShort();
    EXPECT_TRUE(short_version.find(llmap::kVersion) != std::string::npos);
}

TEST(VersionUtilTest, FormatVersionFullContainsCommit) {
    std::string full = llmap::FormatVersionFull();
    EXPECT_TRUE(full.find("commit:") != std::string::npos);
    EXPECT_TRUE(full.find(llmap::kGitCommit) != std::string::npos);
}

TEST(VersionUtilTest, FormatVersionFullContainsBuilt) {
    std::string full = llmap::FormatVersionFull();
    EXPECT_TRUE(full.find("built:") != std::string::npos);
}

TEST(VersionUtilTest, FormatVersionFullContainsType) {
    std::string full = llmap::FormatVersionFull();
    EXPECT_TRUE(full.find("type:") != std::string::npos);
}

TEST(VersionUtilTest, FormatVersionFullContainsCompiler) {
    std::string full = llmap::FormatVersionFull();
    EXPECT_TRUE(full.find("compiler:") != std::string::npos);
}

TEST(VersionUtilTest, FormatVersionFullContainsFeatures) {
    std::string full = llmap::FormatVersionFull();
    EXPECT_TRUE(full.find("features:") != std::string::npos);
}

// ========== Feature Query Functions ==========

TEST(VersionFeatureTest, HasCudaMatchesConstant) {
    EXPECT_EQ(llmap::HasCuda(), llmap::kHasCuda);
}

TEST(VersionFeatureTest, HasOnnxRuntimeMatchesConstant) {
    EXPECT_EQ(llmap::HasOnnxRuntime(), llmap::kHasOnnxRuntime);
}

TEST(VersionFeatureTest, HasFaissMatchesConstant) {
    EXPECT_EQ(llmap::HasFaiss(), llmap::kHasFaiss);
}

TEST(VersionFeatureTest, HasFaissGpuMatchesConstant) {
    EXPECT_EQ(llmap::HasFaissGpu(), llmap::kHasFaissGpu);
}

TEST(VersionFeatureTest, HasClaudeMatchesConstant) {
    EXPECT_EQ(llmap::HasClaude(), llmap::kHasClaude);
}

// ========== CLI --version Output ==========

class VersionCliTest : public ::testing::Test {
protected:
    void SetUp() override {
        llmap_bin_ = GetLlmapBinary();
    }

    std::string llmap_bin_;
};

TEST_F(VersionCliTest, VersionFlagExitCode) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_EQ(result.exit_code, 0);
}

TEST_F(VersionCliTest, VersionFlagContainsLlmap) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find("llmap") != std::string::npos);
}

TEST_F(VersionCliTest, VersionFlagContainsVersionNumber) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find(llmap::kVersion) != std::string::npos);
}

TEST_F(VersionCliTest, VersionFlagContainsCommit) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find("commit:") != std::string::npos);
}

TEST_F(VersionCliTest, VersionFlagContainsBuilt) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find("built:") != std::string::npos);
}

TEST_F(VersionCliTest, VersionFlagContainsType) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find("type:") != std::string::npos);
}

TEST_F(VersionCliTest, VersionFlagContainsCompiler) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find("compiler:") != std::string::npos);
}

TEST_F(VersionCliTest, VersionFlagContainsFeatures) {
    auto result = Exec(llmap_bin_ + " --version");
    EXPECT_TRUE(result.output.find("features:") != std::string::npos);
}

TEST_F(VersionCliTest, ShortVersionFlagWorks) {
    auto result = Exec(llmap_bin_ + " -V");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("llmap") != std::string::npos);
}

TEST_F(VersionCliTest, ShortVersionFlagMatchesLong) {
    auto long_result = Exec(llmap_bin_ + " --version");
    auto short_result = Exec(llmap_bin_ + " -V");

    EXPECT_EQ(long_result.exit_code, short_result.exit_code);
    EXPECT_EQ(long_result.output, short_result.output);
}
