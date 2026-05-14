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

#include "core/version.h"

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
        // Unique directory per test to avoid parallel test conflicts
        auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string test_name = std::string(test_info->test_suite_name()) +
                                "_" + test_info->name();
        test_dir_ = std::filesystem::temp_directory_path() / ("llmap_cli_" + test_name);
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
    EXPECT_TRUE(result.output.find(llmap::kVersion) != std::string::npos);
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

// ========== Index Subcommand ==========

TEST_F(LlmapCliTest, IndexHelp) {
    auto result = Exec(llmap_bin_ + " index --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Build a minimizer index") != std::string::npos);
    EXPECT_TRUE(result.output.find("--reference") != std::string::npos);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos);
    EXPECT_TRUE(result.output.find("--kmer") != std::string::npos);
    EXPECT_TRUE(result.output.find("--window") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexHelpShort) {
    auto result = Exec(llmap_bin_ + " index -h");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--reference") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexMissingReference) {
    auto result = Exec(llmap_bin_ + " index -o /tmp/out.llmi");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--reference") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexMissingOutput) {
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto result = Exec(llmap_bin_ + " index -r " + fasta_path.string());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexFileNotFound) {
    auto result = Exec(llmap_bin_ + " index -r /nonexistent/file.fa -o /tmp/out.llmi");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexUnknownOption) {
    auto result = Exec(llmap_bin_ + " index --unknown-option");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Unknown") != std::string::npos ||
                result.output.find("unknown") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexBasicRun) {
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        // Long enough sequence for minimizer extraction (needs k+w bases)
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        fasta << ">chr2\n";
        fasta << "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
        fasta << "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n";
    }

    auto output_path = test_dir_ / "ref.llmi";

    auto result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                       " -o " + output_path.string() + " -v");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("Index built successfully") != std::string::npos);
    EXPECT_TRUE(result.output.find("Sequences") != std::string::npos);
    EXPECT_TRUE(result.output.find("Minimizers") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST_F(LlmapCliTest, IndexCustomParams) {
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto output_path = test_dir_ / "ref.llmi";

    auto result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                       " -o " + output_path.string() +
                       " -k 15 -w 10 --max-occ 200 -v");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("k-mer size:    15") != std::string::npos);
    EXPECT_TRUE(result.output.find("Window size:   10") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST_F(LlmapCliTest, IndexInvalidKmerSize) {
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto output_path = test_dir_ / "ref.llmi";

    auto result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                       " -o " + output_path.string() + " -k 50");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("k-mer size must be between") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexEmptyFasta) {
    auto fasta_path = test_dir_ / "empty.fa";
    {
        std::ofstream fasta(fasta_path);
        // Empty file
    }

    auto output_path = test_dir_ / "empty.llmi";

    auto result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                       " -o " + output_path.string());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("no sequences") != std::string::npos);
}

TEST_F(LlmapCliTest, IndexVerboseStats) {
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto output_path = test_dir_ / "ref.llmi";

    auto result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                       " -o " + output_path.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Loading reference") != std::string::npos);
    EXPECT_TRUE(result.output.find("Building index") != std::string::npos);
    EXPECT_TRUE(result.output.find("Saving index") != std::string::npos);
}

// ========== Align Subcommand ==========

TEST_F(LlmapCliTest, AlignHelp) {
    auto result = Exec(llmap_bin_ + " align --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Align reads") != std::string::npos);
    EXPECT_TRUE(result.output.find("--reads") != std::string::npos);
    EXPECT_TRUE(result.output.find("--reference") != std::string::npos);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignMissingReads) {
    auto result = Exec(llmap_bin_ + " align --reference /tmp/ref.fa -o /tmp/out.sam");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--reads") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignMissingReference) {
    auto fastq = CreateTestFastq("test.fastq");
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() + " -o /tmp/out.sam");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--reference") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignMissingOutput) {
    auto fastq = CreateTestFastq("test.fastq");
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() + " --reference /tmp/ref.fa");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignFileNotFound) {
    auto result = Exec(llmap_bin_ + " align -r /nonexistent/reads.fastq --reference /tmp/ref.fa -o /tmp/out.sam");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos);
}

// ========== Align --llm Flag ==========

TEST_F(LlmapCliTest, AlignHelpShowsLlmFlag) {
    auto result = Exec(llmap_bin_ + " align --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--llm") != std::string::npos);
    EXPECT_TRUE(result.output.find("LLM-assisted") != std::string::npos ||
                result.output.find("Claude") != std::string::npos);
    EXPECT_TRUE(result.output.find("--llm-api-key") != std::string::npos);
    EXPECT_TRUE(result.output.find("--llm-threshold") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignLlmFlagNoApiKey) {
    // Create test files
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    // Unset env var to ensure no key is found
    unsetenv("ANTHROPIC_API_KEY");

    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --llm -v");

    // Should still complete alignment (LLM is optional), with warning
    // Either success with warning, or success without LLM
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignLlmThreshold) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    // Very high threshold (0.99) - diagnostics should trigger for almost any run
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --llm --llm-threshold 0.99 -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    // Should either trigger diagnostics or warn about missing API key
    EXPECT_TRUE(result.output.find("diagnostics") != std::string::npos ||
                result.output.find("Warning") != std::string::npos ||
                result.output.find("API key") != std::string::npos ||
                result.output.find("threshold") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignLlmWorkDir) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";
    auto llm_dir = test_dir_ / "llm_work";

    // With --llm-work-dir, should create directory
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --llm --llm-work-dir " + llm_dir.string() +
                       " --llm-threshold 0.99 -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

// ========== Align --index Flag (Index Caching) ==========

TEST_F(LlmapCliTest, AlignHelpShowsIndexFlag) {
    auto result = Exec(llmap_bin_ + " align --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--index") != std::string::npos);
    EXPECT_TRUE(result.output.find(".llmi") != std::string::npos);
    EXPECT_TRUE(result.output.find("pre-built") != std::string::npos ||
                result.output.find("Index caching") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignIndexFileNotFound) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --index /nonexistent/index.llmi");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignWithPrebuiltIndex) {
    // Create test files
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto index_path = test_dir_ / "ref.llmi";
    auto output = test_dir_ / "out.sam";

    // First build the index using llmap index
    auto index_result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                             " -o " + index_path.string() + " -k 15 -w 10");
    EXPECT_EQ(index_result.exit_code, 0);
    EXPECT_TRUE(std::filesystem::exists(index_path));

    // Now use the pre-built index for alignment
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " -i " + index_path.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Loading pre-built index") != std::string::npos ||
                result.output.find("Index loaded") != std::string::npos);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignIndexShortFlag) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto index_path = test_dir_ / "ref.llmi";
    auto output = test_dir_ / "out.sam";

    // Build index
    auto index_result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                             " -o " + index_path.string() + " -k 15 -w 10");
    EXPECT_EQ(index_result.exit_code, 0);

    // Test -i short flag
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " -i " + index_path.string());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignIndexVerboseShowsParams) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto index_path = test_dir_ / "ref.llmi";
    auto output = test_dir_ / "out.sam";

    // Build index with specific k and w
    auto index_result = Exec(llmap_bin_ + " index -r " + fasta_path.string() +
                             " -o " + index_path.string() + " -k 19 -w 10");
    EXPECT_EQ(index_result.exit_code, 0);

    // Align with verbose - should show loaded index params
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " -i " + index_path.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    // Should report the k and w from the loaded index
    EXPECT_TRUE(result.output.find("k=19") != std::string::npos);
    EXPECT_TRUE(result.output.find("w=10") != std::string::npos);
}

// ========== Align -x Presets ==========

TEST_F(LlmapCliTest, AlignHelpShowsPresets) {
    auto result = Exec(llmap_bin_ + " align --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("-x PRESET") != std::string::npos);
    EXPECT_TRUE(result.output.find("map-hifi") != std::string::npos);
    EXPECT_TRUE(result.output.find("map-ont") != std::string::npos);
    EXPECT_TRUE(result.output.find("HiFi") != std::string::npos ||
                result.output.find("PacBio") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPresetInvalid) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -x bogus-preset -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Unknown preset") != std::string::npos ||
                result.output.find("bogus-preset") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPresetMapHifi) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -x map-hifi -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    // HiFi preset uses k=19, w=19
    EXPECT_TRUE(result.output.find("k=19") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPresetMapOnt) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -x map-ont -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    // ONT preset uses k=15, w=10
    EXPECT_TRUE(result.output.find("k=15") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPresetSr) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -x sr -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    // SR preset uses k=21, w=11
    EXPECT_TRUE(result.output.find("k=21") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPresetOverrideKmer) {
    // Explicit -k should override preset
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    // map-hifi default is k=19, but we override to k=15
    auto result = Exec(llmap_bin_ + " align -x map-hifi -k 15 -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("k=15") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPresetHifiShortName) {
    // Test "hifi" shorthand works as well as "map-hifi"
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -x hifi -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("k=19") != std::string::npos);
}

// ========== Align --classical-only Mode ==========

TEST_F(LlmapCliTest, AlignHelpShowsClassicalOnlyFlag) {
    auto result = Exec(llmap_bin_ + " align --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--classical-only") != std::string::npos);
    EXPECT_TRUE(result.output.find("seed-chain-extend") != std::string::npos ||
                result.output.find("probabilistic") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignClassicalOnlyBasicRun) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align --classical-only -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(output));
}

TEST_F(LlmapCliTest, AlignClassicalOnlyOverridesLlm) {
    // --classical-only should disable --llm even if both are specified
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    // Unset env var to ensure no key is found
    unsetenv("ANTHROPIC_API_KEY");

    auto result = Exec(llmap_bin_ + " align --classical-only --llm -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    // Should see a note about --llm being ignored
    EXPECT_TRUE(result.output.find("ignored") != std::string::npos ||
                result.output.find("classical-only") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignClassicalOnlyWithPreset) {
    // --classical-only should work in combination with -x presets
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align --classical-only -x map-hifi -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    // Should use HiFi preset k=19
    EXPECT_TRUE(result.output.find("k=19") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignClassicalOnlyReducesMemory) {
    // Running with --classical-only should complete successfully
    // (memory reduction is a behavioral guarantee, not easily testable in unit tests)
    auto fastq = CreateTestFastq("test.fastq", 50);  // More reads
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align --classical-only -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() + " -t 2");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

// ========== Regression: Banner Shows on Empty Args ==========

TEST_F(LlmapCliTest, BannerShowsCorrectTagline) {
    auto result = Exec(llmap_bin_);

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Lossless") != std::string::npos);
    EXPECT_TRUE(result.output.find("Wave-particle") != std::string::npos ||
                result.output.find("wave") != std::string::npos);
}

// ========== sc-paralog-matrix Subcommand ==========

TEST_F(LlmapCliTest, ScParalogMatrixHelp) {
    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("cell") != std::string::npos);
    EXPECT_TRUE(result.output.find("paralog") != std::string::npos);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos);
    EXPECT_TRUE(result.output.find("--cb-tag") != std::string::npos);
    EXPECT_TRUE(result.output.find("--cb-pattern") != std::string::npos);
}

TEST_F(LlmapCliTest, ScParalogMatrixHelpShort) {
    auto result = Exec(llmap_bin_ + " sc-paralog-matrix -h");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos);
}

TEST_F(LlmapCliTest, ScParalogMatrixMissingParquet) {
    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --output /tmp/out.csv");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, ScParalogMatrixMissingOutput) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:AAAA_read1,paralog1,0.8,1.0,0,1,true\n";
    }

    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, ScParalogMatrixFileNotFound) {
    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet /nonexistent/file.csv --output /tmp/out.csv");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos);
}

TEST_F(LlmapCliTest, ScParalogMatrixBasicRun) {
    // Create test CSV with probability entries
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:AAAA_read1,paralog1,0.8,1.0,0,1,true\n";
        csv << "CB:Z:AAAA_read2,paralog1,0.7,1.0,0,1,true\n";
        csv << "CB:Z:AAAA_read3,paralog2,0.9,1.0,0,1,true\n";
        csv << "CB:Z:BBBB_read4,paralog1,0.6,1.0,0,1,true\n";
        csv << "CB:Z:BBBB_read5,paralog2,0.5,1.0,0,1,true\n";
    }

    auto output_path = test_dir_ / "matrix.csv";

    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                       " --output " + output_path.string() + " --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);
    EXPECT_TRUE(result.output.find("Unique cells") != std::string::npos);

    // Check output file exists
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST_F(LlmapCliTest, ScParalogMatrixWithRegex) {
    // Create test CSV with barcode in read name using pattern
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "AAGGCCTT_read1,paralog1,0.8,1.0,0,1,true\n";
        csv << "AAGGCCTT_read2,paralog2,0.7,1.0,0,1,true\n";
        csv << "TTCCAAGG_read3,paralog1,0.9,1.0,0,1,true\n";
    }

    auto output_path = test_dir_ / "matrix.csv";

    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                       " --output " + output_path.string() +
                       " --cb-pattern '([ACGT]{8})_' --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST_F(LlmapCliTest, ScParalogMatrixDenseOutput) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,PARA1,0.8,1.0,0,1,true\n";
        csv << "CB:Z:CELL1_read2,PARA2,0.6,1.0,0,1,true\n";
        csv << "CB:Z:CELL2_read3,PARA1,0.9,1.0,0,1,true\n";
    }

    auto output_path = test_dir_ / "matrix_dense.csv";

    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                       " --output " + output_path.string() + " --dense --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(std::filesystem::exists(output_path));

    // Read and verify dense format (first row is header)
    std::ifstream in(output_path);
    std::string line;
    EXPECT_TRUE(std::getline(in, line));  // Header
    EXPECT_TRUE(line.find("PARA1") != std::string::npos ||
                line.find("PARA2") != std::string::npos);
}

TEST_F(LlmapCliTest, ScParalogMatrixAggregationModes) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,PARA1,0.8,1.0,0,1,true\n";
        csv << "CB:Z:CELL1_read2,PARA1,0.4,1.0,0,1,true\n";  // Same cell+paralog
    }

    auto output_path = test_dir_ / "matrix.csv";

    // Test mean aggregation
    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                       " --output " + output_path.string() +
                       " --aggregation mean --verbose");
    EXPECT_EQ(result.exit_code, 0);

    // Test max aggregation
    result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                  " --output " + output_path.string() +
                  " --aggregation max --verbose");
    EXPECT_EQ(result.exit_code, 0);

    // Test sum aggregation
    result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                  " --output " + output_path.string() +
                  " --aggregation sum --verbose");
    EXPECT_EQ(result.exit_code, 0);
}

TEST_F(LlmapCliTest, ScParalogMatrixNormalize) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,PARA1,0.5,1.0,0,1,true\n";
        csv << "CB:Z:CELL1_read2,PARA2,0.5,1.0,0,1,true\n";
    }

    auto output_path = test_dir_ / "matrix.csv";

    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                       " --output " + output_path.string() +
                       " --normalize --verbose");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST_F(LlmapCliTest, ScParalogMatrixMinReadsFilter) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,PARA1,0.8,1.0,0,1,true\n";
        csv << "CB:Z:CELL1_read2,PARA1,0.7,1.0,0,1,true\n";
        csv << "CB:Z:CELL1_read3,PARA1,0.9,1.0,0,1,true\n";
        csv << "CB:Z:CELL2_read4,PARA2,0.5,1.0,0,1,true\n";  // Only 1 read
    }

    auto output_path = test_dir_ / "matrix.csv";

    // Filter requiring min 2 reads per cell
    auto result = Exec(llmap_bin_ + " sc-paralog-matrix --parquet " + csv_path.string() +
                       " --output " + output_path.string() +
                       " --min-reads 2 --verbose");

    EXPECT_EQ(result.exit_code, 0);
    // CELL2 should be filtered out (only 1 read)
    EXPECT_TRUE(result.output.find("Unique cells") != std::string::npos);
}

// ========== Align --psv-catalog Flag ==========

TEST_F(LlmapCliTest, AlignHelpShowsPsvCatalogFlag) {
    auto result = Exec(llmap_bin_ + " align --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--psv-catalog") != std::string::npos);
    EXPECT_TRUE(result.output.find("--psv-weight") != std::string::npos);
    EXPECT_TRUE(result.output.find("--psv-min-posterior") != std::string::npos);
    EXPECT_TRUE(result.output.find("--psv-only") != std::string::npos);
    EXPECT_TRUE(result.output.find("paralog") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPsvCatalogFileNotFound) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }
    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --psv-catalog /nonexistent/catalog.bed");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos ||
                result.output.find("PSV catalog") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPsvCatalogBasicRun) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto psv_path = test_dir_ / "catalog.bed";
    {
        std::ofstream psv(psv_path);
        psv << "chr1\t10\tA\tparalog1:A,paralog2:G\t0.8\n";
        psv << "chr1\t20\tC\tparalog1:C,paralog2:T\t0.7\n";
    }

    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --psv-catalog " + psv_path.string() +
                       " -v");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
    EXPECT_TRUE(result.output.find("PSV") != std::string::npos ||
                result.output.find("catalog") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPsvWeight) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto psv_path = test_dir_ / "catalog.bed";
    {
        std::ofstream psv(psv_path);
        psv << "chr1\t10\tA\tparalog1:A,paralog2:G\t0.8\n";
    }

    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --psv-catalog " + psv_path.string() +
                       " --psv-weight 0.8 --psv-min-posterior 0.95 -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignPsvOnly) {
    auto fastq = CreateTestFastq("test.fastq");
    auto fasta_path = test_dir_ / "ref.fa";
    {
        std::ofstream fasta(fasta_path);
        fasta << ">chr1\n";
        fasta << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    }

    auto psv_path = test_dir_ / "catalog.bed";
    {
        std::ofstream psv(psv_path);
        psv << "chr1\t10\tA\tparalog1:A,paralog2:G\t0.8\n";
    }

    auto output = test_dir_ / "out.sam";

    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() +
                       " --reference " + fasta_path.string() +
                       " -o " + output.string() +
                       " --psv-catalog " + psv_path.string() +
                       " --psv-only -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("Alignment complete") != std::string::npos);
}

// ========== sc-qc-report Subcommand ==========

TEST_F(LlmapCliTest, ScQcReportHelp) {
    auto result = Exec(llmap_bin_ + " sc-qc-report --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("QC") != std::string::npos);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos);
    EXPECT_TRUE(result.output.find("--qc-json") != std::string::npos);
    EXPECT_TRUE(result.output.find("--qc-tsv") != std::string::npos);
    EXPECT_TRUE(result.output.find("--filtered-matrix") != std::string::npos);
    EXPECT_TRUE(result.output.find("--min-assignment-rate") != std::string::npos);
    EXPECT_TRUE(result.output.find("--min-confidence") != std::string::npos);
    EXPECT_TRUE(result.output.find("--min-reads-per-cell") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportHelpShort) {
    auto result = Exec(llmap_bin_ + " sc-qc-report -h");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportMissingParquet) {
    auto result = Exec(llmap_bin_ + " sc-qc-report --qc-json /tmp/report.json");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--parquet") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportMissingOutput) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:AAAA_read1,paralog1,0.8,1.0,0,1,true\n";
    }

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("output") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportFileNotFound) {
    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet /nonexistent/file.csv --qc-json /tmp/report.json");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("not found") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportBasicJsonOutput) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,paralog1,0.8,0.9,0,1,true\n";
        csv << "CB:Z:CELL1_read2,paralog1,0.7,0.85,0,1,true\n";
        csv << "CB:Z:CELL1_read3,paralog2,0.9,0.95,0,1,true\n";
        csv << "CB:Z:CELL2_read4,paralog1,0.6,0.8,0,1,true\n";
        csv << "CB:Z:CELL2_read5,paralog2,0.5,0.7,0,1,true\n";
    }

    auto json_path = test_dir_ / "report.json";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --qc-json " + json_path.string() + " --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);
    EXPECT_TRUE(result.output.find("Cells passing QC") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(json_path));

    // Verify JSON contains expected fields
    std::ifstream in(json_path);
    std::string json_content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    EXPECT_TRUE(json_content.find("global") != std::string::npos);
    EXPECT_TRUE(json_content.find("cells") != std::string::npos);
    EXPECT_TRUE(json_content.find("paralogs") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportTsvOutput) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,paralog1,0.8,0.9,0,1,true\n";
        csv << "CB:Z:CELL1_read2,paralog2,0.7,0.85,0,1,true\n";
        csv << "CB:Z:CELL2_read3,paralog1,0.6,0.8,0,1,true\n";
    }

    auto tsv_dir = test_dir_ / "qc_report";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --qc-tsv " + tsv_dir.string() + " --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);

    // Check TSV files exist
    EXPECT_TRUE(std::filesystem::exists(tsv_dir / "cells_qc.tsv"));
    EXPECT_TRUE(std::filesystem::exists(tsv_dir / "paralogs_qc.tsv"));
    EXPECT_TRUE(std::filesystem::exists(tsv_dir / "summary_qc.tsv"));
}

TEST_F(LlmapCliTest, ScQcReportFilteredMatrix) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        // CELL1: 3 reads, high confidence - should pass QC
        csv << "CB:Z:CELL1_read1,paralog1,0.9,0.95,0,1,true\n";
        csv << "CB:Z:CELL1_read2,paralog1,0.85,0.9,0,1,true\n";
        csv << "CB:Z:CELL1_read3,paralog2,0.8,0.88,0,1,true\n";
        // CELL2: only 1 read - may fail min reads threshold
        csv << "CB:Z:CELL2_read4,paralog1,0.5,0.6,0,1,true\n";
    }

    auto filtered_path = test_dir_ / "filtered.csv";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --filtered-matrix " + filtered_path.string() +
                       " --min-reads-per-cell 2 --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);
    EXPECT_TRUE(result.output.find("Cells passing QC") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(filtered_path));
}

TEST_F(LlmapCliTest, ScQcReportCustomThresholds) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,paralog1,0.9,0.95,0,1,true\n";
        csv << "CB:Z:CELL1_read2,paralog2,0.8,0.9,0,1,true\n";
    }

    auto json_path = test_dir_ / "report.json";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --qc-json " + json_path.string() +
                       " --min-assignment-rate 0.5 --min-confidence 0.8 "
                       "--min-reads-per-cell 1 --max-entropy 2.0 --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(json_path));
}

TEST_F(LlmapCliTest, ScQcReportSampleId) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,paralog1,0.8,0.9,0,1,true\n";
    }

    auto json_path = test_dir_ / "report.json";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --qc-json " + json_path.string() +
                       " --sample-id test_sample_001 --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(std::filesystem::exists(json_path));

    // Verify sample_id in JSON
    std::ifstream in(json_path);
    std::string json_content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    EXPECT_TRUE(json_content.find("test_sample_001") != std::string::npos);
}

TEST_F(LlmapCliTest, ScQcReportMultipleOutputs) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "CB:Z:CELL1_read1,paralog1,0.9,0.95,0,1,true\n";
        csv << "CB:Z:CELL1_read2,paralog2,0.8,0.9,0,1,true\n";
        csv << "CB:Z:CELL2_read3,paralog1,0.7,0.85,0,1,true\n";
    }

    auto json_path = test_dir_ / "report.json";
    auto tsv_dir = test_dir_ / "qc_tsv";
    auto filtered_path = test_dir_ / "filtered.csv";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --qc-json " + json_path.string() +
                       " --qc-tsv " + tsv_dir.string() +
                       " --filtered-matrix " + filtered_path.string() +
                       " --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(std::filesystem::exists(json_path));
    EXPECT_TRUE(std::filesystem::exists(tsv_dir / "cells_qc.tsv"));
    EXPECT_TRUE(std::filesystem::exists(filtered_path));
}

TEST_F(LlmapCliTest, ScQcReportWithRegexBarcode) {
    auto csv_path = test_dir_ / "probs.csv";
    {
        std::ofstream csv(csv_path);
        csv << "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed\n";
        csv << "AAGGCCTT_read1,paralog1,0.8,0.9,0,1,true\n";
        csv << "AAGGCCTT_read2,paralog2,0.7,0.85,0,1,true\n";
        csv << "TTCCAAGG_read3,paralog1,0.9,0.95,0,1,true\n";
    }

    auto json_path = test_dir_ / "report.json";

    auto result = Exec(llmap_bin_ + " sc-qc-report --parquet " + csv_path.string() +
                       " --qc-json " + json_path.string() +
                       " --cb-pattern '([ACGT]{8})_' --verbose");

    EXPECT_EQ(result.exit_code, 0) << "Output: " << result.output;
    EXPECT_TRUE(result.output.find("complete") != std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(json_path));
}

TEST_F(LlmapCliTest, ScQcReportShowsInHelpMenu) {
    auto result = Exec(llmap_bin_ + " --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("sc-qc-report") != std::string::npos);
}

// ============== Check Command Tests ==============

TEST_F(LlmapCliTest, CheckHelpFlag) {
    auto result = Exec(llmap_bin_ + " check --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("V1.0 readiness check") != std::string::npos);
    EXPECT_TRUE(result.output.find("--verbose") != std::string::npos);
    EXPECT_TRUE(result.output.find("--json") != std::string::npos);
    EXPECT_TRUE(result.output.find("--category") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckBasicRun) {
    auto result = Exec(llmap_bin_ + " check");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("LLmap V1.0 Readiness Check") != std::string::npos);
    EXPECT_TRUE(result.output.find("Result: READY") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckVerboseOutput) {
    auto result = Exec(llmap_bin_ + " check -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("[+] AlignmentRecord") != std::string::npos);
    EXPECT_TRUE(result.output.find("[+] BucketPyramid") != std::string::npos);
    EXPECT_TRUE(result.output.find("[+] WaveState") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckJsonOutput) {
    auto result = Exec(llmap_bin_ + " check -j");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("\"all_passed\": true") != std::string::npos);
    EXPECT_TRUE(result.output.find("\"categories\":") != std::string::npos);
    EXPECT_TRUE(result.output.find("\"checks\":") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckCategoryCoreOnly) {
    auto result = Exec(llmap_bin_ + " check -c Core -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("[Core]") != std::string::npos);
    EXPECT_TRUE(result.output.find("AlignmentRecord") != std::string::npos);
    EXPECT_FALSE(result.output.find("[Foundation]") != std::string::npos);
    EXPECT_FALSE(result.output.find("[Agent]") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckCategoryProductionOnly) {
    auto result = Exec(llmap_bin_ + " check -c Production -v");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("[Production]") != std::string::npos);
    EXPECT_TRUE(result.output.find("Logging") != std::string::npos);
    EXPECT_TRUE(result.output.find("Config") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckInvalidCategory) {
    auto result = Exec(llmap_bin_ + " check -c InvalidCategory");

    EXPECT_EQ(result.exit_code, 1);
    EXPECT_TRUE(result.output.find("Unknown category") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckShowsInMainHelp) {
    auto result = Exec(llmap_bin_ + " --help");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("check") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckReportsAllCategories) {
    auto result = Exec(llmap_bin_ + " check");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("[Core]") != std::string::npos);
    EXPECT_TRUE(result.output.find("[Foundation]") != std::string::npos);
    EXPECT_TRUE(result.output.find("[Classical]") != std::string::npos);
    EXPECT_TRUE(result.output.find("[Production]") != std::string::npos);
}

TEST_F(LlmapCliTest, CheckJsonCategories) {
    auto result = Exec(llmap_bin_ + " check -j");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("\"name\": \"Core\"") != std::string::npos);
    EXPECT_TRUE(result.output.find("\"name\": \"Production\"") != std::string::npos);
}

}  // namespace
