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
    auto result = Exec(llmap_bin_ + " align -x /tmp/ref.fa -o /tmp/out.sam");

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
    auto result = Exec(llmap_bin_ + " align -r " + fastq.string() + " -x /tmp/ref.fa");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_TRUE(result.output.find("--output") != std::string::npos ||
                result.output.find("required") != std::string::npos);
}

TEST_F(LlmapCliTest, AlignFileNotFound) {
    auto result = Exec(llmap_bin_ + " align -r /nonexistent/reads.fastq -x /tmp/ref.fa -o /tmp/out.sam");

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
                       " -x " + fasta_path.string() +
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
                       " -x " + fasta_path.string() +
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
                       " -x " + fasta_path.string() +
                       " -o " + output.string() +
                       " --llm --llm-work-dir " + llm_dir.string() +
                       " --llm-threshold 0.99 -v");

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

}  // namespace
