#include <gtest/gtest.h>

#include "claude_agent/cuda_sandbox.h"

#include <filesystem>
#include <fstream>
#include <string>

namespace llmap::claude_agent {
namespace {

class CudaAnalyzerTest : public ::testing::Test {
protected:
    CudaAnalyzer analyzer_;
};

TEST_F(CudaAnalyzerTest, AnalyzeSafeKernel) {
    const char* source = R"(
__global__ void add_kernel(float* a, float* b, float* c, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        c[idx] = a[idx] + b[idx];
    }
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.safe);
    EXPECT_EQ(result.violations.size(), 0);
    EXPECT_EQ(result.kernel_count, 1);
}

TEST_F(CudaAnalyzerTest, DetectsFileIO) {
    const char* source = R"(
__global__ void bad_kernel() {
    FILE* f = fopen("/etc/passwd", "r");
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);
    EXPECT_GE(result.violations.size(), 1);

    bool found_fileio = false;
    for (const auto& v : result.violations) {
        if (v.type == SecurityViolation::FileIODetected) {
            found_fileio = true;
            break;
        }
    }
    EXPECT_TRUE(found_fileio);
}

TEST_F(CudaAnalyzerTest, DetectsNetworkCalls) {
    const char* source = R"(
__device__ void connect_to_server() {
    int sock = socket(AF_INET, SOCK_STREAM, 0);
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);

    bool found_network = false;
    for (const auto& v : result.violations) {
        if (v.type == SecurityViolation::NetworkDetected) {
            found_network = true;
            break;
        }
    }
    EXPECT_TRUE(found_network);
}

TEST_F(CudaAnalyzerTest, DetectsSyscalls) {
    const char* source = R"(
__global__ void evil_kernel() {
    fork();
    execve("/bin/sh", nullptr, nullptr);
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);

    bool found_syscall = false;
    for (const auto& v : result.violations) {
        if (v.type == SecurityViolation::SyscallDetected) {
            found_syscall = true;
            break;
        }
    }
    EXPECT_TRUE(found_syscall);
}

TEST_F(CudaAnalyzerTest, DetectsInlineAssembly) {
    const char* source = R"(
__global__ void asm_kernel() {
    asm volatile("mov.u32 %0, %%laneid;" : "=r"(x));
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);

    bool found_asm = false;
    for (const auto& v : result.violations) {
        if (v.type == SecurityViolation::AssemblyBlock) {
            found_asm = true;
            break;
        }
    }
    EXPECT_TRUE(found_asm);
}

TEST_F(CudaAnalyzerTest, AllowsCudaIntrinsics) {
    const char* source = R"(
__global__ void sync_kernel(float* data, int n) {
    __shared__ float shared_data[256];
    int idx = threadIdx.x;
    shared_data[idx] = data[idx];
    __syncthreads();
    data[idx] = shared_data[255 - idx];
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.safe);
    EXPECT_EQ(result.violations.size(), 0);
}

TEST_F(CudaAnalyzerTest, AllowsMathFunctions) {
    const char* source = R"(
__global__ void math_kernel(float* input, float* output, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float x = input[idx];
        output[idx] = sqrtf(x) * expf(-x) + sinf(x) * cosf(x);
    }
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.safe);
    EXPECT_EQ(result.violations.size(), 0);
}

TEST_F(CudaAnalyzerTest, AllowsAtomics) {
    const char* source = R"(
__global__ void atomic_kernel(int* counter) {
    atomicAdd(counter, 1);
    atomicMax(counter, 100);
    atomicCAS(counter, 100, 0);
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.safe);
    EXPECT_EQ(result.violations.size(), 0);
}

TEST_F(CudaAnalyzerTest, DetectsExcessiveLines) {
    std::string source = "__global__ void big_kernel() {\n";
    for (int i = 0; i < 6000; ++i) {
        source += "    int x" + std::to_string(i) + " = " + std::to_string(i) + ";\n";
    }
    source += "}\n";

    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);

    bool found_complexity = false;
    for (const auto& v : result.violations) {
        if (v.type == SecurityViolation::ExcessiveComplexity) {
            found_complexity = true;
            break;
        }
    }
    EXPECT_TRUE(found_complexity);
}

TEST_F(CudaAnalyzerTest, CountsKernels) {
    const char* source = R"(
__global__ void kernel1() {}
__global__ void kernel2() {}
__global__ void kernel3() {}
__device__ void device_func() {}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_EQ(result.kernel_count, 3);
}

TEST_F(CudaAnalyzerTest, ExtractsSymbols) {
    const char* source = R"(
__global__ void test_kernel(float* data) {
    int idx = threadIdx.x;
    float val = sqrtf(data[idx]);
    __syncthreads();
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.symbols_used.contains("test_kernel") ||
                result.symbols_used.contains("sqrtf"));
}

TEST_F(CudaAnalyzerTest, IgnoresComments) {
    const char* source = R"(
// fopen is mentioned in a comment
/* fork() is also mentioned here
   system("rm -rf /") is too
*/
__global__ void safe_kernel() {
    // This is fine
    int x = 1;
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.safe);
}

TEST_F(CudaAnalyzerTest, DetectsSystemCall) {
    const char* source = R"(
__host__ void bad_host() {
    system("rm -rf /");
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);
}

TEST_F(CudaAnalyzerTest, DetectsDynamicLinking) {
    const char* source = R"(
__host__ void load_lib() {
    void* handle = dlopen("evil.so", RTLD_NOW);
    void* sym = dlsym(handle, "payload");
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);
}

TEST_F(CudaAnalyzerTest, DetectsMemoryMapping) {
    const char* source = R"(
__host__ void map_memory() {
    void* ptr = mmap(nullptr, 4096, PROT_EXEC, MAP_PRIVATE, -1, 0);
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_FALSE(result.safe);
}

TEST_F(CudaAnalyzerTest, AllowsCooperativeGroups) {
    const char* source = R"(
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

__global__ void cg_kernel() {
    auto block = cg::this_thread_block();
    cg::sync(block);
}
)";
    auto result = analyzer_.Analyze(source);
    EXPECT_TRUE(result.safe);
}

TEST_F(CudaAnalyzerTest, MeasuresAnalysisTime) {
    const char* source = "__global__ void empty() {}";
    auto result = analyzer_.Analyze(source);
    EXPECT_GE(result.analysis_time.count(), 0);
}

class CudaAnalyzerConfigTest : public ::testing::Test {
protected:
    void SetUp() override {
        config_.max_lines = 100;
        config_.max_functions = 5;
        config_.max_kernels = 2;
    }

    AnalyzerConfig config_;
};

TEST_F(CudaAnalyzerConfigTest, EnforcesMaxKernels) {
    CudaAnalyzer analyzer(config_);

    const char* source = R"(
__global__ void k1() {}
__global__ void k2() {}
__global__ void k3() {}
)";
    auto result = analyzer.Analyze(source);
    EXPECT_FALSE(result.safe);

    bool found_complexity = false;
    for (const auto& v : result.violations) {
        if (v.type == SecurityViolation::ExcessiveComplexity &&
            v.description.find("kernels") != std::string::npos) {
            found_complexity = true;
        }
    }
    EXPECT_TRUE(found_complexity);
}

TEST_F(CudaAnalyzerConfigTest, GetConfigReturnsCorrectValues) {
    CudaAnalyzer analyzer(config_);
    EXPECT_EQ(analyzer.GetConfig().max_lines, 100);
    EXPECT_EQ(analyzer.GetConfig().max_functions, 5);
    EXPECT_EQ(analyzer.GetConfig().max_kernels, 2);
}

class SecurityViolationTest : public ::testing::Test {};

TEST_F(SecurityViolationTest, AllNamesAreValid) {
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::None), "none");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::SyscallDetected), "syscall-detected");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::FileIODetected), "file-io-detected");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::NetworkDetected), "network-detected");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::ForbiddenSymbol), "forbidden-symbol");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::MalformedCode), "malformed-code");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::ExcessiveComplexity), "excessive-complexity");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::AssemblyBlock), "assembly-block");
    EXPECT_STREQ(SecurityViolationName(SecurityViolation::UnsafePointer), "unsafe-pointer");
}

class CompileStatusTest : public ::testing::Test {};

TEST_F(CompileStatusTest, AllNamesAreValid) {
    EXPECT_STREQ(CompileStatusName(CompileStatus::Success), "success");
    EXPECT_STREQ(CompileStatusName(CompileStatus::CompileError), "compile-error");
    EXPECT_STREQ(CompileStatusName(CompileStatus::Timeout), "timeout");
    EXPECT_STREQ(CompileStatusName(CompileStatus::SecurityViolation), "security-violation");
    EXPECT_STREQ(CompileStatusName(CompileStatus::SandboxFailure), "sandbox-failure");
    EXPECT_STREQ(CompileStatusName(CompileStatus::LinkError), "link-error");
    EXPECT_STREQ(CompileStatusName(CompileStatus::ResourceExhausted), "resource-exhausted");
}

class LoadStatusTest : public ::testing::Test {};

TEST_F(LoadStatusTest, AllNamesAreValid) {
    EXPECT_STREQ(LoadStatusName(LoadStatus::Success), "success");
    EXPECT_STREQ(LoadStatusName(LoadStatus::FileNotFound), "file-not-found");
    EXPECT_STREQ(LoadStatusName(LoadStatus::InvalidModule), "invalid-module");
    EXPECT_STREQ(LoadStatusName(LoadStatus::SymbolNotFound), "symbol-not-found");
    EXPECT_STREQ(LoadStatusName(LoadStatus::BudgetExceeded), "budget-exceeded");
    EXPECT_STREQ(LoadStatusName(LoadStatus::DeviceError), "device-error");
}

class SourceHashTest : public ::testing::Test {};

TEST_F(SourceHashTest, ProducesConsistentHash) {
    const char* source = "__global__ void test() {}";
    std::string hash1 = ComputeSourceHash(source);
    std::string hash2 = ComputeSourceHash(source);
    EXPECT_EQ(hash1, hash2);
}

TEST_F(SourceHashTest, DifferentSourcesDifferentHashes) {
    std::string hash1 = ComputeSourceHash("__global__ void a() {}");
    std::string hash2 = ComputeSourceHash("__global__ void b() {}");
    EXPECT_NE(hash1, hash2);
}

TEST_F(SourceHashTest, HashIsHexString) {
    std::string hash = ComputeSourceHash("test");
    EXPECT_EQ(hash.size(), 16);
    for (char c : hash) {
        EXPECT_TRUE((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'));
    }
}

class AllowedSymbolsTest : public ::testing::Test {};

TEST_F(AllowedSymbolsTest, ContainsCudaIntrinsics) {
    const auto& symbols = CudaAnalyzer::AllowedSymbols();
    EXPECT_TRUE(symbols.contains("__syncthreads"));
    EXPECT_TRUE(symbols.contains("atomicAdd"));
    EXPECT_TRUE(symbols.contains("__global__"));
}

TEST_F(AllowedSymbolsTest, ContainsMathFunctions) {
    const auto& symbols = CudaAnalyzer::AllowedSymbols();
    EXPECT_TRUE(symbols.contains("sqrtf"));
    EXPECT_TRUE(symbols.contains("expf"));
    EXPECT_TRUE(symbols.contains("sinf"));
}

TEST_F(AllowedSymbolsTest, ContainsThreadIndexing) {
    const auto& symbols = CudaAnalyzer::AllowedSymbols();
    EXPECT_TRUE(symbols.contains("threadIdx"));
    EXPECT_TRUE(symbols.contains("blockIdx"));
    EXPECT_TRUE(symbols.contains("blockDim"));
}

class ForbiddenPatternsTest : public ::testing::Test {};

TEST_F(ForbiddenPatternsTest, ContainsFileIO) {
    const auto& patterns = CudaAnalyzer::ForbiddenPatterns();
    EXPECT_TRUE(patterns.contains("fopen"));
    EXPECT_TRUE(patterns.contains("fread"));
    EXPECT_TRUE(patterns.contains("fwrite"));
}

TEST_F(ForbiddenPatternsTest, ContainsNetworking) {
    const auto& patterns = CudaAnalyzer::ForbiddenPatterns();
    EXPECT_TRUE(patterns.contains("socket"));
    EXPECT_TRUE(patterns.contains("connect"));
    EXPECT_TRUE(patterns.contains("bind"));
}

TEST_F(ForbiddenPatternsTest, ContainsSyscalls) {
    const auto& patterns = CudaAnalyzer::ForbiddenPatterns();
    EXPECT_TRUE(patterns.contains("fork"));
    EXPECT_TRUE(patterns.contains("execve"));
    EXPECT_TRUE(patterns.contains("system"));
}

class CudaLoaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        budget_.max_runtime = std::chrono::milliseconds(1000);
        budget_.max_registers_per_thread = 64;
        budget_.max_shared_memory_bytes = 48 * 1024;
    }

    KernelBudget budget_;
};

TEST_F(CudaLoaderTest, DefaultBudgetValues) {
    CudaLoader loader;
    EXPECT_EQ(loader.GetBudget().max_runtime.count(), 1000);
    EXPECT_EQ(loader.GetBudget().max_registers_per_thread, 64);
}

TEST_F(CudaLoaderTest, CustomBudget) {
    CudaLoader loader(budget_);
    EXPECT_EQ(loader.GetBudget().max_runtime.count(), 1000);
    EXPECT_EQ(loader.GetBudget().max_shared_memory_bytes, 48 * 1024);
}

TEST_F(CudaLoaderTest, InitiallyNoModules) {
    CudaLoader loader;
    EXPECT_EQ(loader.LoadedModuleCount(), 0);
    EXPECT_TRUE(loader.LoadedModuleIds().empty());
}

TEST_F(CudaLoaderTest, LoadNonexistentFile) {
    CudaLoader loader;
    auto result = loader.Load("/nonexistent/path/kernel.cubin");
    EXPECT_EQ(result.status, LoadStatus::FileNotFound);
}

TEST_F(CudaLoaderTest, UnloadAllClearsModules) {
    CudaLoader loader;
    loader.UnloadAll();
    EXPECT_EQ(loader.LoadedModuleCount(), 0);
}

TEST_F(CudaLoaderTest, GetKernelFromUnknownModule) {
    CudaLoader loader;
    auto kernel = loader.GetKernel("unknown_module", "kernel");
    EXPECT_FALSE(kernel.has_value());
}

class CudaCompilerTest : public ::testing::Test {
protected:
    void SetUp() override {
        temp_dir_ = std::filesystem::temp_directory_path() / "llmap_test_cuda";
        std::filesystem::create_directories(temp_dir_);

        config_.sandbox_dir = temp_dir_ / "sandbox";
        config_.output_dir = temp_dir_ / "output";
        config_.nvcc_path = "/usr/local/cuda/bin/nvcc";
        config_.timeout_seconds = 30;
        config_.use_bubblewrap = false;
    }

    void TearDown() override {
        std::filesystem::remove_all(temp_dir_);
    }

    std::filesystem::path temp_dir_;
    CompileConfig config_;
};

TEST_F(CudaCompilerTest, GetConfigReturnsCorrectValues) {
    CudaCompiler compiler(config_);
    EXPECT_EQ(compiler.GetConfig().timeout_seconds, 30);
    EXPECT_EQ(compiler.GetConfig().use_bubblewrap, false);
}

TEST_F(CudaCompilerTest, CreatesSandboxDirectories) {
    CudaCompiler compiler(config_);
    EXPECT_TRUE(std::filesystem::exists(config_.sandbox_dir));
    EXPECT_TRUE(std::filesystem::exists(config_.output_dir));
}

TEST_F(CudaCompilerTest, BuildsBubblewrapCommand) {
    config_.use_bubblewrap = true;
    CudaCompiler compiler(config_);

    auto cmd = compiler.BuildBubblewrapCommand(
        config_.sandbox_dir / "test.cu",
        config_.output_dir / "test.cubin"
    );

    EXPECT_TRUE(cmd.find("bwrap") != std::string::npos);
    EXPECT_TRUE(cmd.find("--unshare-net") != std::string::npos);
    EXPECT_TRUE(cmd.find("--ro-bind") != std::string::npos);
}

TEST_F(CudaCompilerTest, CompileNonexistentFile) {
    CudaCompiler compiler(config_);
    auto result = compiler.CompileFile("/nonexistent/file.cu", "test");
    EXPECT_EQ(result.status, CompileStatus::CompileError);
    EXPECT_TRUE(result.stderr_log.find("does not exist") != std::string::npos);
}

class CudaSandboxConfigTest : public ::testing::Test {
protected:
    void SetUp() override {
        temp_dir_ = std::filesystem::temp_directory_path() / "llmap_sandbox_test";
        std::filesystem::create_directories(temp_dir_);

        config_.compiler.sandbox_dir = temp_dir_ / "sandbox";
        config_.compiler.output_dir = temp_dir_ / "output";
        config_.compiler.nvcc_path = "/usr/local/cuda/bin/nvcc";
        config_.compiler.timeout_seconds = 30;
        config_.compiler.use_bubblewrap = false;
        config_.audit_log_path = temp_dir_ / "audit.log";
        config_.session_id = "test_session";
        config_.require_analysis_pass = true;
    }

    void TearDown() override {
        std::filesystem::remove_all(temp_dir_);
    }

    std::filesystem::path temp_dir_;
    CudaSandbox::Config config_;
};

TEST_F(CudaSandboxConfigTest, GetConfigReturnsCorrectValues) {
    CudaSandbox sandbox(config_);
    EXPECT_EQ(sandbox.GetConfig().session_id, "test_session");
    EXPECT_TRUE(sandbox.GetConfig().require_analysis_pass);
}

TEST_F(CudaSandboxConfigTest, AuditLogInitiallyEmpty) {
    CudaSandbox sandbox(config_);
    EXPECT_TRUE(sandbox.ReadAuditLog().empty());
}

TEST_F(CudaSandboxConfigTest, WriteAndReadAuditEntry) {
    CudaSandbox sandbox(config_);

    AuditEntry entry{
        .timestamp = std::chrono::system_clock::now(),
        .session_id = "test",
        .action = "compile",
        .allowed = true,
        .source_hash = "abc123",
        .details = "test details"
    };

    sandbox.WriteAuditEntry(entry);

    auto log = sandbox.ReadAuditLog();
    EXPECT_EQ(log.size(), 1);
    EXPECT_EQ(log[0].session_id, "test");
    EXPECT_EQ(log[0].action, "compile");
    EXPECT_TRUE(log[0].allowed);
}

TEST_F(CudaSandboxConfigTest, ClearAuditLog) {
    CudaSandbox sandbox(config_);

    AuditEntry entry{
        .timestamp = std::chrono::system_clock::now(),
        .session_id = "test",
        .action = "test",
        .allowed = true
    };

    sandbox.WriteAuditEntry(entry);
    EXPECT_EQ(sandbox.ReadAuditLog().size(), 1);

    sandbox.ClearAuditLog();
    EXPECT_TRUE(sandbox.ReadAuditLog().empty());
}

TEST_F(CudaSandboxConfigTest, RejectsUnsafeCode) {
    CudaSandbox sandbox(config_);

    const char* evil_code = R"(
__global__ void evil_kernel() {
    system("rm -rf /");
}
)";

    auto result = sandbox.CompileAndLoad(evil_code, "evil");
    EXPECT_FALSE(result.success);
    EXPECT_FALSE(result.analysis.safe);
    EXPECT_TRUE(result.error_message.find("analysis failed") != std::string::npos ||
                result.error_message.find("Static analysis failed") != std::string::npos);
}

TEST_F(CudaSandboxConfigTest, GetKernelFromUnloadedModule) {
    CudaSandbox sandbox(config_);
    auto kernel = sandbox.GetKernel("nonexistent");
    EXPECT_FALSE(kernel.has_value());
}

class KernelBudgetTest : public ::testing::Test {};

TEST_F(KernelBudgetTest, DefaultValues) {
    KernelBudget budget;
    EXPECT_EQ(budget.max_runtime.count(), 1000);
    EXPECT_EQ(budget.max_registers_per_thread, 64);
    EXPECT_EQ(budget.max_shared_memory_bytes, 48 * 1024);
    EXPECT_EQ(budget.max_threads_per_block, 1024);
    EXPECT_EQ(budget.max_blocks, 65535);
}

class AnalysisResultTest : public ::testing::Test {};

TEST_F(AnalysisResultTest, DefaultValues) {
    AnalysisResult result;
    EXPECT_FALSE(result.safe);
    EXPECT_TRUE(result.violations.empty());
    EXPECT_EQ(result.lines_analyzed, 0);
    EXPECT_EQ(result.functions_found, 0);
    EXPECT_EQ(result.kernel_count, 0);
}

class CompileResultTest : public ::testing::Test {};

TEST_F(CompileResultTest, DefaultValues) {
    CompileResult result;
    EXPECT_EQ(result.status, CompileStatus::CompileError);
    EXPECT_FALSE(result.output_path.has_value());
    EXPECT_TRUE(result.stdout_log.empty());
    EXPECT_TRUE(result.stderr_log.empty());
}

class LoadResultTest : public ::testing::Test {};

TEST_F(LoadResultTest, DefaultValues) {
    LoadResult result;
    EXPECT_EQ(result.status, LoadStatus::FileNotFound);
    EXPECT_TRUE(result.kernels.empty());
    EXPECT_TRUE(result.module_id.empty());
}

class AnalyzerFileTest : public ::testing::Test {
protected:
    void SetUp() override {
        temp_dir_ = std::filesystem::temp_directory_path() / "llmap_analyzer_test";
        std::filesystem::create_directories(temp_dir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(temp_dir_);
    }

    std::filesystem::path temp_dir_;
    CudaAnalyzer analyzer_;
};

TEST_F(AnalyzerFileTest, AnalyzeValidFile) {
    auto path = temp_dir_ / "safe.cu";
    std::ofstream out(path);
    out << "__global__ void safe_kernel() {}\n";
    out.close();

    auto result = analyzer_.AnalyzeFile(path);
    EXPECT_TRUE(result.safe);
    EXPECT_EQ(result.kernel_count, 1);
}

TEST_F(AnalyzerFileTest, AnalyzeNonexistentFile) {
    auto result = analyzer_.AnalyzeFile("/nonexistent/path.cu");
    EXPECT_FALSE(result.safe);
    EXPECT_FALSE(result.violations.empty());
    EXPECT_EQ(result.violations[0].type, SecurityViolation::MalformedCode);
}

class MoveConstructorTest : public ::testing::Test {};

TEST_F(MoveConstructorTest, CudaAnalyzerMoveWorks) {
    CudaAnalyzer analyzer1;
    CudaAnalyzer analyzer2 = std::move(analyzer1);

    auto result = analyzer2.Analyze("__global__ void test() {}");
    EXPECT_TRUE(result.safe);
}

TEST_F(MoveConstructorTest, CudaLoaderMoveWorks) {
    CudaLoader loader1;
    CudaLoader loader2 = std::move(loader1);

    EXPECT_EQ(loader2.LoadedModuleCount(), 0);
}

}  // namespace
}  // namespace llmap::claude_agent
