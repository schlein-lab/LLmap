#pragma once

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace llmap::claude_agent {

enum class SecurityViolation : uint8_t {
    None = 0,
    SyscallDetected = 1,
    FileIODetected = 2,
    NetworkDetected = 3,
    ForbiddenSymbol = 4,
    MalformedCode = 5,
    ExcessiveComplexity = 6,
    AssemblyBlock = 7,
    UnsafePointer = 8
};

constexpr const char* SecurityViolationName(SecurityViolation v) {
    switch (v) {
        case SecurityViolation::None:               return "none";
        case SecurityViolation::SyscallDetected:    return "syscall-detected";
        case SecurityViolation::FileIODetected:     return "file-io-detected";
        case SecurityViolation::NetworkDetected:    return "network-detected";
        case SecurityViolation::ForbiddenSymbol:    return "forbidden-symbol";
        case SecurityViolation::MalformedCode:      return "malformed-code";
        case SecurityViolation::ExcessiveComplexity:return "excessive-complexity";
        case SecurityViolation::AssemblyBlock:      return "assembly-block";
        case SecurityViolation::UnsafePointer:      return "unsafe-pointer";
    }
    return "unknown";
}

struct AnalysisViolation {
    SecurityViolation type{SecurityViolation::None};
    size_t line{0};
    size_t column{0};
    std::string description;
    std::string context;
};

struct AnalysisResult {
    bool safe{false};
    std::vector<AnalysisViolation> violations;
    size_t lines_analyzed{0};
    size_t functions_found{0};
    size_t kernel_count{0};
    std::set<std::string> symbols_used;
    std::chrono::milliseconds analysis_time{0};
};

struct AnalyzerConfig {
    size_t max_lines{5000};
    size_t max_functions{50};
    size_t max_kernels{10};
    bool allow_device_functions{true};
    bool allow_shared_memory{true};
    bool strict_pointer_checks{true};
};

struct CompileConfig {
    std::filesystem::path nvcc_path;
    std::filesystem::path sandbox_dir;
    std::filesystem::path output_dir;
    size_t timeout_seconds{60};
    size_t memory_limit_mb{512};
    bool use_bubblewrap{true};
    std::string arch{"sm_70"};
    std::vector<std::string> extra_flags;
};

enum class CompileStatus : uint8_t {
    Success = 0,
    CompileError = 1,
    Timeout = 2,
    SecurityViolation = 3,
    SandboxFailure = 4,
    LinkError = 5,
    ResourceExhausted = 6
};

constexpr const char* CompileStatusName(CompileStatus s) {
    switch (s) {
        case CompileStatus::Success:           return "success";
        case CompileStatus::CompileError:      return "compile-error";
        case CompileStatus::Timeout:           return "timeout";
        case CompileStatus::SecurityViolation: return "security-violation";
        case CompileStatus::SandboxFailure:    return "sandbox-failure";
        case CompileStatus::LinkError:         return "link-error";
        case CompileStatus::ResourceExhausted: return "resource-exhausted";
    }
    return "unknown";
}

struct CompileResult {
    CompileStatus status{CompileStatus::CompileError};
    std::optional<std::filesystem::path> output_path;
    std::string stdout_log;
    std::string stderr_log;
    std::chrono::milliseconds compile_time{0};
    size_t memory_used_mb{0};
};

struct KernelBudget {
    std::chrono::milliseconds max_runtime{1000};
    size_t max_registers_per_thread{64};
    size_t max_shared_memory_bytes{48 * 1024};
    size_t max_threads_per_block{1024};
    size_t max_blocks{65535};
};

enum class LoadStatus : uint8_t {
    Success = 0,
    FileNotFound = 1,
    InvalidModule = 2,
    SymbolNotFound = 3,
    BudgetExceeded = 4,
    DeviceError = 5
};

constexpr const char* LoadStatusName(LoadStatus s) {
    switch (s) {
        case LoadStatus::Success:       return "success";
        case LoadStatus::FileNotFound:  return "file-not-found";
        case LoadStatus::InvalidModule: return "invalid-module";
        case LoadStatus::SymbolNotFound:return "symbol-not-found";
        case LoadStatus::BudgetExceeded:return "budget-exceeded";
        case LoadStatus::DeviceError:   return "device-error";
    }
    return "unknown";
}

struct LoadedKernel {
    std::string name;
    void* function_ptr{nullptr};
    size_t registers_per_thread{0};
    size_t shared_memory_bytes{0};
    size_t max_threads_per_block{0};
};

struct LoadResult {
    LoadStatus status{LoadStatus::FileNotFound};
    std::vector<LoadedKernel> kernels;
    std::string module_id;
    std::chrono::milliseconds load_time{0};
};

struct AuditEntry {
    std::chrono::system_clock::time_point timestamp;
    std::string session_id;
    std::string action;
    bool allowed{false};
    std::string source_hash;
    std::string details;
};

class CudaAnalyzer {
public:
    explicit CudaAnalyzer(AnalyzerConfig config = {});
    ~CudaAnalyzer();

    CudaAnalyzer(const CudaAnalyzer&) = delete;
    CudaAnalyzer& operator=(const CudaAnalyzer&) = delete;
    CudaAnalyzer(CudaAnalyzer&&) noexcept;
    CudaAnalyzer& operator=(CudaAnalyzer&&) noexcept;

    AnalysisResult Analyze(std::string_view cuda_source) const;
    AnalysisResult AnalyzeFile(const std::filesystem::path& path) const;

    const AnalyzerConfig& GetConfig() const;
    static const std::set<std::string>& AllowedSymbols();
    static const std::set<std::string>& ForbiddenPatterns();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class CudaCompiler {
public:
    explicit CudaCompiler(CompileConfig config);
    ~CudaCompiler();

    CudaCompiler(const CudaCompiler&) = delete;
    CudaCompiler& operator=(const CudaCompiler&) = delete;
    CudaCompiler(CudaCompiler&&) noexcept;
    CudaCompiler& operator=(CudaCompiler&&) noexcept;

    CompileResult Compile(std::string_view cuda_source, std::string_view output_name);
    CompileResult CompileFile(const std::filesystem::path& source_path,
                              std::string_view output_name);

    const CompileConfig& GetConfig() const;
    bool IsBubblewrapAvailable() const;
    std::string BuildBubblewrapCommand(const std::filesystem::path& source,
                                       const std::filesystem::path& output) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class CudaLoader {
public:
    explicit CudaLoader(KernelBudget budget = {});
    ~CudaLoader();

    CudaLoader(const CudaLoader&) = delete;
    CudaLoader& operator=(const CudaLoader&) = delete;
    CudaLoader(CudaLoader&&) noexcept;
    CudaLoader& operator=(CudaLoader&&) noexcept;

    LoadResult Load(const std::filesystem::path& cubin_path);
    void Unload(std::string_view module_id);
    void UnloadAll();

    std::optional<LoadedKernel> GetKernel(std::string_view module_id,
                                          std::string_view kernel_name) const;

    const KernelBudget& GetBudget() const;
    size_t LoadedModuleCount() const;
    std::vector<std::string> LoadedModuleIds() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class CudaSandbox {
public:
    struct Config {
        AnalyzerConfig analyzer;
        CompileConfig compiler;
        KernelBudget budget;
        std::filesystem::path audit_log_path;
        bool require_analysis_pass{true};
        std::string session_id;
    };

    explicit CudaSandbox(Config config);
    ~CudaSandbox();

    CudaSandbox(const CudaSandbox&) = delete;
    CudaSandbox& operator=(const CudaSandbox&) = delete;
    CudaSandbox(CudaSandbox&&) noexcept;
    CudaSandbox& operator=(CudaSandbox&&) noexcept;

    struct ExecutionResult {
        bool success{false};
        AnalysisResult analysis;
        CompileResult compile;
        LoadResult load;
        std::string error_message;
    };

    ExecutionResult CompileAndLoad(std::string_view cuda_source,
                                   std::string_view kernel_name);

    std::optional<LoadedKernel> GetKernel(std::string_view kernel_name) const;

    void WriteAuditEntry(const AuditEntry& entry);
    std::vector<AuditEntry> ReadAuditLog() const;
    void ClearAuditLog();

    const Config& GetConfig() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::string ComputeSourceHash(std::string_view source);

}  // namespace llmap::claude_agent
