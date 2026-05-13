#include "cuda_sandbox_internal.h"

#include <cctype>
#include <regex>

namespace llmap::claude_agent::internal {

namespace {

const std::set<std::string> kAllowedSymbols = {
    // CUDA intrinsics
    "__syncthreads", "__syncwarp", "__threadfence", "__threadfence_block",
    "__ballot_sync", "__shfl_sync", "__shfl_up_sync", "__shfl_down_sync", "__shfl_xor_sync",
    "atomicAdd", "atomicSub", "atomicMin", "atomicMax", "atomicExch", "atomicCAS",
    "atomicAnd", "atomicOr", "atomicXor", "atomicInc", "atomicDec",

    // CUDA memory
    "__shared__", "__constant__", "__device__", "__global__", "__host__",
    "__restrict__", "__align__",

    // Math functions
    "sqrtf", "sqrt", "rsqrtf", "rsqrt", "cbrtf", "cbrt",
    "expf", "exp", "exp2f", "exp2", "exp10f", "exp10", "expm1f", "expm1",
    "logf", "log", "log2f", "log2", "log10f", "log10", "log1pf", "log1p",
    "powf", "pow", "sinf", "sin", "cosf", "cos", "tanf", "tan",
    "asinf", "asin", "acosf", "acos", "atanf", "atan", "atan2f", "atan2",
    "sinhf", "sinh", "coshf", "cosh", "tanhf", "tanh",
    "asinhf", "asinh", "acoshf", "acosh", "atanhf", "atanh",
    "fabsf", "fabs", "fmaxf", "fmax", "fminf", "fmin", "fmodf", "fmod",
    "ceilf", "ceil", "floorf", "floor", "truncf", "trunc", "roundf", "round",
    "copysignf", "copysign", "fdimf", "fdim",
    "__float2half", "__half2float", "__float2half_rn",
    "__ldg",

    // Thread indexing
    "threadIdx", "blockIdx", "blockDim", "gridDim", "warpSize",

    // Cooperative groups
    "cooperative_groups",

    // Type casts
    "__float_as_int", "__int_as_float", "__double_as_longlong", "__longlong_as_double",
};

const std::set<std::string> kForbiddenPatterns = {
    // System calls
    "syscall", "__syscall", "__NR_", "SYS_",

    // File I/O
    "fopen", "fclose", "fread", "fwrite", "fseek", "ftell", "fflush", "fgets", "fputs",
    "fprintf", "fscanf", "fgetc", "fputc", "open", "close", "read", "write",
    "lseek", "access", "unlink", "remove", "rename", "mkdir", "rmdir",
    "opendir", "readdir", "closedir", "stat", "fstat", "lstat", "chmod", "chown",
    "FILE", "stdin", "stdout", "stderr",

    // Network
    "socket", "connect", "bind", "listen", "accept", "send", "recv", "sendto", "recvfrom",
    "gethostbyname", "getaddrinfo", "inet_", "htons", "htonl", "ntohs", "ntohl",
    "AF_INET", "SOCK_STREAM", "SOCK_DGRAM",

    // Process control
    "fork", "exec", "execv", "execve", "execvp", "system", "popen", "pclose",
    "exit", "_exit", "abort", "kill", "signal", "sigaction", "raise",
    "setjmp", "longjmp",

    // Memory mapping
    "mmap", "munmap", "mprotect", "mremap", "madvise",

    // Dynamic linking
    "dlopen", "dlclose", "dlsym", "dlerror",

    // Shell/environment
    "getenv", "setenv", "putenv", "unsetenv", "environ",

    // Dangerous string ops
    "gets", "sprintf", "vsprintf",

    // CUDA driver API (only runtime API allowed)
    "cuModuleLoad", "cuModuleLoadData", "cuModuleGetFunction",
    "cuLaunchKernel", "cuCtxCreate", "cuCtxDestroy",
};

}  // namespace

const std::set<std::string>& AllowedSymbolSet() {
    return kAllowedSymbols;
}

const std::set<std::string>& ForbiddenPatternSet() {
    return kForbiddenPatterns;
}

std::optional<AnalysisViolation> CheckForbiddenPattern(const LineInfo& line,
                                                       const std::string& pattern) {
    std::string stripped = StripCommentsAndStrings(line.content);
    auto pos = stripped.find(pattern);
    if (pos == std::string::npos) return std::nullopt;

    if (pos > 0) {
        char prev = stripped[pos - 1];
        if (std::isalnum(prev) || prev == '_') return std::nullopt;
    }
    if (pos + pattern.size() < stripped.size()) {
        char next = stripped[pos + pattern.size()];
        if (std::isalnum(next) || next == '_') return std::nullopt;
    }

    SecurityViolation type = SecurityViolation::ForbiddenSymbol;
    if (pattern.find("socket") != std::string::npos ||
        pattern.find("connect") != std::string::npos ||
        pattern.find("AF_INET") != std::string::npos) {
        type = SecurityViolation::NetworkDetected;
    } else if (pattern.find("fopen") != std::string::npos ||
               pattern.find("fread") != std::string::npos ||
               pattern.find("FILE") != std::string::npos ||
               pattern == "open" || pattern == "read" || pattern == "write") {
        type = SecurityViolation::FileIODetected;
    } else if (pattern.find("syscall") != std::string::npos ||
               pattern.find("fork") != std::string::npos ||
               pattern.find("exec") != std::string::npos ||
               pattern.find("system") != std::string::npos) {
        type = SecurityViolation::SyscallDetected;
    }

    return AnalysisViolation{
        .type = type,
        .line = line.line_number,
        .column = pos + 1,
        .description = "Forbidden pattern detected: " + pattern,
        .context = line.trimmed.substr(0, 80)
    };
}

std::optional<AnalysisViolation> CheckInlineAssembly(const LineInfo& line) {
    std::string stripped = StripCommentsAndStrings(line.content);
    static const std::regex asm_pattern(R"(\b(asm|__asm__|__asm)\s*(\(|volatile))");
    std::smatch match;
    if (std::regex_search(stripped, match, asm_pattern)) {
        return AnalysisViolation{
            .type = SecurityViolation::AssemblyBlock,
            .line = line.line_number,
            .column = static_cast<size_t>(match.position() + 1),
            .description = "Inline assembly is not allowed in agent-generated kernels",
            .context = line.trimmed.substr(0, 80)
        };
    }
    return std::nullopt;
}

bool ContainsKernelDeclaration(const std::string& content) {
    static const std::regex kernel_pattern(R"(__global__\s+\w)");
    return std::regex_search(content, kernel_pattern);
}

bool ContainsDeviceFunction(const std::string& content) {
    static const std::regex device_pattern(R"(__device__\s+\w)");
    return std::regex_search(content, device_pattern);
}

size_t CountFunctions(const std::vector<LineInfo>& lines) {
    static const std::regex func_pattern(R"(\b(\w+)\s+(\w+)\s*\([^)]*\)\s*\{?)");
    size_t count = 0;
    for (const auto& line : lines) {
        if (std::regex_search(line.content, func_pattern)) {
            ++count;
        }
    }
    return count;
}

size_t CountKernels(const std::vector<LineInfo>& lines) {
    static const std::regex kernel_pattern(R"(__global__\s+\w+\s+\w+\s*\()");
    size_t count = 0;
    for (const auto& line : lines) {
        if (std::regex_search(line.content, kernel_pattern)) {
            ++count;
        }
    }
    return count;
}

std::set<std::string> ExtractSymbols(const std::vector<LineInfo>& lines) {
    std::set<std::string> symbols;
    static const std::regex symbol_pattern(R"(\b([a-zA-Z_]\w*)\s*\()");

    for (const auto& line : lines) {
        if (IsCommentOrPreprocessor(line.trimmed)) continue;

        std::sregex_iterator iter(line.content.begin(), line.content.end(), symbol_pattern);
        std::sregex_iterator end;

        while (iter != end) {
            symbols.insert((*iter)[1].str());
            ++iter;
        }
    }
    return symbols;
}

}  // namespace llmap::claude_agent::internal
