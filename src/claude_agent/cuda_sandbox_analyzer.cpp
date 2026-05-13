#include "cuda_sandbox.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <regex>
#include <sstream>

namespace llmap::claude_agent {

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

struct LineInfo {
    size_t line_number;
    std::string content;
    std::string trimmed;
};

std::vector<LineInfo> ParseLines(std::string_view source) {
    std::vector<LineInfo> lines;
    std::istringstream stream{std::string(source)};
    std::string line;
    size_t line_num = 0;

    while (std::getline(stream, line)) {
        ++line_num;
        std::string trimmed = line;
        auto start = trimmed.find_first_not_of(" \t");
        if (start != std::string::npos) {
            trimmed = trimmed.substr(start);
        }
        auto end = trimmed.find_last_not_of(" \t\r\n");
        if (end != std::string::npos) {
            trimmed = trimmed.substr(0, end + 1);
        }
        lines.push_back({line_num, std::move(line), std::move(trimmed)});
    }
    return lines;
}

bool IsCommentOrPreprocessor(const std::string& trimmed) {
    return trimmed.empty() ||
           trimmed.starts_with("//") ||
           trimmed.starts_with("/*") ||
           trimmed.starts_with("*") ||
           trimmed.starts_with("#");
}

std::string StripAllComments(std::string_view source) {
    std::string result;
    result.reserve(source.size());
    bool in_line_comment = false;
    bool in_block_comment = false;
    bool in_string = false;
    char string_char = 0;

    for (size_t i = 0; i < source.size(); ++i) {
        char c = source[i];
        char next = (i + 1 < source.size()) ? source[i + 1] : 0;

        if (in_block_comment) {
            if (c == '\n') {
                result += '\n';
            } else {
                result += ' ';
            }
            if (c == '*' && next == '/') {
                in_block_comment = false;
                ++i;
                result += ' ';
            }
            continue;
        }

        if (in_line_comment) {
            if (c == '\n') {
                in_line_comment = false;
                result += '\n';
            } else {
                result += ' ';
            }
            continue;
        }

        if (in_string) {
            result += c;
            if (c == '\\' && i + 1 < source.size()) {
                ++i;
                result += source[i];
                continue;
            }
            if (c == string_char) {
                in_string = false;
            }
            continue;
        }

        if (c == '/' && next == '/') {
            in_line_comment = true;
            result += ' ';
            ++i;
            result += ' ';
            continue;
        }

        if (c == '/' && next == '*') {
            in_block_comment = true;
            result += ' ';
            ++i;
            result += ' ';
            continue;
        }

        if (c == '"' || c == '\'') {
            in_string = true;
            string_char = c;
        }

        result += c;
    }

    return result;
}

std::string StripCommentsAndStrings(const std::string& content) {
    std::string result;
    result.reserve(content.size());
    bool in_line_comment = false;
    bool in_string = false;
    char string_char = 0;

    for (size_t i = 0; i < content.size(); ++i) {
        char c = content[i];

        if (in_line_comment) {
            result += ' ';
            continue;
        }

        if (in_string) {
            result += ' ';
            if (c == '\\' && i + 1 < content.size()) {
                ++i;
                result += ' ';
                continue;
            }
            if (c == string_char) {
                in_string = false;
            }
            continue;
        }

        if (c == '/' && i + 1 < content.size() && content[i + 1] == '/') {
            in_line_comment = true;
            result += ' ';
            continue;
        }

        if (c == '"' || c == '\'') {
            in_string = true;
            string_char = c;
            result += ' ';
            continue;
        }

        result += c;
    }

    return result;
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

}  // namespace

struct CudaAnalyzer::Impl {
    AnalyzerConfig config;

    explicit Impl(AnalyzerConfig cfg) : config(std::move(cfg)) {}
};

CudaAnalyzer::CudaAnalyzer(AnalyzerConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

CudaAnalyzer::~CudaAnalyzer() = default;

CudaAnalyzer::CudaAnalyzer(CudaAnalyzer&&) noexcept = default;
CudaAnalyzer& CudaAnalyzer::operator=(CudaAnalyzer&&) noexcept = default;

AnalysisResult CudaAnalyzer::Analyze(std::string_view cuda_source) const {
    auto start = std::chrono::steady_clock::now();

    AnalysisResult result;
    result.safe = true;

    auto original_lines = ParseLines(cuda_source);
    result.lines_analyzed = original_lines.size();

    if (original_lines.size() > impl_->config.max_lines) {
        result.safe = false;
        result.violations.push_back({
            .type = SecurityViolation::ExcessiveComplexity,
            .line = 0,
            .column = 0,
            .description = "Source exceeds maximum line count (" +
                          std::to_string(impl_->config.max_lines) + ")",
            .context = ""
        });
    }

    std::string stripped_source = StripAllComments(cuda_source);
    auto lines = ParseLines(stripped_source);

    for (size_t i = 0; i < lines.size(); ++i) {
        const auto& line = lines[i];
        if (IsCommentOrPreprocessor(line.trimmed)) continue;

        for (const auto& pattern : kForbiddenPatterns) {
            if (auto violation = CheckForbiddenPattern(line, pattern)) {
                violation->line = line.line_number;
                if (i < original_lines.size()) {
                    violation->context = original_lines[i].trimmed.substr(0, 80);
                }
                result.safe = false;
                result.violations.push_back(std::move(*violation));
            }
        }

        if (auto violation = CheckInlineAssembly(line)) {
            violation->line = line.line_number;
            if (i < original_lines.size()) {
                violation->context = original_lines[i].trimmed.substr(0, 80);
            }
            result.safe = false;
            result.violations.push_back(std::move(*violation));
        }
    }

    result.functions_found = CountFunctions(lines);
    if (result.functions_found > impl_->config.max_functions) {
        result.safe = false;
        result.violations.push_back({
            .type = SecurityViolation::ExcessiveComplexity,
            .line = 0,
            .column = 0,
            .description = "Too many functions (" + std::to_string(result.functions_found) +
                          " > " + std::to_string(impl_->config.max_functions) + ")",
            .context = ""
        });
    }

    result.kernel_count = CountKernels(lines);
    if (result.kernel_count > impl_->config.max_kernels) {
        result.safe = false;
        result.violations.push_back({
            .type = SecurityViolation::ExcessiveComplexity,
            .line = 0,
            .column = 0,
            .description = "Too many kernels (" + std::to_string(result.kernel_count) +
                          " > " + std::to_string(impl_->config.max_kernels) + ")",
            .context = ""
        });
    }

    result.symbols_used = ExtractSymbols(lines);

    auto end = std::chrono::steady_clock::now();
    result.analysis_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    return result;
}

AnalysisResult CudaAnalyzer::AnalyzeFile(const std::filesystem::path& path) const {
    std::ifstream file(path);
    if (!file) {
        return {
            .safe = false,
            .violations = {{
                .type = SecurityViolation::MalformedCode,
                .description = "Failed to open file: " + path.string()
            }}
        };
    }

    std::ostringstream ss;
    ss << file.rdbuf();
    return Analyze(ss.str());
}

const AnalyzerConfig& CudaAnalyzer::GetConfig() const {
    return impl_->config;
}

const std::set<std::string>& CudaAnalyzer::AllowedSymbols() {
    return kAllowedSymbols;
}

const std::set<std::string>& CudaAnalyzer::ForbiddenPatterns() {
    return kForbiddenPatterns;
}

std::string ComputeSourceHash(std::string_view source) {
    uint64_t hash = 0xcbf29ce484222325ULL;
    constexpr uint64_t prime = 0x100000001b3ULL;
    for (char c : source) {
        hash ^= static_cast<uint8_t>(c);
        hash *= prime;
    }
    char buf[17];
    snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(hash));
    return std::string(buf);
}

}  // namespace llmap::claude_agent
