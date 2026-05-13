#include "cuda_sandbox.h"
#include "cuda_sandbox_internal.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace llmap::claude_agent {

using namespace internal;

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

    const auto& forbidden = ForbiddenPatternSet();
    for (size_t i = 0; i < lines.size(); ++i) {
        const auto& line = lines[i];
        if (IsCommentOrPreprocessor(line.trimmed)) continue;

        for (const auto& pattern : forbidden) {
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
    return AllowedSymbolSet();
}

const std::set<std::string>& CudaAnalyzer::ForbiddenPatterns() {
    return ForbiddenPatternSet();
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
