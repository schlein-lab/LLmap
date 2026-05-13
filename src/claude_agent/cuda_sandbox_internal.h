#pragma once

#include "cuda_sandbox.h"

#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::claude_agent::internal {

struct LineInfo {
    size_t line_number;
    std::string content;
    std::string trimmed;
};

std::vector<LineInfo> ParseLines(std::string_view source);

bool IsCommentOrPreprocessor(const std::string& trimmed);

std::string StripAllComments(std::string_view source);

std::string StripCommentsAndStrings(const std::string& content);

const std::set<std::string>& AllowedSymbolSet();
const std::set<std::string>& ForbiddenPatternSet();

std::optional<AnalysisViolation> CheckForbiddenPattern(const LineInfo& line,
                                                       const std::string& pattern);

std::optional<AnalysisViolation> CheckInlineAssembly(const LineInfo& line);

bool ContainsKernelDeclaration(const std::string& content);

bool ContainsDeviceFunction(const std::string& content);

size_t CountFunctions(const std::vector<LineInfo>& lines);

size_t CountKernels(const std::vector<LineInfo>& lines);

std::set<std::string> ExtractSymbols(const std::vector<LineInfo>& lines);

}  // namespace llmap::claude_agent::internal
