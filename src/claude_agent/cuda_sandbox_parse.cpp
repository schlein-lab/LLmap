#include "cuda_sandbox_internal.h"

#include <sstream>

namespace llmap::claude_agent::internal {

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

}  // namespace llmap::claude_agent::internal
