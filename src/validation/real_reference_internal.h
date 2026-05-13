// LLmap — Internal utilities for real_reference module.
// Not part of public API.

#pragma once

#include <chrono>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::validation::internal {

class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    float ElapsedMs() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<float, std::milli>(now - start_).count();
    }

private:
    std::chrono::steady_clock::time_point start_;
};

std::vector<std::string_view> SplitLine(std::string_view line, char delim);

uint64_t ParseUint64(std::string_view s);

std::string TrimWhitespace(std::string_view s);

}  // namespace llmap::validation::internal
