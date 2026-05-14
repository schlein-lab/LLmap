// LLmap — ClassicalPipeline internal header for split implementation.
#pragma once

#include "classical/classical_pipeline.h"

#include <chrono>

namespace llmap::classical {

// Helper to measure time (shared across implementation files)
class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    float ElapsedUs() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<float, std::micro>(now - start_).count();
    }

    float ElapsedMs() const {
        return ElapsedUs() / 1000.0f;
    }

private:
    std::chrono::steady_clock::time_point start_;
};

}  // namespace llmap::classical
