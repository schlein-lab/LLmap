// LLmap — Version utility functions.

#pragma once

#include <string>
#include <string_view>
#include "core/version.h"

namespace llmap {

// Format version as "MAJOR.MINOR.PATCH"
inline std::string FormatVersion() {
    return std::string(kVersion);
}

// Format short version line: "llmap 0.1.0"
inline std::string FormatVersionShort() {
    return "llmap " + std::string(kVersion);
}

// Format full version info for --version output
inline std::string FormatVersionFull() {
    std::string result;
    result.reserve(256);
    result += "llmap ";
    result += kVersion;
    result += "\n  commit:   ";
    result += kGitCommit;
    result += "\n  built:    ";
    result += kBuildDate;
    result += "\n  type:     ";
    result += kBuildType;
    result += "\n  compiler: ";
    result += kCompilerId;
    result += " ";
    result += kCompilerVersion;
    result += "\n  features: ";
    result += kFeatures;
    result += "\n";
    return result;
}

// Check if a specific feature is enabled
inline bool HasCuda() { return kHasCuda; }
inline bool HasOnnxRuntime() { return kHasOnnxRuntime; }
inline bool HasFaiss() { return kHasFaiss; }
inline bool HasFaissGpu() { return kHasFaissGpu; }
inline bool HasClaude() { return kHasClaude; }

}  // namespace llmap
