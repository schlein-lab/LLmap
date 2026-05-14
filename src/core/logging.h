// LLmap — Structured logging framework.
//
// Thread-safe, zero-allocation hot path, configurable sinks.
// Designed for production use: JSON output mode for log aggregation.

#pragma once

#include <cstdint>
#include <cstdio>
#include <chrono>
#include <string>
#include <string_view>
#include <source_location>
#include <atomic>
#include <mutex>
#include <functional>
#include <vector>

namespace llmap {

enum class LogLevel : uint8_t {
    kTrace = 0,
    kDebug = 1,
    kInfo  = 2,
    kWarn  = 3,
    kError = 4,
    kFatal = 5,
    kOff   = 6
};

constexpr std::string_view LogLevelName(LogLevel level) noexcept {
    switch (level) {
        case LogLevel::kTrace: return "TRACE";
        case LogLevel::kDebug: return "DEBUG";
        case LogLevel::kInfo:  return "INFO";
        case LogLevel::kWarn:  return "WARN";
        case LogLevel::kError: return "ERROR";
        case LogLevel::kFatal: return "FATAL";
        case LogLevel::kOff:   return "OFF";
    }
    return "UNKNOWN";
}

constexpr std::string_view LogLevelNameLower(LogLevel level) noexcept {
    switch (level) {
        case LogLevel::kTrace: return "trace";
        case LogLevel::kDebug: return "debug";
        case LogLevel::kInfo:  return "info";
        case LogLevel::kWarn:  return "warn";
        case LogLevel::kError: return "error";
        case LogLevel::kFatal: return "fatal";
        case LogLevel::kOff:   return "off";
    }
    return "unknown";
}

struct LogRecord {
    LogLevel level;
    std::string_view message;
    std::string_view file;
    std::string_view function;
    uint32_t line;
    std::chrono::system_clock::time_point timestamp;
    uint64_t thread_id;
};

enum class LogFormat : uint8_t {
    kText,
    kJson
};

using LogSink = std::function<void(const LogRecord&)>;

class Logger {
public:
    static Logger& Instance() noexcept;

    void SetLevel(LogLevel level) noexcept {
        level_.store(level, std::memory_order_release);
    }

    LogLevel GetLevel() const noexcept {
        return level_.load(std::memory_order_acquire);
    }

    void SetFormat(LogFormat format) noexcept {
        format_.store(format, std::memory_order_release);
    }

    LogFormat GetFormat() const noexcept {
        return format_.load(std::memory_order_acquire);
    }

    void SetOutput(FILE* out) noexcept {
        std::lock_guard<std::mutex> lock(mutex_);
        output_ = out;
    }

    void AddSink(LogSink sink) {
        std::lock_guard<std::mutex> lock(mutex_);
        sinks_.push_back(std::move(sink));
    }

    void ClearSinks() {
        std::lock_guard<std::mutex> lock(mutex_);
        sinks_.clear();
    }

    bool ShouldLog(LogLevel level) const noexcept {
        return level >= level_.load(std::memory_order_acquire);
    }

    void Log(LogLevel level, std::string_view message,
             const std::source_location& loc = std::source_location::current());

    void LogWithContext(LogLevel level, std::string_view message,
                        std::string_view context,
                        const std::source_location& loc = std::source_location::current());

    void Flush();

private:
    Logger();
    ~Logger() = default;

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void WriteText(const LogRecord& record);
    void WriteJson(const LogRecord& record);
    uint64_t GetThreadId() const noexcept;

    std::atomic<LogLevel> level_{LogLevel::kInfo};
    std::atomic<LogFormat> format_{LogFormat::kText};
    FILE* output_{stderr};
    std::mutex mutex_;
    std::vector<LogSink> sinks_;
    char buffer_[2048];
};

// Convenience macros with compile-time level check
#define LLMAP_LOG(level, msg) \
    do { \
        if (::llmap::Logger::Instance().ShouldLog(level)) { \
            ::llmap::Logger::Instance().Log(level, msg); \
        } \
    } while (false)

#define LLMAP_LOG_TRACE(msg) LLMAP_LOG(::llmap::LogLevel::kTrace, msg)
#define LLMAP_LOG_DEBUG(msg) LLMAP_LOG(::llmap::LogLevel::kDebug, msg)
#define LLMAP_LOG_INFO(msg)  LLMAP_LOG(::llmap::LogLevel::kInfo, msg)
#define LLMAP_LOG_WARN(msg)  LLMAP_LOG(::llmap::LogLevel::kWarn, msg)
#define LLMAP_LOG_ERROR(msg) LLMAP_LOG(::llmap::LogLevel::kError, msg)
#define LLMAP_LOG_FATAL(msg) LLMAP_LOG(::llmap::LogLevel::kFatal, msg)

// Context-aware logging for structured diagnostics
#define LLMAP_LOG_CTX(level, msg, ctx) \
    do { \
        if (::llmap::Logger::Instance().ShouldLog(level)) { \
            ::llmap::Logger::Instance().LogWithContext(level, msg, ctx); \
        } \
    } while (false)

// Parse log level from string (case-insensitive)
LogLevel ParseLogLevel(std::string_view str) noexcept;

// Initialize logger from environment variables:
//   LLMAP_LOG_LEVEL: trace|debug|info|warn|error|fatal|off
//   LLMAP_LOG_FORMAT: text|json
void InitLoggerFromEnv();

}  // namespace llmap
