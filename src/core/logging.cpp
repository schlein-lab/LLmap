// LLmap — Structured logging implementation.

#include "core/logging.h"

#include <algorithm>
#include <cstring>
#include <ctime>
#include <thread>

#ifdef __linux__
#include <sys/syscall.h>
#include <unistd.h>
#endif

namespace llmap {

Logger& Logger::Instance() noexcept {
    static Logger instance;
    return instance;
}

Logger::Logger() {
    InitLoggerFromEnv();
}

uint64_t Logger::GetThreadId() const noexcept {
#ifdef __linux__
    return static_cast<uint64_t>(syscall(SYS_gettid));
#else
    return std::hash<std::thread::id>{}(std::this_thread::get_id());
#endif
}

void Logger::Log(LogLevel level, std::string_view message,
                 const std::source_location& loc) {
    LogRecord record{
        .level = level,
        .message = message,
        .file = loc.file_name(),
        .function = loc.function_name(),
        .line = loc.line(),
        .timestamp = std::chrono::system_clock::now(),
        .thread_id = GetThreadId()
    };

    std::lock_guard<std::mutex> lock(mutex_);

    if (format_.load(std::memory_order_acquire) == LogFormat::kJson) {
        WriteJson(record);
    } else {
        WriteText(record);
    }

    for (const auto& sink : sinks_) {
        sink(record);
    }
}

void Logger::LogWithContext(LogLevel level, std::string_view message,
                            std::string_view context,
                            const std::source_location& loc) {
    std::string full_message;
    full_message.reserve(message.size() + context.size() + 16);
    full_message.append(message);
    full_message.append(" [");
    full_message.append(context);
    full_message.append("]");

    Log(level, full_message, loc);
}

void Logger::WriteText(const LogRecord& record) {
    auto time = std::chrono::system_clock::to_time_t(record.timestamp);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        record.timestamp.time_since_epoch()) % 1000;

    std::tm tm_buf{};
    localtime_r(&time, &tm_buf);

    // Extract filename from path
    const char* filename = record.file.data();
    const char* last_slash = std::strrchr(filename, '/');
    if (last_slash) {
        filename = last_slash + 1;
    }

    int len = std::snprintf(
        buffer_, sizeof(buffer_),
        "%04d-%02d-%02d %02d:%02d:%02d.%03d [%s] [%lu] %s:%u: %.*s\n",
        tm_buf.tm_year + 1900, tm_buf.tm_mon + 1, tm_buf.tm_mday,
        tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec,
        static_cast<int>(ms.count()),
        LogLevelName(record.level).data(),
        record.thread_id,
        filename, record.line,
        static_cast<int>(record.message.size()), record.message.data()
    );

    if (len > 0 && output_) {
        std::fwrite(buffer_, 1, static_cast<size_t>(len), output_);
    }
}

void Logger::WriteJson(const LogRecord& record) {
    auto time = std::chrono::system_clock::to_time_t(record.timestamp);
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
        record.timestamp.time_since_epoch()) % 1000000;

    std::tm tm_buf{};
    gmtime_r(&time, &tm_buf);

    // Escape JSON string characters in message
    std::string escaped_message;
    escaped_message.reserve(record.message.size() + 16);
    for (char c : record.message) {
        switch (c) {
            case '"':  escaped_message += "\\\""; break;
            case '\\': escaped_message += "\\\\"; break;
            case '\n': escaped_message += "\\n"; break;
            case '\r': escaped_message += "\\r"; break;
            case '\t': escaped_message += "\\t"; break;
            default:   escaped_message += c; break;
        }
    }

    // Extract filename from path
    const char* filename = record.file.data();
    const char* last_slash = std::strrchr(filename, '/');
    if (last_slash) {
        filename = last_slash + 1;
    }

    int len = std::snprintf(
        buffer_, sizeof(buffer_),
        "{\"ts\":\"%04d-%02d-%02dT%02d:%02d:%02d.%06ldZ\","
        "\"level\":\"%.*s\","
        "\"tid\":%lu,"
        "\"file\":\"%s\","
        "\"line\":%u,"
        "\"msg\":\"%s\"}\n",
        tm_buf.tm_year + 1900, tm_buf.tm_mon + 1, tm_buf.tm_mday,
        tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec,
        static_cast<long>(us.count()),
        static_cast<int>(LogLevelNameLower(record.level).size()),
        LogLevelNameLower(record.level).data(),
        record.thread_id,
        filename, record.line,
        escaped_message.c_str()
    );

    if (len > 0 && output_) {
        std::fwrite(buffer_, 1, static_cast<size_t>(len), output_);
    }
}

void Logger::Flush() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (output_) {
        std::fflush(output_);
    }
}

LogLevel ParseLogLevel(std::string_view str) noexcept {
    std::string lower;
    lower.reserve(str.size());
    for (char c : str) {
        lower += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (lower == "trace") return LogLevel::kTrace;
    if (lower == "debug") return LogLevel::kDebug;
    if (lower == "info")  return LogLevel::kInfo;
    if (lower == "warn" || lower == "warning") return LogLevel::kWarn;
    if (lower == "error") return LogLevel::kError;
    if (lower == "fatal") return LogLevel::kFatal;
    if (lower == "off")   return LogLevel::kOff;

    return LogLevel::kInfo;  // Default
}

void InitLoggerFromEnv() {
    if (const char* level_env = std::getenv("LLMAP_LOG_LEVEL")) {
        Logger::Instance().SetLevel(ParseLogLevel(level_env));
    }

    if (const char* format_env = std::getenv("LLMAP_LOG_FORMAT")) {
        std::string_view fmt{format_env};
        if (fmt == "json" || fmt == "JSON") {
            Logger::Instance().SetFormat(LogFormat::kJson);
        } else {
            Logger::Instance().SetFormat(LogFormat::kText);
        }
    }
}

}  // namespace llmap
