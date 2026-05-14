// LLmap — Error handling framework.
//
// Result<T, E> type for explicit error handling without exceptions.
// LLmapError hierarchy for typed error categorization.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <variant>
#include <optional>
#include <utility>
#include <source_location>
#include <type_traits>
#include <vector>

namespace llmap {

enum class ErrorCode : uint16_t {
    kOk = 0,

    // I/O errors (100-199)
    kIoFileNotFound = 100,
    kIoPermissionDenied = 101,
    kIoReadError = 102,
    kIoWriteError = 103,
    kIoSeekError = 104,
    kIoEndOfFile = 105,
    kIoInvalidPath = 106,

    // Parse errors (200-299)
    kParseInvalidFormat = 200,
    kParseUnexpectedToken = 201,
    kParseMissingField = 202,
    kParseInvalidValue = 203,
    kParseCorruptedData = 204,
    kParseVersionMismatch = 205,

    // Config errors (300-399)
    kConfigMissingRequired = 300,
    kConfigInvalidValue = 301,
    kConfigConflict = 302,
    kConfigUnsupported = 303,

    // Validation errors (400-499)
    kValidateEmpty = 400,
    kValidateOutOfRange = 401,
    kValidateInvariantViolation = 402,
    kValidateChecksum = 403,
    kValidateDuplicate = 404,

    // Resource errors (500-599)
    kResourceOutOfMemory = 500,
    kResourceTimeout = 501,
    kResourceExhausted = 502,
    kResourceUnavailable = 503,
    kResourceBusy = 504,

    // System errors (600-699)
    kSystemOsError = 600,
    kSystemCudaError = 601,
    kSystemThreadError = 602,
    kSystemSignal = 603,

    // Algorithm errors (700-799)
    kAlgoConvergenceFailed = 700,
    kAlgoNumericalInstability = 701,
    kAlgoInvalidState = 702,
    kAlgoNotImplemented = 703,

    // External errors (800-899)
    kExternalApiError = 800,
    kExternalNetworkError = 801,
    kExternalDependencyMissing = 802,

    // Unknown
    kUnknown = 999
};

constexpr std::string_view ErrorCodeName(ErrorCode code) noexcept {
    switch (code) {
        case ErrorCode::kOk: return "OK";
        case ErrorCode::kIoFileNotFound: return "IO_FILE_NOT_FOUND";
        case ErrorCode::kIoPermissionDenied: return "IO_PERMISSION_DENIED";
        case ErrorCode::kIoReadError: return "IO_READ_ERROR";
        case ErrorCode::kIoWriteError: return "IO_WRITE_ERROR";
        case ErrorCode::kIoSeekError: return "IO_SEEK_ERROR";
        case ErrorCode::kIoEndOfFile: return "IO_END_OF_FILE";
        case ErrorCode::kIoInvalidPath: return "IO_INVALID_PATH";
        case ErrorCode::kParseInvalidFormat: return "PARSE_INVALID_FORMAT";
        case ErrorCode::kParseUnexpectedToken: return "PARSE_UNEXPECTED_TOKEN";
        case ErrorCode::kParseMissingField: return "PARSE_MISSING_FIELD";
        case ErrorCode::kParseInvalidValue: return "PARSE_INVALID_VALUE";
        case ErrorCode::kParseCorruptedData: return "PARSE_CORRUPTED_DATA";
        case ErrorCode::kParseVersionMismatch: return "PARSE_VERSION_MISMATCH";
        case ErrorCode::kConfigMissingRequired: return "CONFIG_MISSING_REQUIRED";
        case ErrorCode::kConfigInvalidValue: return "CONFIG_INVALID_VALUE";
        case ErrorCode::kConfigConflict: return "CONFIG_CONFLICT";
        case ErrorCode::kConfigUnsupported: return "CONFIG_UNSUPPORTED";
        case ErrorCode::kValidateEmpty: return "VALIDATE_EMPTY";
        case ErrorCode::kValidateOutOfRange: return "VALIDATE_OUT_OF_RANGE";
        case ErrorCode::kValidateInvariantViolation: return "VALIDATE_INVARIANT_VIOLATION";
        case ErrorCode::kValidateChecksum: return "VALIDATE_CHECKSUM";
        case ErrorCode::kValidateDuplicate: return "VALIDATE_DUPLICATE";
        case ErrorCode::kResourceOutOfMemory: return "RESOURCE_OUT_OF_MEMORY";
        case ErrorCode::kResourceTimeout: return "RESOURCE_TIMEOUT";
        case ErrorCode::kResourceExhausted: return "RESOURCE_EXHAUSTED";
        case ErrorCode::kResourceUnavailable: return "RESOURCE_UNAVAILABLE";
        case ErrorCode::kResourceBusy: return "RESOURCE_BUSY";
        case ErrorCode::kSystemOsError: return "SYSTEM_OS_ERROR";
        case ErrorCode::kSystemCudaError: return "SYSTEM_CUDA_ERROR";
        case ErrorCode::kSystemThreadError: return "SYSTEM_THREAD_ERROR";
        case ErrorCode::kSystemSignal: return "SYSTEM_SIGNAL";
        case ErrorCode::kAlgoConvergenceFailed: return "ALGO_CONVERGENCE_FAILED";
        case ErrorCode::kAlgoNumericalInstability: return "ALGO_NUMERICAL_INSTABILITY";
        case ErrorCode::kAlgoInvalidState: return "ALGO_INVALID_STATE";
        case ErrorCode::kAlgoNotImplemented: return "ALGO_NOT_IMPLEMENTED";
        case ErrorCode::kExternalApiError: return "EXTERNAL_API_ERROR";
        case ErrorCode::kExternalNetworkError: return "EXTERNAL_NETWORK_ERROR";
        case ErrorCode::kExternalDependencyMissing: return "EXTERNAL_DEPENDENCY_MISSING";
        case ErrorCode::kUnknown: return "UNKNOWN";
    }
    return "UNKNOWN";
}

constexpr bool IsIoError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 100 && static_cast<uint16_t>(code) < 200;
}

constexpr bool IsParseError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 200 && static_cast<uint16_t>(code) < 300;
}

constexpr bool IsConfigError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 300 && static_cast<uint16_t>(code) < 400;
}

constexpr bool IsValidateError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 400 && static_cast<uint16_t>(code) < 500;
}

constexpr bool IsResourceError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 500 && static_cast<uint16_t>(code) < 600;
}

constexpr bool IsSystemError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 600 && static_cast<uint16_t>(code) < 700;
}

constexpr bool IsAlgoError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 700 && static_cast<uint16_t>(code) < 800;
}

constexpr bool IsExternalError(ErrorCode code) noexcept {
    return static_cast<uint16_t>(code) >= 800 && static_cast<uint16_t>(code) < 900;
}

struct ErrorLocation {
    std::string_view file{};
    std::string_view function{};
    uint32_t line{0};

    static ErrorLocation Current(
        const std::source_location& loc = std::source_location::current()) noexcept {
        return {loc.file_name(), loc.function_name(), loc.line()};
    }
};

class LLmapError {
public:
    LLmapError() noexcept = default;

    explicit LLmapError(ErrorCode code,
                        const std::source_location& loc = std::source_location::current())
        : code_(code), location_(ErrorLocation::Current(loc)) {}

    LLmapError(ErrorCode code, std::string message,
               const std::source_location& loc = std::source_location::current())
        : code_(code), message_(std::move(message)), location_(ErrorLocation::Current(loc)) {}

    LLmapError(ErrorCode code, std::string message, std::string context,
               const std::source_location& loc = std::source_location::current())
        : code_(code), message_(std::move(message)), context_(std::move(context)),
          location_(ErrorLocation::Current(loc)) {}

    ErrorCode code() const noexcept { return code_; }
    std::string_view message() const noexcept { return message_; }
    std::string_view context() const noexcept { return context_; }
    const ErrorLocation& location() const noexcept { return location_; }

    bool ok() const noexcept { return code_ == ErrorCode::kOk; }
    explicit operator bool() const noexcept { return !ok(); }

    LLmapError& WithContext(std::string ctx) {
        context_ = std::move(ctx);
        return *this;
    }

    LLmapError& WithOsError(int errno_val);

    std::string ToString() const;
    std::string ToJson() const;

    bool operator==(const LLmapError& other) const noexcept {
        return code_ == other.code_;
    }

    bool operator!=(const LLmapError& other) const noexcept {
        return !(*this == other);
    }

private:
    ErrorCode code_{ErrorCode::kOk};
    std::string message_;
    std::string context_;
    ErrorLocation location_;
    int os_errno_{0};
};

// Factory functions for common errors
LLmapError IoError(ErrorCode code, std::string_view path,
                   const std::source_location& loc = std::source_location::current());
LLmapError ParseError(ErrorCode code, std::string_view what, size_t line = 0,
                      const std::source_location& loc = std::source_location::current());
LLmapError ConfigError(ErrorCode code, std::string_view key,
                       const std::source_location& loc = std::source_location::current());
LLmapError ValidationError(ErrorCode code, std::string_view what,
                           const std::source_location& loc = std::source_location::current());

// Result<T, E> type
template <typename T, typename E = LLmapError>
class Result {
public:
    static_assert(!std::is_same_v<T, E>, "T and E must be different types");

    Result() : storage_(T{}) {}

    Result(const T& value) : storage_(value) {}
    Result(T&& value) : storage_(std::move(value)) {}

    Result(const E& error) : storage_(error) {}
    Result(E&& error) : storage_(std::move(error)) {}

    bool ok() const noexcept { return std::holds_alternative<T>(storage_); }
    bool is_err() const noexcept { return !ok(); }
    explicit operator bool() const noexcept { return ok(); }

    T& value() & { return std::get<T>(storage_); }
    const T& value() const& { return std::get<T>(storage_); }
    T&& value() && { return std::move(std::get<T>(storage_)); }
    const T&& value() const&& { return std::move(std::get<T>(storage_)); }

    E& error() & { return std::get<E>(storage_); }
    const E& error() const& { return std::get<E>(storage_); }
    E&& error() && { return std::move(std::get<E>(storage_)); }
    const E&& error() const&& { return std::move(std::get<E>(storage_)); }

    T* operator->() { return &value(); }
    const T* operator->() const { return &value(); }

    T& operator*() & { return value(); }
    const T& operator*() const& { return value(); }
    T&& operator*() && { return std::move(value()); }

    T value_or(T&& default_val) const& {
        return ok() ? value() : std::forward<T>(default_val);
    }

    T value_or(T&& default_val) && {
        return ok() ? std::move(value()) : std::forward<T>(default_val);
    }

    template <typename F>
    auto map(F&& f) -> Result<std::invoke_result_t<F, T>, E> {
        using U = std::invoke_result_t<F, T>;
        if (ok()) {
            return Result<U, E>(std::forward<F>(f)(value()));
        }
        return Result<U, E>(error());
    }

    template <typename F>
    auto map(F&& f) const -> Result<std::invoke_result_t<F, const T&>, E> {
        using U = std::invoke_result_t<F, const T&>;
        if (ok()) {
            return Result<U, E>(std::forward<F>(f)(value()));
        }
        return Result<U, E>(error());
    }

    template <typename F>
    auto and_then(F&& f) -> std::invoke_result_t<F, T> {
        if (ok()) {
            return std::forward<F>(f)(value());
        }
        using RetType = std::invoke_result_t<F, T>;
        return RetType(error());
    }

    template <typename F>
    auto or_else(F&& f) -> Result<T, std::invoke_result_t<F, E>> {
        if (ok()) {
            return *this;
        }
        return std::forward<F>(f)(error());
    }

    template <typename F>
    Result& inspect(F&& f) {
        if (ok()) {
            std::forward<F>(f)(value());
        }
        return *this;
    }

    template <typename F>
    const Result& inspect(F&& f) const {
        if (ok()) {
            std::forward<F>(f)(value());
        }
        return *this;
    }

    template <typename F>
    Result& inspect_err(F&& f) {
        if (is_err()) {
            std::forward<F>(f)(error());
        }
        return *this;
    }

private:
    std::variant<T, E> storage_;
};

// Specialization for void success type
template <typename E>
class Result<void, E> {
public:
    Result() : error_(std::nullopt) {}

    Result(const E& error) : error_(error) {}
    Result(E&& error) : error_(std::move(error)) {}

    bool ok() const noexcept { return !error_.has_value(); }
    bool is_err() const noexcept { return error_.has_value(); }
    explicit operator bool() const noexcept { return ok(); }

    E& error() { return error_.value(); }
    const E& error() const { return error_.value(); }

    template <typename F>
    auto and_then(F&& f) -> std::invoke_result_t<F> {
        if (ok()) {
            return std::forward<F>(f)();
        }
        using RetType = std::invoke_result_t<F>;
        return RetType(error_.value());
    }

private:
    std::optional<E> error_;
};

// Helper for returning success
struct Ok {};

// Helper factory functions
template <typename T>
Result<std::decay_t<T>, LLmapError> MakeOk(T&& value) {
    return Result<std::decay_t<T>, LLmapError>(std::forward<T>(value));
}

inline Result<void, LLmapError> MakeOk() {
    return Result<void, LLmapError>();
}

template <typename T = void>
Result<T, LLmapError> MakeErr(LLmapError error) {
    return Result<T, LLmapError>(std::move(error));
}

template <typename T = void>
Result<T, LLmapError> MakeErr(ErrorCode code, std::string message = "",
                              const std::source_location& loc = std::source_location::current()) {
    return Result<T, LLmapError>(LLmapError(code, std::move(message), loc));
}

// Macro for early return on error
#define LLMAP_TRY(expr)                          \
    do {                                         \
        auto&& _result = (expr);                 \
        if (!_result.ok()) {                     \
            return std::move(_result).error();   \
        }                                        \
    } while (false)

#define LLMAP_TRY_ASSIGN(var, expr)              \
    auto&& _tmp_##var = (expr);                  \
    if (!_tmp_##var.ok()) {                      \
        return std::move(_tmp_##var).error();    \
    }                                            \
    var = std::move(*_tmp_##var)

// Aggregate multiple errors
class ErrorList {
public:
    void Add(LLmapError error) {
        if (error) {
            errors_.push_back(std::move(error));
        }
    }

    bool empty() const noexcept { return errors_.empty(); }
    size_t size() const noexcept { return errors_.size(); }
    bool has_errors() const noexcept { return !empty(); }

    const std::vector<LLmapError>& errors() const noexcept { return errors_; }

    LLmapError first() const {
        if (empty()) {
            return LLmapError();
        }
        return errors_.front();
    }

    Result<void, LLmapError> ToResult() const {
        if (empty()) {
            return MakeOk();
        }
        return MakeErr<void>(first());
    }

    std::string ToString() const;

private:
    std::vector<LLmapError> errors_;
};

}  // namespace llmap
