// LLmap — Error handling framework implementation.

#include "error.h"

#include <cstring>
#include <sstream>

namespace llmap {

LLmapError& LLmapError::WithOsError(int errno_val) {
    os_errno_ = errno_val;
    return *this;
}

std::string LLmapError::ToString() const {
    std::ostringstream oss;
    oss << "[" << ErrorCodeName(code_) << "]";

    if (!message_.empty()) {
        oss << " " << message_;
    }

    if (!context_.empty()) {
        oss << " (context: " << context_ << ")";
    }

    if (os_errno_ != 0) {
        oss << " (errno=" << os_errno_ << ": " << std::strerror(os_errno_) << ")";
    }

    if (location_.line > 0) {
        oss << " at " << location_.file << ":" << location_.line;
    }

    return oss.str();
}

std::string LLmapError::ToJson() const {
    std::ostringstream oss;
    oss << "{";
    oss << "\"code\":" << static_cast<uint16_t>(code_);
    oss << ",\"code_name\":\"" << ErrorCodeName(code_) << "\"";

    if (!message_.empty()) {
        oss << ",\"message\":\"";
        for (char c : message_) {
            if (c == '"') oss << "\\\"";
            else if (c == '\\') oss << "\\\\";
            else if (c == '\n') oss << "\\n";
            else if (c == '\r') oss << "\\r";
            else if (c == '\t') oss << "\\t";
            else oss << c;
        }
        oss << "\"";
    }

    if (!context_.empty()) {
        oss << ",\"context\":\"";
        for (char c : context_) {
            if (c == '"') oss << "\\\"";
            else if (c == '\\') oss << "\\\\";
            else if (c == '\n') oss << "\\n";
            else if (c == '\r') oss << "\\r";
            else if (c == '\t') oss << "\\t";
            else oss << c;
        }
        oss << "\"";
    }

    if (os_errno_ != 0) {
        oss << ",\"errno\":" << os_errno_;
        oss << ",\"errno_str\":\"" << std::strerror(os_errno_) << "\"";
    }

    if (location_.line > 0) {
        oss << ",\"file\":\"" << location_.file << "\"";
        oss << ",\"line\":" << location_.line;
    }

    oss << "}";
    return oss.str();
}

LLmapError IoError(ErrorCode code, std::string_view path,
                   const std::source_location& loc) {
    std::string message;
    switch (code) {
        case ErrorCode::kIoFileNotFound:
            message = "File not found: " + std::string(path);
            break;
        case ErrorCode::kIoPermissionDenied:
            message = "Permission denied: " + std::string(path);
            break;
        case ErrorCode::kIoReadError:
            message = "Read error: " + std::string(path);
            break;
        case ErrorCode::kIoWriteError:
            message = "Write error: " + std::string(path);
            break;
        case ErrorCode::kIoInvalidPath:
            message = "Invalid path: " + std::string(path);
            break;
        default:
            message = "I/O error: " + std::string(path);
            break;
    }
    return LLmapError(code, std::move(message), loc);
}

LLmapError ParseError(ErrorCode code, std::string_view what, size_t line,
                      const std::source_location& loc) {
    std::ostringstream oss;
    switch (code) {
        case ErrorCode::kParseInvalidFormat:
            oss << "Invalid format";
            break;
        case ErrorCode::kParseUnexpectedToken:
            oss << "Unexpected token";
            break;
        case ErrorCode::kParseMissingField:
            oss << "Missing field";
            break;
        case ErrorCode::kParseInvalidValue:
            oss << "Invalid value";
            break;
        case ErrorCode::kParseCorruptedData:
            oss << "Corrupted data";
            break;
        case ErrorCode::kParseVersionMismatch:
            oss << "Version mismatch";
            break;
        default:
            oss << "Parse error";
            break;
    }

    if (!what.empty()) {
        oss << ": " << what;
    }

    if (line > 0) {
        oss << " at line " << line;
    }

    return LLmapError(code, oss.str(), loc);
}

LLmapError ConfigError(ErrorCode code, std::string_view key,
                       const std::source_location& loc) {
    std::string message;
    switch (code) {
        case ErrorCode::kConfigMissingRequired:
            message = "Missing required config: " + std::string(key);
            break;
        case ErrorCode::kConfigInvalidValue:
            message = "Invalid config value: " + std::string(key);
            break;
        case ErrorCode::kConfigConflict:
            message = "Conflicting config: " + std::string(key);
            break;
        case ErrorCode::kConfigUnsupported:
            message = "Unsupported config: " + std::string(key);
            break;
        default:
            message = "Config error: " + std::string(key);
            break;
    }
    return LLmapError(code, std::move(message), loc);
}

LLmapError ValidationError(ErrorCode code, std::string_view what,
                           const std::source_location& loc) {
    std::string message;
    switch (code) {
        case ErrorCode::kValidateEmpty:
            message = "Empty " + std::string(what);
            break;
        case ErrorCode::kValidateOutOfRange:
            message = "Out of range: " + std::string(what);
            break;
        case ErrorCode::kValidateInvariantViolation:
            message = "Invariant violation: " + std::string(what);
            break;
        case ErrorCode::kValidateChecksum:
            message = "Checksum mismatch: " + std::string(what);
            break;
        case ErrorCode::kValidateDuplicate:
            message = "Duplicate: " + std::string(what);
            break;
        default:
            message = "Validation error: " + std::string(what);
            break;
    }
    return LLmapError(code, std::move(message), loc);
}

std::string ErrorList::ToString() const {
    if (empty()) {
        return "No errors";
    }

    std::ostringstream oss;
    oss << errors_.size() << " error(s):\n";
    for (size_t i = 0; i < errors_.size(); ++i) {
        oss << "  " << (i + 1) << ". " << errors_[i].ToString() << "\n";
    }
    return oss.str();
}

}  // namespace llmap
