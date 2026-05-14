// Tests for LLmap error handling framework.

#include <gtest/gtest.h>

#include "core/error.h"

#include <string>
#include <vector>

namespace llmap {
namespace {

// ============================================================================
// ErrorCode tests
// ============================================================================

TEST(ErrorCodeTest, ErrorCodeNameReturnsCorrectStrings) {
    EXPECT_EQ(ErrorCodeName(ErrorCode::kOk), "OK");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kIoFileNotFound), "IO_FILE_NOT_FOUND");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kParseInvalidFormat), "PARSE_INVALID_FORMAT");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kConfigMissingRequired), "CONFIG_MISSING_REQUIRED");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kValidateOutOfRange), "VALIDATE_OUT_OF_RANGE");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kResourceOutOfMemory), "RESOURCE_OUT_OF_MEMORY");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kSystemCudaError), "SYSTEM_CUDA_ERROR");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kAlgoConvergenceFailed), "ALGO_CONVERGENCE_FAILED");
    EXPECT_EQ(ErrorCodeName(ErrorCode::kExternalApiError), "EXTERNAL_API_ERROR");
}

TEST(ErrorCodeTest, IsIoErrorIdentifiesIoErrors) {
    EXPECT_TRUE(IsIoError(ErrorCode::kIoFileNotFound));
    EXPECT_TRUE(IsIoError(ErrorCode::kIoPermissionDenied));
    EXPECT_TRUE(IsIoError(ErrorCode::kIoReadError));
    EXPECT_FALSE(IsIoError(ErrorCode::kOk));
    EXPECT_FALSE(IsIoError(ErrorCode::kParseInvalidFormat));
}

TEST(ErrorCodeTest, IsParseErrorIdentifiesParseErrors) {
    EXPECT_TRUE(IsParseError(ErrorCode::kParseInvalidFormat));
    EXPECT_TRUE(IsParseError(ErrorCode::kParseCorruptedData));
    EXPECT_FALSE(IsParseError(ErrorCode::kIoFileNotFound));
}

TEST(ErrorCodeTest, IsConfigErrorIdentifiesConfigErrors) {
    EXPECT_TRUE(IsConfigError(ErrorCode::kConfigMissingRequired));
    EXPECT_TRUE(IsConfigError(ErrorCode::kConfigInvalidValue));
    EXPECT_FALSE(IsConfigError(ErrorCode::kParseInvalidFormat));
}

TEST(ErrorCodeTest, IsValidateErrorIdentifiesValidateErrors) {
    EXPECT_TRUE(IsValidateError(ErrorCode::kValidateEmpty));
    EXPECT_TRUE(IsValidateError(ErrorCode::kValidateOutOfRange));
    EXPECT_FALSE(IsValidateError(ErrorCode::kConfigMissingRequired));
}

TEST(ErrorCodeTest, IsResourceErrorIdentifiesResourceErrors) {
    EXPECT_TRUE(IsResourceError(ErrorCode::kResourceOutOfMemory));
    EXPECT_TRUE(IsResourceError(ErrorCode::kResourceTimeout));
    EXPECT_FALSE(IsResourceError(ErrorCode::kValidateEmpty));
}

TEST(ErrorCodeTest, IsSystemErrorIdentifiesSystemErrors) {
    EXPECT_TRUE(IsSystemError(ErrorCode::kSystemOsError));
    EXPECT_TRUE(IsSystemError(ErrorCode::kSystemCudaError));
    EXPECT_FALSE(IsSystemError(ErrorCode::kResourceOutOfMemory));
}

TEST(ErrorCodeTest, IsAlgoErrorIdentifiesAlgoErrors) {
    EXPECT_TRUE(IsAlgoError(ErrorCode::kAlgoConvergenceFailed));
    EXPECT_TRUE(IsAlgoError(ErrorCode::kAlgoNotImplemented));
    EXPECT_FALSE(IsAlgoError(ErrorCode::kSystemOsError));
}

TEST(ErrorCodeTest, IsExternalErrorIdentifiesExternalErrors) {
    EXPECT_TRUE(IsExternalError(ErrorCode::kExternalApiError));
    EXPECT_TRUE(IsExternalError(ErrorCode::kExternalNetworkError));
    EXPECT_FALSE(IsExternalError(ErrorCode::kAlgoConvergenceFailed));
}

// ============================================================================
// LLmapError tests
// ============================================================================

TEST(LLmapErrorTest, DefaultConstructorCreatesOkError) {
    LLmapError err;
    EXPECT_TRUE(err.ok());
    EXPECT_EQ(err.code(), ErrorCode::kOk);
}

TEST(LLmapErrorTest, ConstructWithCodeOnly) {
    LLmapError err(ErrorCode::kIoFileNotFound);
    EXPECT_FALSE(err.ok());
    EXPECT_TRUE(static_cast<bool>(err));
    EXPECT_EQ(err.code(), ErrorCode::kIoFileNotFound);
}

TEST(LLmapErrorTest, ConstructWithCodeAndMessage) {
    LLmapError err(ErrorCode::kParseInvalidFormat, "Bad JSON");
    EXPECT_FALSE(err.ok());
    EXPECT_EQ(err.code(), ErrorCode::kParseInvalidFormat);
    EXPECT_EQ(err.message(), "Bad JSON");
}

TEST(LLmapErrorTest, ConstructWithCodeMessageAndContext) {
    LLmapError err(ErrorCode::kConfigMissingRequired, "Missing key", "config.toml");
    EXPECT_EQ(err.code(), ErrorCode::kConfigMissingRequired);
    EXPECT_EQ(err.message(), "Missing key");
    EXPECT_EQ(err.context(), "config.toml");
}

TEST(LLmapErrorTest, WithContextAddsContext) {
    LLmapError err(ErrorCode::kValidateOutOfRange, "Value too large");
    err.WithContext("field: age");
    EXPECT_EQ(err.context(), "field: age");
}

TEST(LLmapErrorTest, ToStringFormatsError) {
    LLmapError err(ErrorCode::kIoReadError, "Failed to read");
    std::string s = err.ToString();
    EXPECT_TRUE(s.find("IO_READ_ERROR") != std::string::npos);
    EXPECT_TRUE(s.find("Failed to read") != std::string::npos);
}

TEST(LLmapErrorTest, ToStringIncludesContext) {
    LLmapError err(ErrorCode::kIoWriteError, "Write failed", "output.bam");
    std::string s = err.ToString();
    EXPECT_TRUE(s.find("output.bam") != std::string::npos);
}

TEST(LLmapErrorTest, ToJsonFormatsAsJson) {
    LLmapError err(ErrorCode::kParseInvalidValue, "Invalid number");
    std::string json = err.ToJson();
    EXPECT_TRUE(json.find("\"code\":") != std::string::npos);
    EXPECT_TRUE(json.find("PARSE_INVALID_VALUE") != std::string::npos);
    EXPECT_TRUE(json.find("\"message\":\"Invalid number\"") != std::string::npos);
}

TEST(LLmapErrorTest, ToJsonEscapesSpecialChars) {
    LLmapError err(ErrorCode::kParseInvalidFormat, "Bad \"quote\" and \\ backslash");
    std::string json = err.ToJson();
    EXPECT_TRUE(json.find("\\\"quote\\\"") != std::string::npos);
    EXPECT_TRUE(json.find("\\\\") != std::string::npos);
}

TEST(LLmapErrorTest, EqualityComparesCode) {
    LLmapError err1(ErrorCode::kIoFileNotFound);
    LLmapError err2(ErrorCode::kIoFileNotFound, "different message");
    LLmapError err3(ErrorCode::kIoReadError);

    EXPECT_EQ(err1, err2);
    EXPECT_NE(err1, err3);
}

TEST(LLmapErrorTest, CapturesSourceLocation) {
    LLmapError err(ErrorCode::kAlgoConvergenceFailed, "EM did not converge");
    const auto& loc = err.location();
    EXPECT_TRUE(loc.line > 0);
    EXPECT_TRUE(std::string(loc.file).find("test_error.cpp") != std::string::npos);
}

// ============================================================================
// Error factory function tests
// ============================================================================

TEST(ErrorFactoryTest, IoErrorCreatesIoError) {
    auto err = IoError(ErrorCode::kIoFileNotFound, "/path/to/file.fa");
    EXPECT_EQ(err.code(), ErrorCode::kIoFileNotFound);
    EXPECT_TRUE(std::string(err.message()).find("/path/to/file.fa") != std::string::npos);
}

TEST(ErrorFactoryTest, ParseErrorCreatesParseError) {
    auto err = ParseError(ErrorCode::kParseInvalidFormat, "expected number", 42);
    EXPECT_EQ(err.code(), ErrorCode::kParseInvalidFormat);
    EXPECT_TRUE(std::string(err.message()).find("line 42") != std::string::npos);
}

TEST(ErrorFactoryTest, ConfigErrorCreatesConfigError) {
    auto err = ConfigError(ErrorCode::kConfigMissingRequired, "api_key");
    EXPECT_EQ(err.code(), ErrorCode::kConfigMissingRequired);
    EXPECT_TRUE(std::string(err.message()).find("api_key") != std::string::npos);
}

TEST(ErrorFactoryTest, ValidationErrorCreatesValidationError) {
    auto err = ValidationError(ErrorCode::kValidateOutOfRange, "quality score");
    EXPECT_EQ(err.code(), ErrorCode::kValidateOutOfRange);
    EXPECT_TRUE(std::string(err.message()).find("quality score") != std::string::npos);
}

// ============================================================================
// Result<T, E> tests
// ============================================================================

TEST(ResultTest, DefaultConstructorCreatesOkResult) {
    Result<int> r;
    EXPECT_TRUE(r.ok());
    EXPECT_EQ(r.value(), 0);
}

TEST(ResultTest, ConstructFromValueCreatesOkResult) {
    Result<int> r(42);
    EXPECT_TRUE(r.ok());
    EXPECT_FALSE(r.is_err());
    EXPECT_EQ(r.value(), 42);
}

TEST(ResultTest, ConstructFromErrorCreatesErrResult) {
    Result<int> r(LLmapError(ErrorCode::kIoReadError));
    EXPECT_FALSE(r.ok());
    EXPECT_TRUE(r.is_err());
    EXPECT_EQ(r.error().code(), ErrorCode::kIoReadError);
}

TEST(ResultTest, BoolConversionWorksCorrectly) {
    Result<int> ok_result(10);
    Result<int> err_result(LLmapError(ErrorCode::kAlgoInvalidState));

    EXPECT_TRUE(static_cast<bool>(ok_result));
    EXPECT_FALSE(static_cast<bool>(err_result));
}

TEST(ResultTest, OperatorArrowAccessesValue) {
    struct Point { int x, y; };
    Result<Point> r(Point{3, 4});
    EXPECT_EQ(r->x, 3);
    EXPECT_EQ(r->y, 4);
}

TEST(ResultTest, OperatorStarAccessesValue) {
    Result<std::string> r("hello");
    EXPECT_EQ(*r, "hello");
}

TEST(ResultTest, ValueOrReturnsValueWhenOk) {
    Result<int> r(42);
    EXPECT_EQ(r.value_or(99), 42);
}

TEST(ResultTest, ValueOrReturnsDefaultWhenErr) {
    Result<int> r(LLmapError(ErrorCode::kResourceTimeout));
    EXPECT_EQ(r.value_or(99), 99);
}

TEST(ResultTest, MapTransformsValue) {
    Result<int> r(10);
    auto r2 = r.map([](int x) { return x * 2; });
    EXPECT_TRUE(r2.ok());
    EXPECT_EQ(r2.value(), 20);
}

TEST(ResultTest, MapPreservesError) {
    Result<int> r(LLmapError(ErrorCode::kAlgoConvergenceFailed));
    auto r2 = r.map([](int x) { return x * 2; });
    EXPECT_FALSE(r2.ok());
    EXPECT_EQ(r2.error().code(), ErrorCode::kAlgoConvergenceFailed);
}

TEST(ResultTest, AndThenChainsOnSuccess) {
    auto parse_int = [](const std::string& s) -> Result<int> {
        try {
            return std::stoi(s);
        } catch (...) {
            return LLmapError(ErrorCode::kParseInvalidValue, "Not a number");
        }
    };

    Result<std::string> r("123");
    auto r2 = r.and_then([&](const std::string& s) { return parse_int(s); });
    EXPECT_TRUE(r2.ok());
    EXPECT_EQ(r2.value(), 123);
}

TEST(ResultTest, AndThenShortCircuitsOnError) {
    auto never_called = [](const std::string&) -> Result<int> {
        ADD_FAILURE() << "Should not be called";
        return 0;
    };

    Result<std::string> r(LLmapError(ErrorCode::kIoReadError));
    auto r2 = r.and_then(never_called);
    EXPECT_FALSE(r2.ok());
}

TEST(ResultTest, InspectCallsOnSuccess) {
    bool called = false;
    Result<int> r(42);
    r.inspect([&](int x) {
        called = true;
        EXPECT_EQ(x, 42);
    });
    EXPECT_TRUE(called);
}

TEST(ResultTest, InspectSkipsOnError) {
    bool called = false;
    Result<int> r(LLmapError(ErrorCode::kUnknown));
    r.inspect([&](int) { called = true; });
    EXPECT_FALSE(called);
}

TEST(ResultTest, InspectErrCallsOnError) {
    bool called = false;
    Result<int> r(LLmapError(ErrorCode::kResourceBusy));
    r.inspect_err([&](const LLmapError& e) {
        called = true;
        EXPECT_EQ(e.code(), ErrorCode::kResourceBusy);
    });
    EXPECT_TRUE(called);
}

// ============================================================================
// Result<void, E> specialization tests
// ============================================================================

TEST(ResultVoidTest, DefaultConstructorCreatesOkResult) {
    Result<void> r;
    EXPECT_TRUE(r.ok());
}

TEST(ResultVoidTest, ConstructFromErrorCreatesErrResult) {
    Result<void> r(LLmapError(ErrorCode::kIoWriteError));
    EXPECT_FALSE(r.ok());
    EXPECT_TRUE(r.is_err());
}

TEST(ResultVoidTest, AndThenChainsOnSuccess) {
    bool called = false;
    Result<void> r;
    auto r2 = r.and_then([&]() -> Result<int> {
        called = true;
        return 42;
    });
    EXPECT_TRUE(called);
    EXPECT_TRUE(r2.ok());
    EXPECT_EQ(r2.value(), 42);
}

// ============================================================================
// MakeOk / MakeErr tests
// ============================================================================

TEST(MakeResultTest, MakeOkCreatesOkResult) {
    auto r = MakeOk(42);
    EXPECT_TRUE(r.ok());
    EXPECT_EQ(r.value(), 42);
}

TEST(MakeResultTest, MakeOkVoidCreatesOkVoidResult) {
    auto r = MakeOk();
    EXPECT_TRUE(r.ok());
}

TEST(MakeResultTest, MakeErrCreatesErrResult) {
    auto r = MakeErr<int>(LLmapError(ErrorCode::kResourceExhausted));
    EXPECT_FALSE(r.ok());
    EXPECT_EQ(r.error().code(), ErrorCode::kResourceExhausted);
}

TEST(MakeResultTest, MakeErrWithCodeAndMessage) {
    auto r = MakeErr<std::string>(ErrorCode::kAlgoNotImplemented, "Feature X");
    EXPECT_FALSE(r.ok());
    EXPECT_EQ(r.error().code(), ErrorCode::kAlgoNotImplemented);
    EXPECT_EQ(r.error().message(), "Feature X");
}

// ============================================================================
// ErrorList tests
// ============================================================================

TEST(ErrorListTest, DefaultConstructorCreatesEmptyList) {
    ErrorList list;
    EXPECT_TRUE(list.empty());
    EXPECT_EQ(list.size(), 0u);
    EXPECT_FALSE(list.has_errors());
}

TEST(ErrorListTest, AddErrorIncreasesSize) {
    ErrorList list;
    list.Add(LLmapError(ErrorCode::kIoReadError));
    list.Add(LLmapError(ErrorCode::kIoWriteError));
    EXPECT_EQ(list.size(), 2u);
    EXPECT_TRUE(list.has_errors());
}

TEST(ErrorListTest, AddOkErrorDoesNotIncreaseSize) {
    ErrorList list;
    list.Add(LLmapError());
    EXPECT_TRUE(list.empty());
}

TEST(ErrorListTest, FirstReturnsFirstError) {
    ErrorList list;
    list.Add(LLmapError(ErrorCode::kIoReadError, "first"));
    list.Add(LLmapError(ErrorCode::kIoWriteError, "second"));
    EXPECT_EQ(list.first().code(), ErrorCode::kIoReadError);
}

TEST(ErrorListTest, FirstReturnsOkOnEmptyList) {
    ErrorList list;
    EXPECT_TRUE(list.first().ok());
}

TEST(ErrorListTest, ToResultSucceedsOnEmptyList) {
    ErrorList list;
    auto r = list.ToResult();
    EXPECT_TRUE(r.ok());
}

TEST(ErrorListTest, ToResultFailsOnNonEmptyList) {
    ErrorList list;
    list.Add(LLmapError(ErrorCode::kValidateChecksum));
    auto r = list.ToResult();
    EXPECT_FALSE(r.ok());
    EXPECT_EQ(r.error().code(), ErrorCode::kValidateChecksum);
}

TEST(ErrorListTest, ToStringFormatsErrors) {
    ErrorList list;
    list.Add(LLmapError(ErrorCode::kIoReadError, "read fail"));
    list.Add(LLmapError(ErrorCode::kIoWriteError, "write fail"));
    std::string s = list.ToString();
    EXPECT_TRUE(s.find("2 error(s)") != std::string::npos);
    EXPECT_TRUE(s.find("read fail") != std::string::npos);
    EXPECT_TRUE(s.find("write fail") != std::string::npos);
}

// ============================================================================
// LLMAP_TRY macro tests
// ============================================================================

namespace {

Result<int> divide(int a, int b) {
    if (b == 0) {
        return LLmapError(ErrorCode::kAlgoNumericalInstability, "Division by zero");
    }
    return a / b;
}

Result<int> compute_with_try(int x) {
    LLMAP_TRY(divide(x, 0));
    return 42;
}

Result<int> compute_success(int x) {
    auto result = divide(x, 2);
    LLMAP_TRY(result);
    return result.value() + 10;
}

}  // namespace

TEST(TryMacroTest, TryReturnsOnError) {
    auto r = compute_with_try(10);
    EXPECT_FALSE(r.ok());
    EXPECT_EQ(r.error().code(), ErrorCode::kAlgoNumericalInstability);
}

TEST(TryMacroTest, TryContinuesOnSuccess) {
    auto r = compute_success(20);
    EXPECT_TRUE(r.ok());
    EXPECT_EQ(r.value(), 20);
}

// ============================================================================
// Complex type tests
// ============================================================================

TEST(ResultComplexTest, WorksWithMoveOnlyTypes) {
    Result<std::unique_ptr<int>> r(std::make_unique<int>(42));
    EXPECT_TRUE(r.ok());
    auto ptr = std::move(r).value();
    EXPECT_EQ(*ptr, 42);
}

TEST(ResultComplexTest, WorksWithVectors) {
    Result<std::vector<int>> r(std::vector<int>{1, 2, 3});
    EXPECT_TRUE(r.ok());
    EXPECT_EQ(r->size(), 3u);
}

TEST(ResultComplexTest, MapChangesType) {
    Result<std::string> r("hello");
    auto r2 = r.map([](const std::string& s) { return s.size(); });
    EXPECT_TRUE(r2.ok());
    EXPECT_EQ(r2.value(), 5u);
}

}  // namespace
}  // namespace llmap
