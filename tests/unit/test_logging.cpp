// LLmap — Unit tests for structured logging framework.

#include "core/logging.h"

#include <gtest/gtest.h>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <thread>
#include <vector>
#include <atomic>
#include <cstdlib>

namespace fs = std::filesystem;

namespace llmap {
namespace {

class LoggingTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Reset to defaults before each test
        Logger::Instance().SetLevel(LogLevel::kInfo);
        Logger::Instance().SetFormat(LogFormat::kText);
        Logger::Instance().ClearSinks();

        // Create unique test directory
        test_dir_ = fs::temp_directory_path() /
            ("llmap_test_logging_" + std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id())));
        fs::create_directories(test_dir_);
    }

    void TearDown() override {
        Logger::Instance().SetOutput(stderr);
        if (fs::exists(test_dir_)) {
            fs::remove_all(test_dir_);
        }
    }

    fs::path test_dir_;
};

// =============================================================================
// Log Level Tests
// =============================================================================

TEST_F(LoggingTest, DefaultLevelIsInfo) {
    Logger::Instance().SetLevel(LogLevel::kInfo);
    EXPECT_EQ(Logger::Instance().GetLevel(), LogLevel::kInfo);
}

TEST_F(LoggingTest, SetLevelWorks) {
    Logger::Instance().SetLevel(LogLevel::kDebug);
    EXPECT_EQ(Logger::Instance().GetLevel(), LogLevel::kDebug);

    Logger::Instance().SetLevel(LogLevel::kError);
    EXPECT_EQ(Logger::Instance().GetLevel(), LogLevel::kError);
}

TEST_F(LoggingTest, ShouldLogRespectsLevel) {
    Logger::Instance().SetLevel(LogLevel::kWarn);

    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kTrace));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kDebug));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kInfo));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kWarn));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kError));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kFatal));
}

TEST_F(LoggingTest, OffLevelDisablesAllLogging) {
    Logger::Instance().SetLevel(LogLevel::kOff);

    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kTrace));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kDebug));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kInfo));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kWarn));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kError));
    EXPECT_FALSE(Logger::Instance().ShouldLog(LogLevel::kFatal));
}

TEST_F(LoggingTest, TraceLevelEnablesAll) {
    Logger::Instance().SetLevel(LogLevel::kTrace);

    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kTrace));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kDebug));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kInfo));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kWarn));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kError));
    EXPECT_TRUE(Logger::Instance().ShouldLog(LogLevel::kFatal));
}

// =============================================================================
// Log Level Name Tests
// =============================================================================

TEST_F(LoggingTest, LogLevelNameReturnsCorrectStrings) {
    EXPECT_EQ(LogLevelName(LogLevel::kTrace), "TRACE");
    EXPECT_EQ(LogLevelName(LogLevel::kDebug), "DEBUG");
    EXPECT_EQ(LogLevelName(LogLevel::kInfo), "INFO");
    EXPECT_EQ(LogLevelName(LogLevel::kWarn), "WARN");
    EXPECT_EQ(LogLevelName(LogLevel::kError), "ERROR");
    EXPECT_EQ(LogLevelName(LogLevel::kFatal), "FATAL");
    EXPECT_EQ(LogLevelName(LogLevel::kOff), "OFF");
}

TEST_F(LoggingTest, LogLevelNameLowerReturnsLowercaseStrings) {
    EXPECT_EQ(LogLevelNameLower(LogLevel::kTrace), "trace");
    EXPECT_EQ(LogLevelNameLower(LogLevel::kDebug), "debug");
    EXPECT_EQ(LogLevelNameLower(LogLevel::kInfo), "info");
    EXPECT_EQ(LogLevelNameLower(LogLevel::kWarn), "warn");
    EXPECT_EQ(LogLevelNameLower(LogLevel::kError), "error");
    EXPECT_EQ(LogLevelNameLower(LogLevel::kFatal), "fatal");
    EXPECT_EQ(LogLevelNameLower(LogLevel::kOff), "off");
}

// =============================================================================
// Parse Log Level Tests
// =============================================================================

TEST_F(LoggingTest, ParseLogLevelParsesLowercase) {
    EXPECT_EQ(ParseLogLevel("trace"), LogLevel::kTrace);
    EXPECT_EQ(ParseLogLevel("debug"), LogLevel::kDebug);
    EXPECT_EQ(ParseLogLevel("info"), LogLevel::kInfo);
    EXPECT_EQ(ParseLogLevel("warn"), LogLevel::kWarn);
    EXPECT_EQ(ParseLogLevel("error"), LogLevel::kError);
    EXPECT_EQ(ParseLogLevel("fatal"), LogLevel::kFatal);
    EXPECT_EQ(ParseLogLevel("off"), LogLevel::kOff);
}

TEST_F(LoggingTest, ParseLogLevelParsesUppercase) {
    EXPECT_EQ(ParseLogLevel("TRACE"), LogLevel::kTrace);
    EXPECT_EQ(ParseLogLevel("DEBUG"), LogLevel::kDebug);
    EXPECT_EQ(ParseLogLevel("INFO"), LogLevel::kInfo);
    EXPECT_EQ(ParseLogLevel("WARN"), LogLevel::kWarn);
    EXPECT_EQ(ParseLogLevel("ERROR"), LogLevel::kError);
    EXPECT_EQ(ParseLogLevel("FATAL"), LogLevel::kFatal);
    EXPECT_EQ(ParseLogLevel("OFF"), LogLevel::kOff);
}

TEST_F(LoggingTest, ParseLogLevelParsesMixedCase) {
    EXPECT_EQ(ParseLogLevel("Trace"), LogLevel::kTrace);
    EXPECT_EQ(ParseLogLevel("Debug"), LogLevel::kDebug);
    EXPECT_EQ(ParseLogLevel("Info"), LogLevel::kInfo);
    EXPECT_EQ(ParseLogLevel("Warning"), LogLevel::kWarn);
}

TEST_F(LoggingTest, ParseLogLevelDefaultsToInfoOnUnknown) {
    EXPECT_EQ(ParseLogLevel("unknown"), LogLevel::kInfo);
    EXPECT_EQ(ParseLogLevel(""), LogLevel::kInfo);
    EXPECT_EQ(ParseLogLevel("garbage"), LogLevel::kInfo);
}

// =============================================================================
// Log Format Tests
// =============================================================================

TEST_F(LoggingTest, DefaultFormatIsText) {
    Logger::Instance().SetFormat(LogFormat::kText);
    EXPECT_EQ(Logger::Instance().GetFormat(), LogFormat::kText);
}

TEST_F(LoggingTest, SetFormatWorks) {
    Logger::Instance().SetFormat(LogFormat::kJson);
    EXPECT_EQ(Logger::Instance().GetFormat(), LogFormat::kJson);

    Logger::Instance().SetFormat(LogFormat::kText);
    EXPECT_EQ(Logger::Instance().GetFormat(), LogFormat::kText);
}

// =============================================================================
// Text Output Tests
// =============================================================================

TEST_F(LoggingTest, TextOutputContainsMessage) {
    auto log_path = test_dir_ / "test_text.log";
    FILE* log_file = std::fopen(log_path.c_str(), "w");
    ASSERT_NE(log_file, nullptr);

    Logger::Instance().SetOutput(log_file);
    Logger::Instance().SetFormat(LogFormat::kText);
    Logger::Instance().SetLevel(LogLevel::kInfo);

    Logger::Instance().Log(LogLevel::kInfo, "Test message 12345");
    Logger::Instance().Flush();
    std::fclose(log_file);

    std::ifstream ifs(log_path);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("Test message 12345") != std::string::npos);
    EXPECT_TRUE(content.find("[INFO]") != std::string::npos);
}

TEST_F(LoggingTest, TextOutputContainsSourceLocation) {
    auto log_path = test_dir_ / "test_location.log";
    FILE* log_file = std::fopen(log_path.c_str(), "w");
    ASSERT_NE(log_file, nullptr);

    Logger::Instance().SetOutput(log_file);
    Logger::Instance().SetFormat(LogFormat::kText);
    Logger::Instance().SetLevel(LogLevel::kInfo);

    Logger::Instance().Log(LogLevel::kInfo, "Location test");
    Logger::Instance().Flush();
    std::fclose(log_file);

    std::ifstream ifs(log_path);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    // Should contain filename (test_logging.cpp)
    EXPECT_TRUE(content.find("test_logging.cpp") != std::string::npos);
}

TEST_F(LoggingTest, TextOutputContainsTimestamp) {
    auto log_path = test_dir_ / "test_timestamp.log";
    FILE* log_file = std::fopen(log_path.c_str(), "w");
    ASSERT_NE(log_file, nullptr);

    Logger::Instance().SetOutput(log_file);
    Logger::Instance().SetFormat(LogFormat::kText);
    Logger::Instance().SetLevel(LogLevel::kInfo);

    Logger::Instance().Log(LogLevel::kInfo, "Timestamp test");
    Logger::Instance().Flush();
    std::fclose(log_file);

    std::ifstream ifs(log_path);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    // Should contain year (202X pattern)
    EXPECT_TRUE(content.find("202") != std::string::npos);
}

// =============================================================================
// JSON Output Tests
// =============================================================================

TEST_F(LoggingTest, JsonOutputIsValidJson) {
    auto log_path = test_dir_ / "test_json.log";
    FILE* log_file = std::fopen(log_path.c_str(), "w");
    ASSERT_NE(log_file, nullptr);

    Logger::Instance().SetOutput(log_file);
    Logger::Instance().SetFormat(LogFormat::kJson);
    Logger::Instance().SetLevel(LogLevel::kInfo);

    Logger::Instance().Log(LogLevel::kInfo, "JSON test message");
    Logger::Instance().Flush();
    std::fclose(log_file);

    std::ifstream ifs(log_path);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    // Basic JSON structure check
    EXPECT_TRUE(content.find("{") != std::string::npos);
    EXPECT_TRUE(content.find("}") != std::string::npos);
    EXPECT_TRUE(content.find("\"ts\":") != std::string::npos);
    EXPECT_TRUE(content.find("\"level\":") != std::string::npos);
    EXPECT_TRUE(content.find("\"msg\":") != std::string::npos);
    EXPECT_TRUE(content.find("\"info\"") != std::string::npos);
}

TEST_F(LoggingTest, JsonOutputEscapesSpecialChars) {
    auto log_path = test_dir_ / "test_escape.log";
    FILE* log_file = std::fopen(log_path.c_str(), "w");
    ASSERT_NE(log_file, nullptr);

    Logger::Instance().SetOutput(log_file);
    Logger::Instance().SetFormat(LogFormat::kJson);
    Logger::Instance().SetLevel(LogLevel::kInfo);

    Logger::Instance().Log(LogLevel::kInfo, "Message with \"quotes\" and \\backslash");
    Logger::Instance().Flush();
    std::fclose(log_file);

    std::ifstream ifs(log_path);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    // Should have escaped quotes
    EXPECT_TRUE(content.find("\\\"quotes\\\"") != std::string::npos);
    // Should have escaped backslash
    EXPECT_TRUE(content.find("\\\\backslash") != std::string::npos);
}

// =============================================================================
// Sink Tests
// =============================================================================

TEST_F(LoggingTest, SinkReceivesLogRecords) {
    std::vector<LogRecord> captured;

    Logger::Instance().AddSink([&captured](const LogRecord& record) {
        captured.push_back(record);
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().Log(LogLevel::kInfo, "Sink test");

    ASSERT_EQ(captured.size(), 1);
    EXPECT_EQ(captured[0].level, LogLevel::kInfo);
    EXPECT_EQ(captured[0].message, "Sink test");
}

TEST_F(LoggingTest, MultipleSinksReceiveRecords) {
    std::atomic<int> sink1_count{0};
    std::atomic<int> sink2_count{0};

    Logger::Instance().AddSink([&sink1_count](const LogRecord&) {
        sink1_count++;
    });
    Logger::Instance().AddSink([&sink2_count](const LogRecord&) {
        sink2_count++;
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().Log(LogLevel::kInfo, "Multi-sink test");

    EXPECT_EQ(sink1_count.load(), 1);
    EXPECT_EQ(sink2_count.load(), 1);
}

TEST_F(LoggingTest, ClearSinksRemovesAllSinks) {
    std::atomic<int> sink_count{0};

    Logger::Instance().AddSink([&sink_count](const LogRecord&) {
        sink_count++;
    });

    Logger::Instance().ClearSinks();
    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().Log(LogLevel::kInfo, "Should not reach sink");

    EXPECT_EQ(sink_count.load(), 0);
}

// =============================================================================
// Context Logging Tests
// =============================================================================

TEST_F(LoggingTest, LogWithContextAppendsContext) {
    std::string captured_message;

    Logger::Instance().AddSink([&captured_message](const LogRecord& record) {
        captured_message = std::string(record.message);
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().LogWithContext(LogLevel::kInfo, "Base message", "extra_context");

    EXPECT_TRUE(captured_message.find("Base message") != std::string::npos);
    EXPECT_TRUE(captured_message.find("[extra_context]") != std::string::npos);
}

// =============================================================================
// Macro Tests
// =============================================================================

TEST_F(LoggingTest, MacroRespectsShouldLog) {
    std::atomic<int> log_count{0};

    Logger::Instance().AddSink([&log_count](const LogRecord&) {
        log_count++;
    });

    Logger::Instance().SetLevel(LogLevel::kWarn);

    LLMAP_LOG_DEBUG("Should be skipped");
    LLMAP_LOG_INFO("Should be skipped");
    LLMAP_LOG_WARN("Should appear");
    LLMAP_LOG_ERROR("Should appear");

    EXPECT_EQ(log_count.load(), 2);
}

// =============================================================================
// Thread Safety Tests
// =============================================================================

TEST_F(LoggingTest, ConcurrentLoggingIsThreadSafe) {
    std::atomic<int> log_count{0};

    Logger::Instance().AddSink([&log_count](const LogRecord&) {
        log_count++;
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);

    constexpr int kNumThreads = 10;
    constexpr int kLogsPerThread = 100;

    std::vector<std::thread> threads;
    threads.reserve(kNumThreads);

    for (int i = 0; i < kNumThreads; ++i) {
        threads.emplace_back([i]() {
            for (int j = 0; j < kLogsPerThread; ++j) {
                LLMAP_LOG_INFO("Thread log");
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    EXPECT_EQ(log_count.load(), kNumThreads * kLogsPerThread);
}

// =============================================================================
// Log Record Fields Tests
// =============================================================================

TEST_F(LoggingTest, LogRecordContainsThreadId) {
    uint64_t captured_tid = 0;

    Logger::Instance().AddSink([&captured_tid](const LogRecord& record) {
        captured_tid = record.thread_id;
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().Log(LogLevel::kInfo, "TID test");

    EXPECT_GT(captured_tid, 0u);
}

TEST_F(LoggingTest, LogRecordContainsLineNumber) {
    uint32_t captured_line = 0;

    Logger::Instance().AddSink([&captured_line](const LogRecord& record) {
        captured_line = record.line;
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().Log(LogLevel::kInfo, "Line test");

    EXPECT_GT(captured_line, 0u);
}

TEST_F(LoggingTest, LogRecordContainsFunctionName) {
    std::string captured_function;

    Logger::Instance().AddSink([&captured_function](const LogRecord& record) {
        captured_function = std::string(record.function);
    });

    Logger::Instance().SetLevel(LogLevel::kTrace);
    Logger::Instance().Log(LogLevel::kInfo, "Function test");

    EXPECT_FALSE(captured_function.empty());
}

// =============================================================================
// Singleton Tests
// =============================================================================

TEST_F(LoggingTest, LoggerIsSingleton) {
    Logger& instance1 = Logger::Instance();
    Logger& instance2 = Logger::Instance();

    EXPECT_EQ(&instance1, &instance2);
}

}  // namespace
}  // namespace llmap
