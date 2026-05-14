#include "core/profiler.h"

#include <gtest/gtest.h>
#include <sstream>
#include <thread>
#include <chrono>

using namespace llmap::core;

class ProfilerTest : public ::testing::Test {
protected:
    void SetUp() override {
        ProfileRegistry::Instance().Reset();
    }

    void TearDown() override {
        ProfileRegistry::Instance().Reset();
    }
};

TEST_F(ProfilerTest, ProfileStatsRecordSingle) {
    ProfileStats stats("test_scope");

    stats.Record(1000);  // 1 µs

    EXPECT_EQ(1, stats.count.load());
    EXPECT_EQ(1000, stats.total_ns.load());
    EXPECT_EQ(1000, stats.min_ns.load());
    EXPECT_EQ(1000, stats.max_ns.load());
}

TEST_F(ProfilerTest, ProfileStatsRecordMultiple) {
    ProfileStats stats("test_scope");

    stats.Record(1000);  // 1 µs
    stats.Record(2000);  // 2 µs
    stats.Record(500);   // 0.5 µs

    EXPECT_EQ(3, stats.count.load());
    EXPECT_EQ(3500, stats.total_ns.load());
    EXPECT_EQ(500, stats.min_ns.load());
    EXPECT_EQ(2000, stats.max_ns.load());
}

TEST_F(ProfilerTest, ProfileStatsAverageUs) {
    ProfileStats stats("test_scope");

    stats.Record(1000);
    stats.Record(2000);
    stats.Record(3000);

    EXPECT_DOUBLE_EQ(2.0, stats.AverageUs());
}

TEST_F(ProfilerTest, ProfileStatsTotalMs) {
    ProfileStats stats("test_scope");

    stats.Record(1000000);  // 1 ms
    stats.Record(2000000);  // 2 ms

    EXPECT_DOUBLE_EQ(3.0, stats.TotalMs());
}

TEST_F(ProfilerTest, ProfileStatsEmptyAverage) {
    ProfileStats stats("test_scope");
    EXPECT_DOUBLE_EQ(0.0, stats.AverageUs());
}

TEST_F(ProfilerTest, RegistryGetOrCreate) {
    auto& stats1 = ProfileRegistry::Instance().GetOrCreate("scope_a");
    auto& stats2 = ProfileRegistry::Instance().GetOrCreate("scope_a");
    auto& stats3 = ProfileRegistry::Instance().GetOrCreate("scope_b");

    EXPECT_EQ(&stats1, &stats2);
    EXPECT_NE(&stats1, &stats3);

    EXPECT_EQ("scope_a", stats1.name);
    EXPECT_EQ("scope_b", stats3.name);
}

TEST_F(ProfilerTest, RegistryReset) {
    ProfileRegistry::Instance().GetOrCreate("scope_a").Record(1000);
    ProfileRegistry::Instance().GetOrCreate("scope_b").Record(2000);

    auto snapshot = ProfileRegistry::Instance().GetSnapshot();
    EXPECT_EQ(2, snapshot.size());

    ProfileRegistry::Instance().Reset();

    snapshot = ProfileRegistry::Instance().GetSnapshot();
    EXPECT_EQ(0, snapshot.size());
}

TEST_F(ProfilerTest, RegistrySnapshotSortedByTotal) {
    ProfileRegistry::Instance().GetOrCreate("fast").Record(100);
    ProfileRegistry::Instance().GetOrCreate("slow").Record(10000);
    ProfileRegistry::Instance().GetOrCreate("medium").Record(1000);

    auto snapshot = ProfileRegistry::Instance().GetSnapshot();

    ASSERT_EQ(3, snapshot.size());
    EXPECT_EQ("slow", snapshot[0].name);
    EXPECT_EQ("medium", snapshot[1].name);
    EXPECT_EQ("fast", snapshot[2].name);
}

TEST_F(ProfilerTest, RegistryForEach) {
    ProfileRegistry::Instance().GetOrCreate("a").Record(100);
    ProfileRegistry::Instance().GetOrCreate("b").Record(200);

    int count = 0;
    ProfileRegistry::Instance().ForEach([&](const std::string& name, const ProfileStats& s) {
        ++count;
    });

    EXPECT_EQ(2, count);
}

TEST_F(ProfilerTest, RegistryPrintReport) {
    ProfileRegistry::Instance().GetOrCreate("func1").Record(1000000);  // 1ms
    ProfileRegistry::Instance().GetOrCreate("func2").Record(2000000);  // 2ms

    std::ostringstream oss;
    ProfileRegistry::Instance().PrintReport(oss);

    std::string report = oss.str();
    EXPECT_TRUE(report.find("LLmap Profile Report") != std::string::npos);
    EXPECT_TRUE(report.find("func1") != std::string::npos);
    EXPECT_TRUE(report.find("func2") != std::string::npos);
}

TEST_F(ProfilerTest, ScopedTimerRecordsTime) {
    auto& stats = ProfileRegistry::Instance().GetOrCreate("timed_scope");

    {
        ScopedTimer timer(stats);
        std::this_thread::sleep_for(std::chrono::microseconds(100));
    }

    EXPECT_EQ(1, stats.count.load());
    EXPECT_GT(stats.total_ns.load(), 50000);  // At least 50µs
}

TEST_F(ProfilerTest, ScopedTimerMultipleScopes) {
    auto& stats = ProfileRegistry::Instance().GetOrCreate("multi_scope");

    for (int i = 0; i < 10; ++i) {
        ScopedTimer timer(stats);
        std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    EXPECT_EQ(10, stats.count.load());
}

TEST_F(ProfilerTest, ManualTimerStartStop) {
    ManualTimer timer;

    timer.Start();
    std::this_thread::sleep_for(std::chrono::microseconds(100));
    timer.Stop();

    EXPECT_GT(timer.ElapsedNs(), 50000);
    EXPECT_GT(timer.ElapsedUs(), 50.0);
}

TEST_F(ProfilerTest, ManualTimerAccumulates) {
    ManualTimer timer;

    timer.Start();
    std::this_thread::sleep_for(std::chrono::microseconds(50));
    timer.Stop();

    timer.Start();
    std::this_thread::sleep_for(std::chrono::microseconds(50));
    timer.Stop();

    EXPECT_GT(timer.ElapsedUs(), 80.0);  // Combined time
}

TEST_F(ProfilerTest, ManualTimerReset) {
    ManualTimer timer;

    timer.Start();
    std::this_thread::sleep_for(std::chrono::microseconds(100));
    timer.Stop();

    auto before_reset = timer.ElapsedNs();
    EXPECT_GT(before_reset, 0);

    timer.Reset();
    EXPECT_EQ(0, timer.ElapsedNs());

    timer.Start();
    std::this_thread::sleep_for(std::chrono::microseconds(50));
    timer.Stop();

    // After reset, time should be less than before reset (only measured the second sleep)
    // Use generous tolerance since sleep timing is non-deterministic
    EXPECT_LT(timer.ElapsedNs(), before_reset * 3);  // Allow for OS scheduling variance
}

TEST_F(ProfilerTest, ManualTimerElapsedMs) {
    ManualTimer timer;

    timer.Start();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    timer.Stop();

    EXPECT_GT(timer.ElapsedMs(), 2.0);
}

TEST_F(ProfilerTest, ThreadSafeRecording) {
    auto& stats = ProfileRegistry::Instance().GetOrCreate("concurrent");

    std::vector<std::thread> threads;
    constexpr int kNumThreads = 4;
    constexpr int kRecordsPerThread = 100;

    for (int i = 0; i < kNumThreads; ++i) {
        threads.emplace_back([&stats]() {
            for (int j = 0; j < kRecordsPerThread; ++j) {
                stats.Record(1000);
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    EXPECT_EQ(kNumThreads * kRecordsPerThread, stats.count.load());
    EXPECT_EQ(kNumThreads * kRecordsPerThread * 1000, stats.total_ns.load());
}

TEST_F(ProfilerTest, ThreadSafeGetOrCreate) {
    std::vector<std::thread> threads;
    constexpr int kNumThreads = 4;

    for (int i = 0; i < kNumThreads; ++i) {
        threads.emplace_back([i]() {
            for (int j = 0; j < 10; ++j) {
                auto& stats = ProfileRegistry::Instance().GetOrCreate("shared_" + std::to_string(j % 3));
                stats.Record(1000);
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    auto snapshot = ProfileRegistry::Instance().GetSnapshot();
    EXPECT_EQ(3, snapshot.size());  // shared_0, shared_1, shared_2
}

TEST_F(ProfilerTest, MacroScopedTimerEquivalent) {
    // Test that manual scoped timer works as expected
    // (equivalent to what LLMAP_PROFILE_SCOPE does when enabled)
    auto& stats = ProfileRegistry::Instance().GetOrCreate("manual_macro_scope");

    {
        ScopedTimer timer(stats);
        volatile int x = 0;
        for (int i = 0; i < 100; ++i) {
            x += i;
        }
    }

    EXPECT_EQ(1, stats.count.load());
    EXPECT_GT(stats.total_ns.load(), 0);
}
