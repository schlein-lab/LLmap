// LLmap — ThreadPool unit tests

#include "core/thread_pool.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

namespace llmap::core {
namespace {

class ThreadPoolTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Use a fixed number of threads for reproducible tests
        pool_ = std::make_unique<ThreadPool>(4);
    }

    void TearDown() override {
        pool_.reset();
    }

    std::unique_ptr<ThreadPool> pool_;
};

// Basic construction tests

TEST_F(ThreadPoolTest, ConstructsWithDefaultThreads) {
    ThreadPool pool;
    EXPECT_GT(pool.NumThreads(), 0);
}

TEST_F(ThreadPoolTest, ConstructsWithSpecifiedThreads) {
    ThreadPool pool(8);
    EXPECT_EQ(pool.NumThreads(), 8);
}

TEST_F(ThreadPoolTest, ZeroThreadsBecomesOne) {
    ThreadPool pool(0);
    EXPECT_EQ(pool.NumThreads(), 1);
}

TEST_F(ThreadPoolTest, StartsIdle) {
    EXPECT_TRUE(pool_->IsIdle());
}

// Basic task submission

TEST_F(ThreadPoolTest, SubmitReturnsCorrectResult) {
    auto future = pool_->Submit([]() { return 42; });
    EXPECT_EQ(future.get(), 42);
}

TEST_F(ThreadPoolTest, SubmitWithArgs) {
    auto future = pool_->Submit([](int a, int b) { return a + b; }, 10, 20);
    EXPECT_EQ(future.get(), 30);
}

TEST_F(ThreadPoolTest, ExecuteRunsTask) {
    std::atomic<bool> ran{false};
    pool_->Execute([&ran]() { ran = true; });
    pool_->WaitAll();
    EXPECT_TRUE(ran.load());
}

TEST_F(ThreadPoolTest, MultipleTasksComplete) {
    std::vector<std::future<int>> futures;
    for (int i = 0; i < 100; ++i) {
        futures.push_back(pool_->Submit([i]() { return i * 2; }));
    }

    for (int i = 0; i < 100; ++i) {
        EXPECT_EQ(futures[i].get(), i * 2);
    }
}

// Thread safety tests

TEST_F(ThreadPoolTest, ConcurrentSubmits) {
    std::atomic<int> counter{0};
    std::vector<std::thread> submitters;

    for (int t = 0; t < 4; ++t) {
        submitters.emplace_back([this, &counter]() {
            for (int i = 0; i < 250; ++i) {
                pool_->Execute([&counter]() {
                    counter.fetch_add(1, std::memory_order_relaxed);
                });
            }
        });
    }

    for (auto& t : submitters) {
        t.join();
    }

    pool_->WaitAll();
    EXPECT_EQ(counter.load(), 1000);
}

TEST_F(ThreadPoolTest, TasksCanEnqueueMoreTasks) {
    std::atomic<int> counter{0};

    pool_->Execute([this, &counter]() {
        for (int i = 0; i < 10; ++i) {
            pool_->Execute([&counter]() {
                counter.fetch_add(1, std::memory_order_relaxed);
            });
        }
    });

    pool_->WaitAll();
    EXPECT_EQ(counter.load(), 10);
}

// WaitAll tests

TEST_F(ThreadPoolTest, WaitAllBlocksUntilComplete) {
    std::atomic<int> completed{0};

    for (int i = 0; i < 50; ++i) {
        pool_->Execute([&completed]() {
            std::this_thread::sleep_for(std::chrono::microseconds(100));
            completed.fetch_add(1, std::memory_order_relaxed);
        });
    }

    pool_->WaitAll();
    EXPECT_EQ(completed.load(), 50);
}

TEST_F(ThreadPoolTest, WaitAllOnEmptyPoolReturnsImmediately) {
    auto start = std::chrono::high_resolution_clock::now();
    pool_->WaitAll();
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    EXPECT_LT(duration.count(), 10);
}

// Statistics tests

TEST_F(ThreadPoolTest, StatsTrackSubmitted) {
    for (int i = 0; i < 10; ++i) {
        pool_->Execute([]() {});
    }
    pool_->WaitAll();

    EXPECT_EQ(pool_->Stats().tasks_submitted.load(), 10);
}

TEST_F(ThreadPoolTest, StatsTrackCompleted) {
    for (int i = 0; i < 10; ++i) {
        pool_->Execute([]() {});
    }
    pool_->WaitAll();

    EXPECT_EQ(pool_->Stats().tasks_completed.load(), 10);
}

TEST_F(ThreadPoolTest, ResetStatsClears) {
    pool_->Execute([]() {});
    pool_->WaitAll();

    EXPECT_GT(pool_->Stats().tasks_submitted.load(), 0);

    pool_->ResetStats();

    EXPECT_EQ(pool_->Stats().tasks_submitted.load(), 0);
    EXPECT_EQ(pool_->Stats().tasks_completed.load(), 0);
}

// ParallelFor tests

TEST_F(ThreadPoolTest, ParallelForEmpty) {
    std::atomic<int> counter{0};
    ParallelFor(*pool_, 0, [&counter](size_t) {
        counter.fetch_add(1, std::memory_order_relaxed);
    });
    EXPECT_EQ(counter.load(), 0);
}

TEST_F(ThreadPoolTest, ParallelForSingleElement) {
    std::atomic<int> counter{0};
    ParallelFor(*pool_, 1, [&counter](size_t i) {
        EXPECT_EQ(i, 0);
        counter.fetch_add(1, std::memory_order_relaxed);
    });
    EXPECT_EQ(counter.load(), 1);
}

TEST_F(ThreadPoolTest, ParallelForAllIndices) {
    constexpr size_t N = 1000;
    std::vector<std::atomic<bool>> visited(N);
    for (auto& v : visited) v = false;

    ParallelFor(*pool_, N, [&visited](size_t i) {
        visited[i].store(true, std::memory_order_relaxed);
    });

    for (size_t i = 0; i < N; ++i) {
        EXPECT_TRUE(visited[i].load()) << "Index " << i << " not visited";
    }
}

TEST_F(ThreadPoolTest, ParallelForWithChunkSize) {
    constexpr size_t N = 100;
    std::atomic<int> counter{0};

    ParallelFor(*pool_, N, [&counter](size_t) {
        counter.fetch_add(1, std::memory_order_relaxed);
    }, 10);  // Chunk size of 10

    EXPECT_EQ(counter.load(), 100);
}

TEST_F(ThreadPoolTest, ParallelForCompute) {
    constexpr size_t N = 10000;
    std::vector<int> data(N);
    std::iota(data.begin(), data.end(), 0);

    std::vector<int> results(N);

    ParallelFor(*pool_, N, [&data, &results](size_t i) {
        results[i] = data[i] * data[i];
    });

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(results[i], static_cast<int>(i * i));
    }
}

TEST_F(ThreadPoolTest, ParallelForDefaultPool) {
    constexpr size_t N = 100;
    std::atomic<int> counter{0};

    ParallelFor(N, [&counter](size_t) {
        counter.fetch_add(1, std::memory_order_relaxed);
    });

    EXPECT_EQ(counter.load(), 100);
}

// ParallelMap tests

TEST_F(ThreadPoolTest, ParallelMapEmpty) {
    std::vector<int> inputs;
    auto results = ParallelMap(*pool_, inputs, [](int x) { return x * 2; });
    EXPECT_TRUE(results.empty());
}

TEST_F(ThreadPoolTest, ParallelMapTransforms) {
    std::vector<int> inputs = {1, 2, 3, 4, 5};
    auto results = ParallelMap(*pool_, inputs, [](int x) { return x * x; });

    EXPECT_EQ(results.size(), 5);
    EXPECT_EQ(results[0], 1);
    EXPECT_EQ(results[1], 4);
    EXPECT_EQ(results[2], 9);
    EXPECT_EQ(results[3], 16);
    EXPECT_EQ(results[4], 25);
}

TEST_F(ThreadPoolTest, ParallelMapPreservesOrder) {
    constexpr size_t N = 1000;
    std::vector<size_t> inputs(N);
    std::iota(inputs.begin(), inputs.end(), 0);

    auto results = ParallelMap(*pool_, inputs, [](size_t x) { return x * 3; });

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(results[i], i * 3);
    }
}

TEST_F(ThreadPoolTest, ParallelMapStrings) {
    std::vector<std::string> inputs = {"hello", "world", "test"};
    auto results = ParallelMap(*pool_, inputs, [](const std::string& s) {
        return s.size();
    });

    EXPECT_EQ(results.size(), 3);
    EXPECT_EQ(results[0], 5);
    EXPECT_EQ(results[1], 5);
    EXPECT_EQ(results[2], 4);
}

// BatchProcessor tests

TEST_F(ThreadPoolTest, BatchProcessorBasic) {
    BatchProcessor<int, int> processor(*pool_);
    processor.SetProcessor([](int x) { return x * 2; });

    std::vector<int> inputs = {1, 2, 3, 4, 5};
    auto results = processor.Process(inputs);

    EXPECT_EQ(results.size(), 5);
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_EQ(results[i], inputs[i] * 2);
    }
}

TEST_F(ThreadPoolTest, BatchProcessorProgress) {
    BatchProcessor<int, int> processor(*pool_);

    std::atomic<size_t> last_completed{0};
    processor.SetProcessor([](int x) { return x; })
             .SetProgress([&last_completed](size_t completed, size_t /*total*/) {
                 last_completed.store(completed, std::memory_order_relaxed);
             }, 10);

    std::vector<int> inputs(100);
    std::iota(inputs.begin(), inputs.end(), 0);

    processor.Process(inputs);

    EXPECT_EQ(last_completed.load(), 100);
}

TEST_F(ThreadPoolTest, BatchProcessorChunkSize) {
    BatchProcessor<int, int> processor(*pool_);
    processor.SetProcessor([](int x) { return x; })
             .SetChunkSize(5);

    std::vector<int> inputs(50);
    std::iota(inputs.begin(), inputs.end(), 0);

    auto results = processor.Process(inputs);

    EXPECT_EQ(results.size(), 50);
    for (size_t i = 0; i < 50; ++i) {
        EXPECT_EQ(results[i], static_cast<int>(i));
    }
}

// Default pool tests

TEST(ThreadPoolDefaultTest, DefaultPoolExists) {
    ThreadPool& pool = ThreadPool::Default();
    EXPECT_GT(pool.NumThreads(), 0);
}

TEST(ThreadPoolDefaultTest, DefaultPoolSameInstance) {
    ThreadPool& pool1 = ThreadPool::Default();
    ThreadPool& pool2 = ThreadPool::Default();
    EXPECT_EQ(&pool1, &pool2);
}

TEST(ThreadPoolDefaultTest, DefaultPoolWorks) {
    auto future = ThreadPool::Default().Submit([]() { return 123; });
    EXPECT_EQ(future.get(), 123);
}

// Performance sanity checks

TEST_F(ThreadPoolTest, ParallelFasterThanSequential) {
    constexpr size_t N = 10000;
    std::vector<double> data(N);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (auto& x : data) x = dist(rng);

    auto compute = [](double x) {
        // Some computation
        double result = x;
        for (int i = 0; i < 100; ++i) {
            result = std::sin(result) + std::cos(result);
        }
        return result;
    };

    // Sequential
    std::vector<double> seq_results(N);
    auto seq_start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; ++i) {
        seq_results[i] = compute(data[i]);
    }
    auto seq_end = std::chrono::high_resolution_clock::now();
    auto seq_time = std::chrono::duration_cast<std::chrono::microseconds>(seq_end - seq_start).count();

    // Parallel
    std::vector<double> par_results(N);
    auto par_start = std::chrono::high_resolution_clock::now();
    ParallelFor(*pool_, N, [&data, &par_results, &compute](size_t i) {
        par_results[i] = compute(data[i]);
    });
    auto par_end = std::chrono::high_resolution_clock::now();
    auto par_time = std::chrono::duration_cast<std::chrono::microseconds>(par_end - par_start).count();

    // Verify correctness
    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(seq_results[i], par_results[i]);
    }

    // Parallel should generally be faster (with 4 threads, expect >1.5x speedup)
    // But on loaded systems, thread scheduling can cause variance
    // Use a lenient check: parallel shouldn't be more than 2x slower
    if (pool_->NumThreads() > 1) {
        EXPECT_LT(par_time, seq_time * 2) << "Parallel (" << par_time << "us) should not be more than 2x slower than sequential (" << seq_time << "us)";
    }
}

// Edge cases

TEST_F(ThreadPoolTest, LargeNumberOfSmallTasks) {
    std::atomic<int> counter{0};

    for (int i = 0; i < 10000; ++i) {
        pool_->Execute([&counter]() {
            counter.fetch_add(1, std::memory_order_relaxed);
        });
    }

    pool_->WaitAll();
    EXPECT_EQ(counter.load(), 10000);
}

TEST_F(ThreadPoolTest, VaryingTaskDurations) {
    std::atomic<int> counter{0};
    std::mt19937 rng(42);

    for (int i = 0; i < 100; ++i) {
        int delay = rng() % 1000;  // 0-999 microseconds
        pool_->Execute([&counter, delay]() {
            std::this_thread::sleep_for(std::chrono::microseconds(delay));
            counter.fetch_add(1, std::memory_order_relaxed);
        });
    }

    pool_->WaitAll();
    EXPECT_EQ(counter.load(), 100);
}

}  // namespace
}  // namespace llmap::core
