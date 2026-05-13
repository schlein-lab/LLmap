// LLmap — ThreadPool: work-stealing thread pool for parallel batch processing.
//
// Provides:
//   - ThreadPool — Fixed-size pool with work-stealing for load balancing
//   - ParallelFor — Parallel for-loop with automatic chunking
//   - ParallelMap — Transform inputs to outputs in parallel
//
// Thread safety: ThreadPool is fully thread-safe. Tasks can enqueue more tasks.
// Designed for batch sequence processing where work per item varies.

#pragma once

#include <atomic>
#include <cassert>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <thread>
#include <type_traits>
#include <vector>

namespace llmap::core {

// Statistics for thread pool operations
struct ThreadPoolStats {
    std::atomic<uint64_t> tasks_submitted{0};
    std::atomic<uint64_t> tasks_completed{0};
    std::atomic<uint64_t> tasks_stolen{0};
    std::atomic<uint64_t> total_wait_ns{0};

    ThreadPoolStats() = default;

    ThreadPoolStats(const ThreadPoolStats& other)
        : tasks_submitted(other.tasks_submitted.load()),
          tasks_completed(other.tasks_completed.load()),
          tasks_stolen(other.tasks_stolen.load()),
          total_wait_ns(other.total_wait_ns.load()) {}

    void Reset() {
        tasks_submitted.store(0);
        tasks_completed.store(0);
        tasks_stolen.store(0);
        total_wait_ns.store(0);
    }
};

// Work-stealing thread pool optimized for batch processing
class ThreadPool {
public:
    // Create pool with hardware_concurrency threads by default
    ThreadPool();

    // Create pool with specified number of threads
    explicit ThreadPool(size_t num_threads);

    ~ThreadPool();

    // Non-copyable, non-movable (owns threads)
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool& operator=(ThreadPool&&) = delete;

    // Submit a task and get a future for the result
    template <typename F, typename... Args>
    auto Submit(F&& func, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>> {
        using ReturnType = std::invoke_result_t<F, Args...>;

        auto task = std::make_shared<std::packaged_task<ReturnType()>>(
            std::bind(std::forward<F>(func), std::forward<Args>(args)...)
        );

        std::future<ReturnType> result = task->get_future();

        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            if (shutdown_) {
                throw std::runtime_error("ThreadPool: cannot submit to stopped pool");
            }
            tasks_.emplace_back([task]() { (*task)(); });
            stats_.tasks_submitted.fetch_add(1, std::memory_order_relaxed);
        }

        cv_.notify_one();
        return result;
    }

    // Submit a task without caring about the result
    template <typename F>
    void Execute(F&& func) {
        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            if (shutdown_) {
                throw std::runtime_error("ThreadPool: cannot submit to stopped pool");
            }
            tasks_.emplace_back(std::forward<F>(func));
            stats_.tasks_submitted.fetch_add(1, std::memory_order_relaxed);
        }
        cv_.notify_one();
    }

    // Wait for all submitted tasks to complete
    void WaitAll();

    // Get number of worker threads
    size_t NumThreads() const { return threads_.size(); }

    // Get statistics
    const ThreadPoolStats& Stats() const { return stats_; }

    // Reset statistics
    void ResetStats() { stats_.Reset(); }

    // Check if pool is idle (no pending tasks)
    bool IsIdle() const {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        return tasks_.empty() && active_tasks_ == 0;
    }

    // Get global default thread pool (created lazily)
    static ThreadPool& Default();

private:
    void WorkerLoop(size_t worker_id);
    std::optional<std::function<void()>> TrySteal();

    std::vector<std::thread> threads_;
    mutable std::mutex queue_mutex_;
    std::condition_variable cv_;
    std::condition_variable done_cv_;
    std::deque<std::function<void()>> tasks_;
    std::atomic<bool> shutdown_{false};
    size_t active_tasks_ = 0;
    ThreadPoolStats stats_;
};

// Parallel for-loop: execute func(i) for i in [0, count)
//
// Parameters:
//   pool: Thread pool to use
//   count: Number of iterations
//   func: Function taking size_t index
//   chunk_size: Minimum iterations per task (default: auto-computed)
//
// Blocks until all iterations complete.
template <typename F>
void ParallelFor(ThreadPool& pool, size_t count, F&& func, size_t chunk_size = 0) {
    if (count == 0) return;

    size_t num_threads = pool.NumThreads();
    if (num_threads == 0) num_threads = 1;

    // Auto chunk size: balance overhead vs load balancing
    if (chunk_size == 0) {
        chunk_size = std::max<size_t>(1, count / (num_threads * 4));
    }

    std::atomic<size_t> next_chunk{0};
    std::atomic<size_t> tasks_done{0};
    size_t num_chunks = (count + chunk_size - 1) / chunk_size;

    std::mutex done_mutex;
    std::condition_variable done_cv;

    auto worker = [&, func = std::forward<F>(func)]() mutable {
        while (true) {
            size_t chunk_idx = next_chunk.fetch_add(1, std::memory_order_relaxed);
            if (chunk_idx >= num_chunks) break;

            size_t start = chunk_idx * chunk_size;
            size_t end = std::min(start + chunk_size, count);

            for (size_t i = start; i < end; ++i) {
                func(i);
            }
        }

        if (tasks_done.fetch_add(1, std::memory_order_acq_rel) + 1 == num_threads) {
            std::lock_guard<std::mutex> lock(done_mutex);
            done_cv.notify_all();
        }
    };

    // Submit workers (one per thread in pool)
    for (size_t t = 0; t < num_threads - 1; ++t) {
        pool.Execute(worker);
    }

    // Use calling thread as worker too
    worker();

    // Wait for all workers to finish
    std::unique_lock<std::mutex> lock(done_mutex);
    done_cv.wait(lock, [&] {
        return tasks_done.load(std::memory_order_acquire) >= num_threads;
    });
}

// Parallel for using default pool
template <typename F>
void ParallelFor(size_t count, F&& func, size_t chunk_size = 0) {
    ParallelFor(ThreadPool::Default(), count, std::forward<F>(func), chunk_size);
}

// Parallel map: transform inputs to outputs in parallel
//
// Parameters:
//   pool: Thread pool to use
//   inputs: Input range
//   func: Function taking input element, returning output
//
// Returns: vector of outputs in same order as inputs
template <typename InputRange, typename F>
auto ParallelMap(ThreadPool& pool, const InputRange& inputs, F&& func)
    -> std::vector<std::invoke_result_t<F, decltype(*std::begin(inputs))>> {
    using InputType = decltype(*std::begin(inputs));
    using OutputType = std::invoke_result_t<F, InputType>;

    size_t count = std::distance(std::begin(inputs), std::end(inputs));
    std::vector<OutputType> results(count);

    ParallelFor(pool, count, [&, func = std::forward<F>(func)](size_t i) mutable {
        auto input_it = std::begin(inputs);
        std::advance(input_it, i);
        results[i] = func(*input_it);
    });

    return results;
}

// Parallel map using default pool
template <typename InputRange, typename F>
auto ParallelMap(const InputRange& inputs, F&& func)
    -> std::vector<std::invoke_result_t<F, decltype(*std::begin(inputs))>> {
    return ParallelMap(ThreadPool::Default(), inputs, std::forward<F>(func));
}

// Parallel batch processor for sequences
//
// Provides optimized batch processing for read alignment workloads:
// - Automatic chunking based on sequence lengths
// - Per-thread scratch space
// - Progress callbacks
template <typename Input, typename Output>
class BatchProcessor {
public:
    using ProcessFunc = std::function<Output(const Input&)>;
    using ProgressFunc = std::function<void(size_t completed, size_t total)>;

    explicit BatchProcessor(ThreadPool& pool) : pool_(pool) {}

    // Set the processing function
    BatchProcessor& SetProcessor(ProcessFunc func) {
        processor_ = std::move(func);
        return *this;
    }

    // Set progress callback (called periodically, not per-item)
    BatchProcessor& SetProgress(ProgressFunc func, size_t interval = 1000) {
        progress_ = std::move(func);
        progress_interval_ = interval;
        return *this;
    }

    // Set chunk size (0 = auto)
    BatchProcessor& SetChunkSize(size_t chunk_size) {
        chunk_size_ = chunk_size;
        return *this;
    }

    // Process a batch of inputs
    std::vector<Output> Process(std::span<const Input> inputs) {
        assert(processor_ && "Must set processor before calling Process");

        size_t count = inputs.size();
        if (count == 0) return {};

        std::vector<Output> results(count);
        std::atomic<size_t> completed{0};

        ParallelFor(pool_, count, [&](size_t i) {
            results[i] = processor_(inputs[i]);

            size_t done = completed.fetch_add(1, std::memory_order_relaxed) + 1;
            if (progress_ && done % progress_interval_ == 0) {
                progress_(done, count);
            }
        }, chunk_size_);

        if (progress_) {
            progress_(count, count);
        }

        return results;
    }

private:
    ThreadPool& pool_;
    ProcessFunc processor_;
    ProgressFunc progress_;
    size_t progress_interval_ = 1000;
    size_t chunk_size_ = 0;
};

}  // namespace llmap::core
