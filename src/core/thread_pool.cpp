// LLmap — ThreadPool implementation

#include "core/thread_pool.h"

#include <chrono>

namespace llmap::core {

ThreadPool::ThreadPool() : ThreadPool(std::thread::hardware_concurrency()) {}

ThreadPool::ThreadPool(size_t num_threads) {
    if (num_threads == 0) {
        num_threads = 1;
    }

    threads_.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        threads_.emplace_back(&ThreadPool::WorkerLoop, this, i);
    }
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        shutdown_ = true;
    }
    cv_.notify_all();

    for (auto& thread : threads_) {
        if (thread.joinable()) {
            thread.join();
        }
    }
}

void ThreadPool::WorkerLoop(size_t /*worker_id*/) {
    while (true) {
        std::function<void()> task;

        {
            std::unique_lock<std::mutex> lock(queue_mutex_);

            auto wait_start = std::chrono::high_resolution_clock::now();

            cv_.wait(lock, [this] {
                return shutdown_ || !tasks_.empty();
            });

            auto wait_end = std::chrono::high_resolution_clock::now();
            auto wait_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                wait_end - wait_start).count();
            stats_.total_wait_ns.fetch_add(wait_ns, std::memory_order_relaxed);

            if (shutdown_ && tasks_.empty()) {
                return;
            }

            if (tasks_.empty()) {
                continue;
            }

            task = std::move(tasks_.front());
            tasks_.pop_front();
            ++active_tasks_;
        }

        task();

        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            --active_tasks_;
            stats_.tasks_completed.fetch_add(1, std::memory_order_relaxed);
        }
        done_cv_.notify_all();
    }
}

void ThreadPool::WaitAll() {
    std::unique_lock<std::mutex> lock(queue_mutex_);
    done_cv_.wait(lock, [this] {
        return tasks_.empty() && active_tasks_ == 0;
    });
}

std::optional<std::function<void()>> ThreadPool::TrySteal() {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    if (tasks_.empty()) {
        return std::nullopt;
    }

    auto task = std::move(tasks_.back());
    tasks_.pop_back();
    stats_.tasks_stolen.fetch_add(1, std::memory_order_relaxed);
    return task;
}

ThreadPool& ThreadPool::Default() {
    static ThreadPool instance;
    return instance;
}

}  // namespace llmap::core
