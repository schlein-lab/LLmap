// LLmap — Profiler: lightweight profiling infrastructure for hot paths.
//
// Provides:
//   - LLMAP_PROFILE_SCOPE(name) — Scoped timer for a block
//   - LLMAP_PROFILE_FUNCTION()  — Scoped timer using function name
//   - ProfileRegistry — Singleton to collect and report timings
//   - ScopedTimer — RAII timer that records to registry
//
// Thread-safe by default. Zero overhead when LLMAP_ENABLE_PROFILING is off.

#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::core {

// Statistics for a single profiled scope
struct ProfileStats {
    std::string name;
    std::atomic<uint64_t> count{0};
    std::atomic<uint64_t> total_ns{0};
    std::atomic<uint64_t> min_ns{UINT64_MAX};
    std::atomic<uint64_t> max_ns{0};

    ProfileStats() = default;
    explicit ProfileStats(std::string_view n) : name(n) {}

    // Copy/move for vector operations (needed for sorting/snapshot)
    ProfileStats(const ProfileStats& other)
        : name(other.name),
          count(other.count.load()),
          total_ns(other.total_ns.load()),
          min_ns(other.min_ns.load()),
          max_ns(other.max_ns.load()) {}

    ProfileStats(ProfileStats&& other) noexcept
        : name(std::move(other.name)),
          count(other.count.load()),
          total_ns(other.total_ns.load()),
          min_ns(other.min_ns.load()),
          max_ns(other.max_ns.load()) {}

    ProfileStats& operator=(const ProfileStats& other) {
        if (this != &other) {
            name = other.name;
            count.store(other.count.load());
            total_ns.store(other.total_ns.load());
            min_ns.store(other.min_ns.load());
            max_ns.store(other.max_ns.load());
        }
        return *this;
    }

    ProfileStats& operator=(ProfileStats&& other) noexcept {
        if (this != &other) {
            name = std::move(other.name);
            count.store(other.count.load());
            total_ns.store(other.total_ns.load());
            min_ns.store(other.min_ns.load());
            max_ns.store(other.max_ns.load());
        }
        return *this;
    }

    void Record(uint64_t ns) {
        count.fetch_add(1, std::memory_order_relaxed);
        total_ns.fetch_add(ns, std::memory_order_relaxed);

        uint64_t old_min = min_ns.load(std::memory_order_relaxed);
        while (ns < old_min && !min_ns.compare_exchange_weak(
            old_min, ns, std::memory_order_relaxed)) {}

        uint64_t old_max = max_ns.load(std::memory_order_relaxed);
        while (ns > old_max && !max_ns.compare_exchange_weak(
            old_max, ns, std::memory_order_relaxed)) {}
    }

    double AverageUs() const {
        uint64_t c = count.load();
        if (c == 0) return 0.0;
        return static_cast<double>(total_ns.load()) / c / 1000.0;
    }

    double TotalMs() const {
        return static_cast<double>(total_ns.load()) / 1e6;
    }
};

// Global registry for profile statistics
class ProfileRegistry {
public:
    static ProfileRegistry& Instance() {
        static ProfileRegistry instance;
        return instance;
    }

    ProfileStats& GetOrCreate(std::string_view name) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = stats_.find(std::string(name));
        if (it == stats_.end()) {
            auto [iter, _] = stats_.emplace(std::string(name), ProfileStats(name));
            return iter->second;
        }
        return it->second;
    }

    void Reset() {
        std::lock_guard<std::mutex> lock(mutex_);
        stats_.clear();
    }

    // Get a snapshot of all stats (sorted by total time descending)
    std::vector<ProfileStats> GetSnapshot() const {
        std::lock_guard<std::mutex> lock(mutex_);
        std::vector<ProfileStats> result;
        result.reserve(stats_.size());
        for (const auto& [_, stats] : stats_) {
            result.push_back(stats);
        }
        std::sort(result.begin(), result.end(),
            [](const ProfileStats& a, const ProfileStats& b) {
                return a.total_ns.load() > b.total_ns.load();
            });
        return result;
    }

    // Print report to stream
    void PrintReport(std::ostream& os) const {
        auto snapshot = GetSnapshot();
        os << "\n=== LLmap Profile Report ===\n";
        os << "Scope                                  | Count      | Total (ms) | Avg (µs)   | Min (µs)   | Max (µs)\n";
        os << "---------------------------------------|------------|------------|------------|------------|------------\n";

        for (const auto& s : snapshot) {
            char buf[256];
            uint64_t cnt = s.count.load();
            double total_ms = s.TotalMs();
            double avg_us = s.AverageUs();
            double min_us = s.min_ns.load() / 1000.0;
            double max_us = s.max_ns.load() / 1000.0;

            if (cnt == 0) {
                min_us = 0.0;
                max_us = 0.0;
            }

            snprintf(buf, sizeof(buf), "%-38s | %10lu | %10.2f | %10.2f | %10.2f | %10.2f",
                     s.name.substr(0, 38).c_str(),
                     static_cast<unsigned long>(cnt),
                     total_ms, avg_us, min_us, max_us);
            os << buf << "\n";
        }
        os << "============================\n\n";
    }

    // Iterate over stats
    template<typename Func>
    void ForEach(Func&& func) const {
        std::lock_guard<std::mutex> lock(mutex_);
        for (const auto& [name, stats] : stats_) {
            func(name, stats);
        }
    }

private:
    ProfileRegistry() = default;
    mutable std::mutex mutex_;
    std::map<std::string, ProfileStats> stats_;
};

// RAII scoped timer
class ScopedTimer {
public:
    explicit ScopedTimer(ProfileStats& stats)
        : stats_(stats),
          start_(std::chrono::high_resolution_clock::now()) {}

    ~ScopedTimer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end - start_).count();
        stats_.Record(static_cast<uint64_t>(ns));
    }

    // Non-copyable, non-movable
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&) = delete;
    ScopedTimer& operator=(ScopedTimer&&) = delete;

private:
    ProfileStats& stats_;
    std::chrono::high_resolution_clock::time_point start_;
};

// Manual timer for spans that cross scope boundaries
class ManualTimer {
public:
    ManualTimer() : running_(false), elapsed_ns_(0) {}

    void Start() {
        start_ = std::chrono::high_resolution_clock::now();
        running_ = true;
    }

    void Stop() {
        if (running_) {
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start_).count();
            running_ = false;
        }
    }

    void Reset() {
        elapsed_ns_ = 0;
        running_ = false;
    }

    uint64_t ElapsedNs() const { return elapsed_ns_; }
    double ElapsedUs() const { return elapsed_ns_ / 1000.0; }
    double ElapsedMs() const { return elapsed_ns_ / 1e6; }

private:
    std::chrono::high_resolution_clock::time_point start_;
    bool running_;
    uint64_t elapsed_ns_;
};

}  // namespace llmap::core

// Profiling macros
#ifdef LLMAP_ENABLE_PROFILING

#define LLMAP_PROFILE_SCOPE(name) \
    static auto& LLMAP_PROFILE_STATS_##__LINE__ = \
        ::llmap::core::ProfileRegistry::Instance().GetOrCreate(name); \
    ::llmap::core::ScopedTimer LLMAP_PROFILE_TIMER_##__LINE__(LLMAP_PROFILE_STATS_##__LINE__)

#define LLMAP_PROFILE_FUNCTION() \
    LLMAP_PROFILE_SCOPE(__FUNCTION__)

#define LLMAP_PROFILE_REPORT(os) \
    ::llmap::core::ProfileRegistry::Instance().PrintReport(os)

#define LLMAP_PROFILE_RESET() \
    ::llmap::core::ProfileRegistry::Instance().Reset()

#else  // LLMAP_ENABLE_PROFILING not defined

#define LLMAP_PROFILE_SCOPE(name) ((void)0)
#define LLMAP_PROFILE_FUNCTION() ((void)0)
#define LLMAP_PROFILE_REPORT(os) ((void)0)
#define LLMAP_PROFILE_RESET() ((void)0)

#endif  // LLMAP_ENABLE_PROFILING
