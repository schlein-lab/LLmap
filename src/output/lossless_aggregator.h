// LLmap — Lossless invariant accountant.
//
// Single point of truth for "are we still lossless?" across the entire
// pipeline. Every AlignmentRecord emitted by the alignment pipeline
// flows through Observe(); at the end of the run Finalize() writes a
// summary JSON that downstream tools (and the test harness) can read to
// verify the contract:
//
//     n_records_emitted == n_input_reads
//     all records satisfy is_lossless_consistent()
//
// If either invariant is violated the run is poisoned — the aggregator
// flips a flag and the JSON is written with `"lossless_invariant_ok":
// false` plus a list of offending read_ids (capped to avoid runaway
// memory if everything is broken).
//
// Cost: O(1) per record (counter increment); thread-safe via internal
// mutex so the pipeline can call it from worker threads.
//
// Why a separate module rather than counters inside BamWriter? Because
// the lossless contract is independent of any specific output sink.
// Parquet, BAM, GAF, and the eventual streaming-API output all want the
// same accountancy.

#pragma once

#include "core/alignment_record.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace llmap::output {

/// Snapshot of counters; cheap to copy.
struct LosslessCounters {
    std::uint64_t n_input_reads_declared{0};       // what the run expected
    std::uint64_t n_records_emitted{0};            // what we actually saw

    // Per-status counters; index by AlignmentStatus integer value.
    std::array<std::uint64_t, 10> by_status{};

    // Per-rejection-reason counters; only populated for Unmapped/Tentative.
    std::array<std::uint64_t, 8> by_rejection{};

    // Per-TranscriptKind counters (uint16_t enum; we keep a small fixed
    // table covering the documented kinds and lump the rest under
    // 'other'). Index 0 is Unknown.
    std::unordered_map<std::uint16_t, std::uint64_t> by_kind;

    // Lossless invariant tracking
    bool invariant_ok{true};
    std::vector<std::string> offending_read_ids;   // capped at kMaxOffenders
};

class LosslessAggregator {
public:
    /// Cap on remembered offender IDs. Beyond this we just count.
    static constexpr std::size_t kMaxOffenders = 100;

    /// Tell the aggregator how many input reads the run will process.
    /// Allows the final `n_records == n_input` check to be meaningful.
    /// Optional — pass 0 if streaming and the input count is unknown.
    void SetExpectedInputCount(std::uint64_t n) noexcept;

    /// Record one emitted AlignmentRecord. Thread-safe.
    /// Returns the result of is_lossless_consistent() for the record,
    /// so callers can early-out if they wish (the aggregator already
    /// flags it internally either way).
    bool Observe(const AlignmentRecord& rec);

    /// Immutable view of current counters. Atomic snapshot.
    [[nodiscard]] LosslessCounters Snapshot() const;

    /// Write the summary JSON to `path`. Returns true on success.
    /// Format documented in src/output/lossless_aggregator.cpp.
    bool WriteSummary(const std::filesystem::path& path) const;

    /// Reset all counters and invariant state. Mainly for tests.
    void Reset();

private:
    mutable std::mutex mu_;
    LosslessCounters c_;
};

}  // namespace llmap::output
