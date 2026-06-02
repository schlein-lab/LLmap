// LLmap — LosslessAggregator tests.
//
// The aggregator's job is to prove "n_records == n_input" and that every
// emitted record passes is_lossless_consistent(). These tests exercise:
//   - basic counter increment
//   - per-status / per-rejection / per-kind breakdown
//   - invariant violation detection
//   - JSON summary round-trip via filesystem
//   - thread safety (concurrent Observe() under a barrier)

#include "output/lossless_aggregator.h"
#include "core/alignment_record.h"

#include <gtest/gtest.h>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>

using namespace llmap;
using llmap::output::LosslessAggregator;

namespace {

AlignmentHit hit(std::string target = "chr14") {
    AlignmentHit h;
    h.target_id = std::move(target);
    h.start = 1'000;
    h.end = 11'000;
    h.cigar.ops = "10000M";
    h.score = 9500;
    h.nm = 50;
    return h;
}

TentativeTarget target() {
    TentativeTarget t;
    t.target_id = "chr14";
    t.approx_start = 1'000;
    t.approx_end = 11'000;
    t.n_seeds = 12;
    t.final_probability = 0.5f;
    return t;
}

}  // namespace

// Reading back the JSON summary as text — minimal helper, we don't need a
// JSON parser dep just to assert presence of substrings.
static std::string ReadFile(const std::filesystem::path& p) {
    std::ifstream in(p);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

TEST(LosslessAggregator, EmptyObservesYieldsZeroSummary) {
    LosslessAggregator agg;
    auto snap = agg.Snapshot();
    EXPECT_EQ(snap.n_records_emitted, 0u);
    EXPECT_TRUE(snap.invariant_ok);
    EXPECT_TRUE(snap.offending_read_ids.empty());
}

TEST(LosslessAggregator, ObserveIncrementsCounters) {
    LosslessAggregator agg;
    agg.SetExpectedInputCount(3);

    EXPECT_TRUE(agg.Observe(make_mapped("r1", 1000, hit())));
    EXPECT_TRUE(agg.Observe(
        make_tentative("r2", 1000, {target()}, RejectionReason::DidNotConverge)));
    EXPECT_TRUE(agg.Observe(make_unmapped("r3", 1000, RejectionReason::NoSeeds)));

    auto snap = agg.Snapshot();
    EXPECT_EQ(snap.n_records_emitted, 3u);
    EXPECT_EQ(snap.n_input_reads_declared, 3u);
    EXPECT_TRUE(snap.invariant_ok);
    EXPECT_EQ(snap.by_status[static_cast<size_t>(AlignmentStatus::Mapped)],    1u);
    EXPECT_EQ(snap.by_status[static_cast<size_t>(AlignmentStatus::Tentative)], 1u);
    EXPECT_EQ(snap.by_status[static_cast<size_t>(AlignmentStatus::Unmapped)],  1u);
    EXPECT_EQ(snap.by_rejection[static_cast<size_t>(RejectionReason::DidNotConverge)], 1u);
    EXPECT_EQ(snap.by_rejection[static_cast<size_t>(RejectionReason::NoSeeds)],         1u);
}

TEST(LosslessAggregator, InvariantViolationFlagsRecord) {
    LosslessAggregator agg;

    // Manually construct an inconsistent record (skipping factory).
    AlignmentRecord bad;
    bad.read_id = "broken";
    bad.status = AlignmentStatus::Mapped;
    // primary intentionally unset

    EXPECT_FALSE(agg.Observe(bad));

    auto snap = agg.Snapshot();
    EXPECT_FALSE(snap.invariant_ok);
    ASSERT_EQ(snap.offending_read_ids.size(), 1u);
    EXPECT_EQ(snap.offending_read_ids[0], "broken");
}

TEST(LosslessAggregator, OffenderListCappedAtMaxOffenders) {
    LosslessAggregator agg;
    constexpr std::size_t too_many = LosslessAggregator::kMaxOffenders + 20;

    for (std::size_t i = 0; i < too_many; ++i) {
        AlignmentRecord bad;
        bad.read_id = "bad_" + std::to_string(i);
        bad.status = AlignmentStatus::Mapped;  // missing primary → invalid
        agg.Observe(bad);
    }

    auto snap = agg.Snapshot();
    EXPECT_FALSE(snap.invariant_ok);
    EXPECT_EQ(snap.offending_read_ids.size(), LosslessAggregator::kMaxOffenders)
        << "Offender list must stop growing past kMaxOffenders.";
    EXPECT_EQ(snap.n_records_emitted, too_many)
        << "Counter must keep counting past the cap.";
}

TEST(LosslessAggregator, JsonSummaryRoundTrip) {
    LosslessAggregator agg;
    agg.SetExpectedInputCount(2);
    agg.Observe(make_mapped("r1", 1000, hit()));
    agg.Observe(make_mapped("r2", 1000, hit()));

    auto tmp = std::filesystem::temp_directory_path() / "llmap_lossless_test.json";
    ASSERT_TRUE(agg.WriteSummary(tmp));
    std::string body = ReadFile(tmp);

    EXPECT_NE(body.find("\"n_input_reads_declared\": 2"), std::string::npos);
    EXPECT_NE(body.find("\"n_records_emitted\":     2"), std::string::npos);
    EXPECT_NE(body.find("\"counts_match\":           true"), std::string::npos);
    EXPECT_NE(body.find("\"lossless_invariant_ok\": true"), std::string::npos);
    EXPECT_NE(body.find("\"MAPPED\": 2"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(LosslessAggregator, JsonSummaryListsOffendersOnBreakage) {
    LosslessAggregator agg;
    AlignmentRecord bad;
    bad.read_id = "violator_42";
    bad.status = AlignmentStatus::Mapped;  // no primary set
    agg.Observe(bad);
    agg.Observe(make_mapped("r_ok", 100, hit()));

    auto tmp = std::filesystem::temp_directory_path() / "llmap_lossless_bad.json";
    ASSERT_TRUE(agg.WriteSummary(tmp));
    std::string body = ReadFile(tmp);

    EXPECT_NE(body.find("\"lossless_invariant_ok\": false"), std::string::npos);
    EXPECT_NE(body.find("violator_42"), std::string::npos);
    std::filesystem::remove(tmp);
}

TEST(LosslessAggregator, ThreadSafeConcurrentObserve) {
    LosslessAggregator agg;
    constexpr int kThreads = 8;
    constexpr int kPerThread = 5000;

    std::atomic<bool> go{false};
    std::vector<std::thread> ts;
    for (int t = 0; t < kThreads; ++t) {
        ts.emplace_back([&, t] {
            while (!go.load(std::memory_order_acquire)) { /* spin */ }
            for (int i = 0; i < kPerThread; ++i) {
                agg.Observe(make_mapped(
                    "t" + std::to_string(t) + "_r" + std::to_string(i),
                    1000, hit()));
            }
        });
    }
    go.store(true, std::memory_order_release);
    for (auto& th : ts) th.join();

    auto snap = agg.Snapshot();
    EXPECT_EQ(snap.n_records_emitted,
              static_cast<std::uint64_t>(kThreads) * kPerThread);
    EXPECT_TRUE(snap.invariant_ok);
}

TEST(LosslessAggregator, ResetClearsState) {
    LosslessAggregator agg;
    agg.SetExpectedInputCount(5);
    agg.Observe(make_mapped("r", 100, hit()));

    agg.Reset();
    auto snap = agg.Snapshot();
    EXPECT_EQ(snap.n_records_emitted, 0u);
    EXPECT_EQ(snap.n_input_reads_declared, 0u);
    EXPECT_TRUE(snap.invariant_ok);
}
