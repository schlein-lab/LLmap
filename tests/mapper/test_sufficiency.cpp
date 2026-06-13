// LLmap — sufficiency (early-exit cascade criterion) tests.

#include "mapper/sufficiency.h"

#include <gtest/gtest.h>

#include <vector>

namespace {

using llmap::mapper::ChainSummary;
using llmap::mapper::IsProvablyResolved;
using llmap::mapper::SufficiencyConfig;

ChainSummary Chain(std::int32_t score, std::uint32_t qs, std::uint32_t qe,
                   bool intron = false) {
    return ChainSummary{score, qs, qe, intron};
}

TEST(Sufficiency, SingleDominantFullCoverageResolved) {
    std::vector<ChainSummary> c = {Chain(1000, 0, 1000)};  // covers whole read
    const auto r = IsProvablyResolved(c, 1000);
    EXPECT_TRUE(r.resolved);
    EXPECT_EQ(r.best_idx, 0u);
    EXPECT_STREQ(r.reason, "single-dominant");
}

TEST(Sufficiency, DominantOverRunnerUpResolved) {
    // best 1000 vs second 300: 1000 >= 1.5*300 → dominant; full coverage.
    std::vector<ChainSummary> c = {Chain(300, 0, 980), Chain(1000, 0, 1000)};
    const auto r = IsProvablyResolved(c, 1000);
    EXPECT_TRUE(r.resolved);
    EXPECT_EQ(r.best_idx, 1u);
}

TEST(Sufficiency, MoleculeZeroFragmentEscalates) {
    // Real-data case: one chain covers 91 bp of a 4419 bp read (molecule/0,
    // CIGAR 1254S 91M 3074S) → partial coverage → escalate, NOT early-exit.
    std::vector<ChainSummary> c = {Chain(91, 1254, 1345)};
    const auto r = IsProvablyResolved(c, 4419);
    EXPECT_FALSE(r.resolved);
    EXPECT_STREQ(r.reason, "partial-coverage");
}

TEST(Sufficiency, CompetingLocusEscalates) {
    // Two comparable full-coverage chains → genuinely ambiguous → escalate.
    std::vector<ChainSummary> c = {Chain(1000, 0, 1000), Chain(950, 0, 1000)};
    const auto r = IsProvablyResolved(c, 1000);
    EXPECT_FALSE(r.resolved);
    EXPECT_STREQ(r.reason, "competing-locus");
}

TEST(Sufficiency, IntronSignatureEscalatesToSplice) {
    // A dominant, full-coverage chain BUT with an intron gap → spliced read,
    // must escalate to the splice path (never early-exit as "easy").
    std::vector<ChainSummary> c = {Chain(1000, 0, 1000, /*intron=*/true)};
    const auto r = IsProvablyResolved(c, 1000);
    EXPECT_FALSE(r.resolved);
    EXPECT_STREQ(r.reason, "intron-signature");
}

TEST(Sufficiency, EmptyAndZeroLength) {
    EXPECT_FALSE(IsProvablyResolved({}, 1000).resolved);
    std::vector<ChainSummary> c = {Chain(100, 0, 100)};
    EXPECT_FALSE(IsProvablyResolved(c, 0).resolved);  // no query length
}

}  // namespace
