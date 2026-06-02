// LLmap — Unit tests for the SegDup fallback hierarchy.
//
// Each test wires a minimal StageProbes table with deterministic
// mocks, builds a FallbackChain from StageParams, and verifies the
// chain resolves at the expected rung.
//
// Test matrix:
//   1. Stage 1 success — 3-mismatch read at relaxed mm=4.
//   2. Stage 2 success — extension diverges, chain-only takes over.
//   3. Stage 3 success — two tied chains; multi-position reports both.
//   4. Stage 4 skipped — --llm-fallback not set; chain skips to stage 5.
//   5. Stage 5 flag    — no stage matched; novel_haplotype emitted.

#include "mapper/fallback_chain.h"
#include "mapper/fallback_stages.h"
#include "mapper/fallback_types.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

using namespace llmap::mapper::fallback;

namespace {

/// Build the canonical 5-stage chain matching the IGHG4 catalog.
std::vector<StageParams> CanonicalParams() {
    std::vector<StageParams> p(5);
    p[0].id = StageId::relaxed_mismatch;
    p[0].kmer_size = 25;
    p[0].max_mismatch = 4;

    p[1].id = StageId::chain_only;
    p[1].use_extension = false;

    p[2].id = StageId::multi_position;
    p[2].report_all_top_k = true;
    p[2].top_k = 4;

    p[3].id = StageId::llm_checkpoint;
    p[3].require_llm_flag = true;
    p[3].opt_in_flag = "--llm-fallback";

    p[4].id = StageId::novel_haplotype_flag;
    p[4].emit_warning = true;
    p[4].warning_tag = "novel_haplotype";
    return p;
}

/// Convenience: a fake read sequence — content is opaque to the stages
/// because we inject all probe behaviour via mocks.
std::vector<uint8_t> FakeRead(std::size_t n = 100) {
    return std::vector<uint8_t>(n, 0x42);
}

CandidatePlacement MakePlacement(int32_t score, const char* tag = "canonical") {
    CandidatePlacement p;
    p.ref_id = 14;
    p.ref_start = 105625772;
    p.ref_end = 105625872;
    p.score = score;
    p.mapq = 60;
    p.copy_tag = tag;
    return p;
}

/// Test fixture that resets the global probe table to all-no-op
/// before each test, so cross-test state never leaks.
class FallbackChainTest : public ::testing::Test {
protected:
    void SetUp() override {
        SetGlobalStageProbes(DefaultStageProbes());  // all no-ops
    }
};

}  // namespace

// ---------------------------------------------------------------------------
// Test 1: Stage 1 (relaxed_mismatch) succeeds — 3 mismatches in 25-mer.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage1RelaxedMismatchSucceeds) {
    StageProbes probes = DefaultStageProbes();
    probes.relaxed_mismatch_probe =
        [](std::span<const uint8_t>, uint32_t k, uint32_t mm) {
            // Simulate: primary (mm=2) would have failed; relaxed mm=4 hits.
            EXPECT_EQ(k, 25u);
            EXPECT_GE(mm, 3u);
            return std::vector<CandidatePlacement>{MakePlacement(95)};
        };
    SetGlobalStageProbes(std::move(probes));

    auto params = CanonicalParams();
    auto chain = FallbackChain::FromParams(params);
    ASSERT_EQ(chain.size(), 5u);

    auto seq = FakeRead();
    MappingContext ctx;
    ctx.read_id = "test_stage1";
    ctx.locus_id = "IGHG4_chimdup_tandem";
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_TRUE(outcome.success);
    EXPECT_EQ(outcome.resolved_by, StageId::relaxed_mismatch);
    EXPECT_EQ(outcome.stages_run, 1u);
    ASSERT_EQ(ctx.placements.size(), 1u);
    EXPECT_EQ(ctx.placements[0].score, 95);
    EXPECT_TRUE(ctx.placements[0].extension_applied);
}

// ---------------------------------------------------------------------------
// Test 2: Stage 2 (chain_only) succeeds after stage 1 fails.
//         Extension diverges -> chain-only finds the chain.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage2ChainOnlySucceedsAfterStage1Fail) {
    StageProbes probes = DefaultStageProbes();
    // Stage 1: relaxed mismatch finds nothing (extension would have
    // exploded in real pipeline).
    probes.relaxed_mismatch_probe =
        [](std::span<const uint8_t>, uint32_t, uint32_t) {
            return std::vector<CandidatePlacement>{};
        };
    // Stage 2: chain-only finds the chain.
    probes.chain_only_probe = [](std::span<const uint8_t>) {
        return std::vector<CandidatePlacement>{MakePlacement(72, "dup")};
    };
    SetGlobalStageProbes(std::move(probes));

    auto chain = FallbackChain::FromParams(CanonicalParams());
    auto seq = FakeRead();
    MappingContext ctx;
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_TRUE(outcome.success);
    EXPECT_EQ(outcome.resolved_by, StageId::chain_only);
    EXPECT_EQ(outcome.stages_run, 2u);
    ASSERT_EQ(ctx.placements.size(), 1u);
    // Chain-only stage must clear the extension_applied flag.
    EXPECT_FALSE(ctx.placements[0].extension_applied);
    EXPECT_EQ(ctx.placements[0].copy_tag, "dup");
}

// ---------------------------------------------------------------------------
// Test 3: Stage 3 (multi_position) — two tied chains; both reported.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage3MultiPositionReportsTiedChains) {
    StageProbes probes = DefaultStageProbes();
    probes.relaxed_mismatch_probe =
        [](std::span<const uint8_t>, uint32_t, uint32_t) {
            return std::vector<CandidatePlacement>{};
        };
    probes.chain_only_probe = [](std::span<const uint8_t>) {
        return std::vector<CandidatePlacement>{};
    };
    probes.multi_position_probe =
        [](std::span<const uint8_t>, uint32_t k, float delta) {
            EXPECT_EQ(k, 4u);
            EXPECT_GT(delta, 0.0f);
            return std::vector<CandidatePlacement>{
                MakePlacement(88, "canonical"),
                MakePlacement(88, "dup"),
            };
        };
    SetGlobalStageProbes(std::move(probes));

    auto chain = FallbackChain::FromParams(CanonicalParams());
    auto seq = FakeRead();
    MappingContext ctx;
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_TRUE(outcome.success);
    EXPECT_EQ(outcome.resolved_by, StageId::multi_position);
    EXPECT_EQ(outcome.stages_run, 3u);
    ASSERT_EQ(ctx.placements.size(), 2u);
    EXPECT_EQ(ctx.placements[0].copy_tag, "canonical");
    EXPECT_EQ(ctx.placements[1].copy_tag, "dup");
}

// ---------------------------------------------------------------------------
// Test 3b: Stage 3 disabled by catalog (report_all_top_k = false)
//          -> falls through to stage 5.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage3SkippedWhenCatalogDisablesMultiPosition) {
    auto params = CanonicalParams();
    params[2].report_all_top_k = false;

    StageProbes probes = DefaultStageProbes();  // all no-op
    SetGlobalStageProbes(std::move(probes));

    auto chain = FallbackChain::FromParams(params);
    auto seq = FakeRead();
    MappingContext ctx;
    ctx.llm_fallback_enabled = false;
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_FALSE(outcome.success);
    EXPECT_EQ(outcome.stages_run, 5u);
    ASSERT_EQ(outcome.flags.size(), 1u);
    EXPECT_EQ(outcome.flags[0], "novel_haplotype");
}

// ---------------------------------------------------------------------------
// Test 4: Stage 4 (llm_checkpoint) skipped when --llm-fallback not set.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage4SkippedWithoutLlmFallbackFlag) {
    // Install probes that all fail except an LLM probe that WOULD
    // succeed if invoked. We expect the chain to skip it.
    StageProbes probes = DefaultStageProbes();
    bool llm_probe_called = false;
    probes.llm_checkpoint_probe =
        [&](std::span<const uint8_t>, std::string_view) {
            llm_probe_called = true;
            return std::vector<CandidatePlacement>{MakePlacement(50)};
        };
    SetGlobalStageProbes(std::move(probes));

    auto chain = FallbackChain::FromParams(CanonicalParams());
    auto seq = FakeRead();
    MappingContext ctx;
    ctx.llm_fallback_enabled = false;  // explicit
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_FALSE(outcome.success);
    EXPECT_FALSE(llm_probe_called);
    EXPECT_EQ(outcome.stages_run, 5u);  // all five attempted
    ASSERT_EQ(outcome.flags.size(), 1u);
    EXPECT_EQ(outcome.flags[0], "novel_haplotype");
}

// ---------------------------------------------------------------------------
// Test 4b: Stage 4 invoked when --llm-fallback is set; stub returns
//          empty so we still fall through to stage 5.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage4InvokedWithLlmFallbackFlagButStubFails) {
    bool llm_probe_called = false;
    StageProbes probes = DefaultStageProbes();
    probes.llm_checkpoint_probe =
        [&](std::span<const uint8_t>, std::string_view) {
            llm_probe_called = true;
            return std::vector<CandidatePlacement>{};  // stub: no result
        };
    SetGlobalStageProbes(std::move(probes));

    auto chain = FallbackChain::FromParams(CanonicalParams());
    auto seq = FakeRead();
    MappingContext ctx;
    ctx.llm_fallback_enabled = true;
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_FALSE(outcome.success);
    EXPECT_TRUE(llm_probe_called);  // was invoked this time
    EXPECT_EQ(outcome.stages_run, 5u);
}

// ---------------------------------------------------------------------------
// Test 5: Stage 5 — all prior stages fail; novel_haplotype flag emitted.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, Stage5FlagEmittedWhenAllPriorFail) {
    SetGlobalStageProbes(DefaultStageProbes());  // all no-op

    auto chain = FallbackChain::FromParams(CanonicalParams());
    auto seq = FakeRead();
    MappingContext ctx;
    ctx.read_seq = std::span<const uint8_t>(seq);

    auto outcome = chain.run(ctx);
    EXPECT_FALSE(outcome.success);
    EXPECT_EQ(outcome.stages_run, 5u);
    ASSERT_EQ(outcome.flags.size(), 1u);
    EXPECT_EQ(outcome.flags[0], "novel_haplotype");
    EXPECT_EQ(ctx.last_stage_attempted, StageId::novel_haplotype_flag);
}

// ---------------------------------------------------------------------------
// Test 6: Factory contract — MakeStage returns the correct concrete type
//         for every known id; nullptr semantics for catalog mistakes.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, MakeStageBuildsEveryKnownId) {
    for (auto id : {StageId::relaxed_mismatch,
                    StageId::chain_only,
                    StageId::multi_position,
                    StageId::llm_checkpoint,
                    StageId::novel_haplotype_flag}) {
        StageParams p;
        p.id = id;
        auto stage = MakeStage(p);
        ASSERT_NE(stage, nullptr) << "MakeStage failed for "
                                  << StageIdName(id);
        EXPECT_EQ(stage->id(), id);
    }
}

// ---------------------------------------------------------------------------
// Test 7: Stage names round-trip. Catalog JSON uses the string names, so
//         the loader will rely on this mapping being stable.
// ---------------------------------------------------------------------------

TEST_F(FallbackChainTest, StageNamesAreStable) {
    EXPECT_STREQ(StageIdName(StageId::relaxed_mismatch),     "relaxed_mismatch");
    EXPECT_STREQ(StageIdName(StageId::chain_only),           "chain_only");
    EXPECT_STREQ(StageIdName(StageId::multi_position),       "multi_position");
    EXPECT_STREQ(StageIdName(StageId::llm_checkpoint),       "llm_checkpoint");
    EXPECT_STREQ(StageIdName(StageId::novel_haplotype_flag), "novel_haplotype_flag");
}
