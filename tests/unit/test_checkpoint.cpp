// LLmap -- tests for the Layer 3 checkpoint scaffold.
//
// Covers:
//   - Off mode short-circuits (no agent, no cache lookup, fallback_used=true)
//   - Cache hits work end-to-end
//   - Auto mode falls back silently when the agent is null
//   - Required mode throws when the agent is null and cache misses
//   - Dispatcher stats accumulate correctly across multiple checkpoints
//   - Cache key is deterministic and order-independent on candidates

#include "checkpoint/checkpoint_cache.h"
#include "checkpoint/checkpoint_dispatcher.h"
#include "checkpoint/checkpoint_prompts.h"
#include "checkpoint/checkpoint_tools.h"
#include "checkpoint/checkpoint_types.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <string>

namespace cp = llmap::checkpoint;

namespace {

std::filesystem::path MakeTmpCacheDir(const std::string& tag) {
    auto root = std::filesystem::temp_directory_path() /
                ("llmap_cp_" + tag + "_" +
                 std::to_string(::getpid()) + "_" +
                 std::to_string(std::rand()));
    std::filesystem::create_directories(root);
    return root;
}

cp::CheckpointContext MakeContext(const std::string& read_id = "r1") {
    cp::CheckpointContext ctx;
    ctx.read_id = read_id;
    ctx.ref_id = 7;
    ctx.region_name = "ighg_paralog_cluster";
    ctx.candidate_positions = {{1000u, 42}, {2000u, 41}, {3000u, 40}};
    return ctx;
}

}  // namespace

TEST(CheckpointDispatcher, OffModeShortCircuits) {
    auto dir = MakeTmpCacheDir("off");
    cp::CheckpointCache cache(dir);
    cp::CheckpointDispatcher d(cp::LlmMode::Off, &cache, nullptr);

    auto ctx = MakeContext();
    auto decision = d.Consult(cp::CheckpointType::AmbiguousChain, ctx);

    EXPECT_FALSE(decision.consulted);
    EXPECT_TRUE(decision.fallback_used);
    EXPECT_EQ(d.GetStats().fallback_used, 1u);
    EXPECT_EQ(d.GetStats().cache_hits, 0u);
    EXPECT_EQ(d.GetStats().cache_misses, 0u);
    EXPECT_EQ(d.GetStats().by_type[size_t(cp::CheckpointType::AmbiguousChain)], 1u);

    std::filesystem::remove_all(dir);
}

TEST(CheckpointDispatcher, AutoModeFallsBackOnAgentError) {
    auto dir = MakeTmpCacheDir("auto");
    cp::CheckpointCache cache(dir);
    cp::CheckpointDispatcher d(cp::LlmMode::Auto, &cache, /*agent=*/nullptr);

    auto ctx = MakeContext();
    auto decision = d.Consult(cp::CheckpointType::UnknownRegion, ctx);

    EXPECT_FALSE(decision.consulted);
    EXPECT_TRUE(decision.fallback_used);
    EXPECT_EQ(d.GetStats().cache_misses, 1u);
    EXPECT_EQ(d.GetStats().consulted, 0u);
    EXPECT_EQ(d.GetStats().fallback_used, 1u);

    std::filesystem::remove_all(dir);
}

TEST(CheckpointDispatcher, RequiredModeThrowsWhenAgentMissing) {
    auto dir = MakeTmpCacheDir("required");
    cp::CheckpointCache cache(dir);
    cp::CheckpointDispatcher d(cp::LlmMode::Required, &cache, /*agent=*/nullptr);

    EXPECT_THROW(d.Consult(cp::CheckpointType::ParalogDisambiguation,
                           MakeContext()),
                 cp::CheckpointAgentUnavailable);

    std::filesystem::remove_all(dir);
}

TEST(CheckpointCache, RoundTripsAgentDecision) {
    auto dir = MakeTmpCacheDir("rt");
    cp::CheckpointCache cache(dir);

    cp::AgentDecision in;
    in.consulted = true;
    in.fallback_used = false;
    in.reasoning = "looks like centromeric SD";
    in.special_finding = "novel-sd";
    in.wave = {{1000u, 0.6f}, {2000u, 0.4f}};
    in.override.lambda_scale = 1.5f;
    in.override.report_multi_position = true;

    auto ctx = MakeContext("rt-1");
    auto key = cp::MakeCacheKey(cp::CheckpointType::SDExpansion, ctx);
    ASSERT_TRUE(cache.Store(key, in));

    auto out = cache.Lookup(key);
    ASSERT_TRUE(out.has_value());
    EXPECT_TRUE(out->consulted);
    EXPECT_FALSE(out->fallback_used);
    EXPECT_EQ(out->reasoning, in.reasoning);
    EXPECT_EQ(out->special_finding, in.special_finding);
    ASSERT_EQ(out->wave.size(), 2u);
    EXPECT_EQ(out->wave[0].first, 1000u);
    EXPECT_NEAR(out->wave[0].second, 0.6f, 1e-5);
    EXPECT_EQ(out->wave[1].first, 2000u);
    EXPECT_NEAR(out->wave[1].second, 0.4f, 1e-5);

    std::filesystem::remove_all(dir);
}

TEST(CheckpointDispatcher, CacheHitsAreReturnedWithoutAgent) {
    auto dir = MakeTmpCacheDir("hit");
    cp::CheckpointCache cache(dir);

    auto ctx = MakeContext("cached-read");
    auto key = cp::MakeCacheKey(cp::CheckpointType::AmbiguousChain, ctx);

    cp::AgentDecision stored;
    stored.consulted = true;
    stored.fallback_used = false;
    stored.reasoning = "pre-recorded";
    stored.wave = {{1000u, 1.0f}};
    ASSERT_TRUE(cache.Store(key, stored));

    cp::CheckpointDispatcher d(cp::LlmMode::Auto, &cache, /*agent=*/nullptr);
    auto out = d.Consult(cp::CheckpointType::AmbiguousChain, ctx);

    EXPECT_TRUE(out.consulted);
    EXPECT_FALSE(out.fallback_used);
    EXPECT_EQ(out.reasoning, "pre-recorded");
    EXPECT_EQ(d.GetStats().cache_hits, 1u);
    EXPECT_EQ(d.GetStats().cache_misses, 0u);
    EXPECT_EQ(d.GetStats().consulted, 1u);
    EXPECT_EQ(d.GetStats().fallback_used, 0u);

    std::filesystem::remove_all(dir);
}

TEST(CheckpointDispatcher, StatsAccumulateAcrossTypes) {
    auto dir = MakeTmpCacheDir("stats");
    cp::CheckpointCache cache(dir);
    cp::CheckpointDispatcher d(cp::LlmMode::Auto, &cache, nullptr);

    d.Consult(cp::CheckpointType::AmbiguousChain, MakeContext("a"));
    d.Consult(cp::CheckpointType::AmbiguousChain, MakeContext("b"));
    d.Consult(cp::CheckpointType::NovelInsertion, MakeContext("c"));

    const auto& s = d.GetStats();
    EXPECT_EQ(s.by_type[size_t(cp::CheckpointType::AmbiguousChain)], 2u);
    EXPECT_EQ(s.by_type[size_t(cp::CheckpointType::NovelInsertion)], 1u);
    EXPECT_EQ(s.cache_misses, 3u);
    EXPECT_EQ(s.fallback_used, 3u);

    std::filesystem::remove_all(dir);
}

TEST(CheckpointCache, KeyIsDeterministicAndOrderIndependent) {
    cp::CheckpointContext a = MakeContext("same");
    cp::CheckpointContext b = MakeContext("same");
    // Reverse the candidate order: cache key must be identical because
    // MakeCacheKey sorts candidates internally.
    std::reverse(b.candidate_positions.begin(), b.candidate_positions.end());

    auto k1 = cp::MakeCacheKey(cp::CheckpointType::AmbiguousChain, a);
    auto k2 = cp::MakeCacheKey(cp::CheckpointType::AmbiguousChain, b);
    EXPECT_EQ(k1, k2);

    // Different checkpoint type must yield a different key.
    auto k3 = cp::MakeCacheKey(cp::CheckpointType::NovelInsertion, a);
    EXPECT_NE(k1, k3);
}

TEST(CheckpointTools, RegistryHasFiveStubs) {
    cp::CheckpointToolRegistry reg;
    auto tools = reg.List();
    EXPECT_EQ(tools.size(), 5u);

    auto res = reg.Invoke("region_lookup", {{"region", "chr14_ighg"}});
    EXPECT_EQ(res.status, "not_implemented");
    EXPECT_NE(res.output.find("chr14_ighg"), std::string::npos);

    auto missing = reg.Invoke("region_lookup", {});
    EXPECT_EQ(missing.status, "error");

    auto unknown = reg.Invoke("nope", {});
    EXPECT_EQ(unknown.status, "error");

    EXPECT_FALSE(reg.HasNonStubTools());
}

TEST(CheckpointPrompts, BuildsForEveryType) {
    auto ctx = MakeContext();
    for (int i = 0; i < cp::kCheckpointTypeCount; ++i) {
        auto type = static_cast<cp::CheckpointType>(i);
        auto p = cp::BuildPrompt(type, ctx);
        EXPECT_FALSE(p.empty()) << "empty prompt for type " << i;
        // Each prompt embeds the checkpoint kind tag.
        EXPECT_NE(p.find("CHECKPOINT:"), std::string::npos);
    }
    EXPECT_FALSE(cp::SystemPrompt().empty());
}
