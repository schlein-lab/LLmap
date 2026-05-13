// LLmap Phase 0 — BucketPyramid tests.
//
// Tests the 3-level bucket hierarchy, parent indices, biology hints,
// and serialization round-trip equality.

#include "core/bucket_pyramid.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <cstdio>

using namespace llmap;

namespace {

class BucketPyramidTest : public ::testing::Test {
protected:
    std::filesystem::path temp_file_;

    void SetUp() override {
        temp_file_ = std::filesystem::temp_directory_path() / "test_bucket_pyramid.bin";
    }

    void TearDown() override {
        if (std::filesystem::exists(temp_file_)) {
            std::filesystem::remove(temp_file_);
        }
    }

    BucketPyramid build_sample_pyramid() {
        BucketPyramid p;

        // L0: 2 chromosomes
        L0Bucket l0_chr14;
        l0_chr14.name = "chr14";
        l0_chr14.total_span = 107'000'000;
        p.add_l0_bucket(std::move(l0_chr14));

        L0Bucket l0_chr15;
        l0_chr15.name = "chr15";
        l0_chr15.total_span = 101'000'000;
        p.add_l0_bucket(std::move(l0_chr15));

        // L1: 3 regions (2 in chr14, 1 in chr15)
        L1Bucket l1_0;
        l1_0.target_id = "chr14";
        l1_0.start = 0;
        l1_0.end = 5'000'000;
        p.add_l1_bucket(std::move(l1_0), 0);  // parent = chr14

        L1Bucket l1_1;
        l1_1.target_id = "chr14";
        l1_1.start = 105'000'000;
        l1_1.end = 107'000'000;
        p.add_l1_bucket(std::move(l1_1), 0);  // parent = chr14, IGH region

        L1Bucket l1_2;
        l1_2.target_id = "chr15";
        l1_2.start = 0;
        l1_2.end = 5'000'000;
        p.add_l1_bucket(std::move(l1_2), 1);  // parent = chr15

        // L2: 4 fine-grained buckets
        L2Bucket l2_0;
        l2_0.target_id = "chr14";
        l2_0.start = 0;
        l2_0.end = 50'000;
        p.add_l2_bucket(std::move(l2_0), 0);  // parent = l1_0

        L2Bucket l2_1;
        l2_1.target_id = "chr14";
        l2_1.start = 105'000'000;
        l2_1.end = 105'050'000;
        p.add_l2_bucket(std::move(l2_1), 1);  // parent = l1_1 (IGH)

        L2Bucket l2_2;
        l2_2.target_id = "chr14";
        l2_2.start = 105'050'000;
        l2_2.end = 105'100'000;
        p.add_l2_bucket(std::move(l2_2), 1);  // parent = l1_1 (IGH)

        L2Bucket l2_3;
        l2_3.target_id = "chr15";
        l2_3.start = 0;
        l2_3.end = 50'000;
        p.add_l2_bucket(std::move(l2_3), 2);  // parent = l1_2

        // Biology hints for IGH buckets
        BiologyHint hint1;
        hint1.prior_weight = 1.5f;
        hint1.annotation = "IGHG-Constants";
        hint1.paralog_partner_bucket = 2;
        hint1.expected_coverage_multiplier = 2.0f;
        p.set_biology_hint(1, std::move(hint1));

        BiologyHint hint2;
        hint2.prior_weight = 1.5f;
        hint2.annotation = "IGHG-Constants-Dup";
        hint2.paralog_partner_bucket = 1;
        hint2.expected_coverage_multiplier = 2.0f;
        p.set_biology_hint(2, std::move(hint2));

        return p;
    }
};

}  // namespace

TEST_F(BucketPyramidTest, EmptyPyramidIsValid) {
    BucketPyramid p;
    EXPECT_TRUE(p.validate());
    EXPECT_EQ(p.l0_count(), 0u);
    EXPECT_EQ(p.l1_count(), 0u);
    EXPECT_EQ(p.l2_count(), 0u);
}

TEST_F(BucketPyramidTest, BuildsCorrectHierarchy) {
    auto p = build_sample_pyramid();
    EXPECT_EQ(p.l0_count(), 2u);
    EXPECT_EQ(p.l1_count(), 3u);
    EXPECT_EQ(p.l2_count(), 4u);
    EXPECT_TRUE(p.validate());
}

TEST_F(BucketPyramidTest, ParentIndicesCorrect) {
    auto p = build_sample_pyramid();

    // L1 parents
    EXPECT_EQ(p.l1_parent(0), 0u);  // l1_0 → chr14
    EXPECT_EQ(p.l1_parent(1), 0u);  // l1_1 → chr14
    EXPECT_EQ(p.l1_parent(2), 1u);  // l1_2 → chr15

    // L2 parents
    EXPECT_EQ(p.l2_parent(0), 0u);  // l2_0 → l1_0
    EXPECT_EQ(p.l2_parent(1), 1u);  // l2_1 → l1_1
    EXPECT_EQ(p.l2_parent(2), 1u);  // l2_2 → l1_1
    EXPECT_EQ(p.l2_parent(3), 2u);  // l2_3 → l1_2
}

TEST_F(BucketPyramidTest, BiologyHintsRetrieval) {
    auto p = build_sample_pyramid();

    auto hint1 = p.get_biology_hint(1);
    ASSERT_TRUE(hint1.has_value());
    EXPECT_FLOAT_EQ(hint1->prior_weight, 1.5f);
    EXPECT_EQ(hint1->annotation, "IGHG-Constants");
    ASSERT_TRUE(hint1->paralog_partner_bucket.has_value());
    EXPECT_EQ(*hint1->paralog_partner_bucket, 2u);
    EXPECT_FLOAT_EQ(hint1->expected_coverage_multiplier, 2.0f);

    auto hint_none = p.get_biology_hint(0);
    EXPECT_FALSE(hint_none.has_value());

    auto hint_invalid = p.get_biology_hint(999);
    EXPECT_FALSE(hint_invalid.has_value());
}

TEST_F(BucketPyramidTest, InvalidParentIndexDetected) {
    BucketPyramid p;

    L0Bucket l0;
    l0.name = "chr14";
    p.add_l0_bucket(std::move(l0));

    L1Bucket l1;
    l1.target_id = "chr14";
    p.add_l1_bucket(std::move(l1), 99);  // invalid parent

    EXPECT_FALSE(p.validate());
}

TEST_F(BucketPyramidTest, SerializeDeserializeRoundTrip) {
    auto original = build_sample_pyramid();
    ASSERT_TRUE(original.validate());

    original.serialize(temp_file_);
    ASSERT_TRUE(std::filesystem::exists(temp_file_));

    auto loaded = BucketPyramid::deserialize(temp_file_);
    EXPECT_TRUE(loaded.validate());
    EXPECT_TRUE(original == loaded);
}

TEST_F(BucketPyramidTest, EmptyPyramidRoundTrip) {
    BucketPyramid empty;
    empty.serialize(temp_file_);

    auto loaded = BucketPyramid::deserialize(temp_file_);
    EXPECT_TRUE(loaded.validate());
    EXPECT_EQ(loaded.l0_count(), 0u);
    EXPECT_EQ(loaded.l1_count(), 0u);
    EXPECT_EQ(loaded.l2_count(), 0u);
    EXPECT_TRUE(empty == loaded);
}

TEST_F(BucketPyramidTest, SpanAccessorsReturnCorrectData) {
    auto p = build_sample_pyramid();

    auto l0_span = p.l0_buckets();
    ASSERT_EQ(l0_span.size(), 2u);
    EXPECT_EQ(l0_span[0].name, "chr14");
    EXPECT_EQ(l0_span[1].name, "chr15");

    auto l1_span = p.l1_buckets();
    ASSERT_EQ(l1_span.size(), 3u);
    EXPECT_EQ(l1_span[0].target_id, "chr14");
    EXPECT_EQ(l1_span[0].end, 5'000'000u);

    auto l2_span = p.l2_buckets();
    ASSERT_EQ(l2_span.size(), 4u);
    EXPECT_EQ(l2_span[1].start, 105'000'000u);
}

TEST_F(BucketPyramidTest, ParentOutOfRangeThrows) {
    auto p = build_sample_pyramid();
    EXPECT_THROW(p.l1_parent(99), std::out_of_range);
    EXPECT_THROW(p.l2_parent(99), std::out_of_range);
}

TEST_F(BucketPyramidTest, DeserializeInvalidMagicThrows) {
    // Write garbage to file
    {
        std::ofstream ofs(temp_file_, std::ios::binary);
        std::uint32_t bad_magic = 0xDEADBEEF;
        ofs.write(reinterpret_cast<const char*>(&bad_magic), sizeof(bad_magic));
    }
    EXPECT_THROW(BucketPyramid::deserialize(temp_file_), std::runtime_error);
}

TEST_F(BucketPyramidTest, HierarchyChainWorks) {
    auto p = build_sample_pyramid();

    // Trace l2_bucket[2] back to l0
    std::uint32_t l2_id = 2;  // IGHG-Constants-Dup bucket
    std::uint32_t l1_id = p.l2_parent(l2_id);
    std::uint32_t l0_id = p.l1_parent(l1_id);

    auto l0_span = p.l0_buckets();
    EXPECT_EQ(l0_span[l0_id].name, "chr14");

    auto l1_span = p.l1_buckets();
    EXPECT_EQ(l1_span[l1_id].start, 105'000'000u);  // IGH region
}
