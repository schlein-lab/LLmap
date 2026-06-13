// LLmap — mei_detector tests (Block 2 (i)).
//
// Non-reference MEI (novel Alu/L1/SVA insertion) → Layer-3 `bio:mei` + finding,
// distinguished from reference-TE confusion (Layer-1, deferred). Diagnostic
// signature: split + TE-consensus clip + (poly-A | canonical TSD).

#include "provenance/mei_detector.h"

#include <gtest/gtest.h>

namespace {

using llmap::provenance::ClassifyMei;
using llmap::provenance::MeiEvidence;
using llmap::provenance::TeFamily;
using llmap::provenance::TeFamilyTag;

TEST(MeiDetector, NovelAluInsertionWithPolyAAndTsd) {
    MeiEvidence e;
    e.is_split = true;
    e.clip_matches_te = true;
    e.te_family = TeFamily::Alu;
    e.has_polyA = true;
    e.tsd_length = 12;  // canonical TSD
    const auto c = ClassifyMei(e);
    EXPECT_TRUE(c.is_mei);
    EXPECT_EQ(c.te_family, TeFamily::Alu);
    EXPECT_TRUE(c.has_polyA);
    EXPECT_TRUE(c.has_tsd);
    EXPECT_GE(c.confidence, 0.85f);
}

TEST(MeiDetector, BareTeClipDefersToReferenceConfusion) {
    // Split + TE clip but NO insertion hallmark (no poly-A, no TSD) → not a novel
    // MEI; the Layer-1 reference-TE-confusion path handles it.
    MeiEvidence e;
    e.is_split = true;
    e.clip_matches_te = true;
    e.te_family = TeFamily::L1;
    const auto c = ClassifyMei(e);
    EXPECT_FALSE(c.is_mei);
}

TEST(MeiDetector, NotSplitIsNotMei) {
    MeiEvidence e;
    e.clip_matches_te = true;
    e.has_polyA = true;
    e.tsd_length = 10;
    // is_split = false
    EXPECT_FALSE(ClassifyMei(e).is_mei);
}

TEST(MeiDetector, StrongHallmarksOverrideReferenceTeDiscount) {
    // At an annotated reference TE, but with BOTH poly-A and a canonical TSD →
    // a genuine novel insertion overrides the reference-confusion discount.
    MeiEvidence e;
    e.is_split = true;
    e.clip_matches_te = true;
    e.te_family = TeFamily::L1;
    e.has_polyA = true;
    e.tsd_length = 15;
    e.five_prime_truncated = true;
    e.at_reference_te = true;
    const auto c = ClassifyMei(e);
    EXPECT_TRUE(c.is_mei);  // 0.4+0.3+0.2+0.1-0.15 = 0.85 ≥ 0.60
}

TEST(MeiDetector, ReferenceTeBareClipNotMei) {
    MeiEvidence e;
    e.is_split = true;
    e.clip_matches_te = true;
    e.at_reference_te = true;  // bare clip at a ref TE → confusion, not novel MEI
    EXPECT_FALSE(ClassifyMei(e).is_mei);
}

TEST(MeiDetector, TsdOutOfWindowNotCounted) {
    MeiEvidence e;
    e.is_split = true;
    e.clip_matches_te = true;
    e.tsd_length = 100;  // far above the canonical 4-20 bp window → not a TSD
    const auto c = ClassifyMei(e);
    EXPECT_FALSE(c.has_tsd);
    EXPECT_FALSE(c.is_mei);  // no poly-A, no valid TSD
}

TEST(MeiDetector, FamilyTags) {
    EXPECT_STREQ(TeFamilyTag(TeFamily::Alu), "alu");
    EXPECT_STREQ(TeFamilyTag(TeFamily::L1), "l1");
    EXPECT_STREQ(TeFamilyTag(TeFamily::Sva), "sva");
    EXPECT_STREQ(TeFamilyTag(TeFamily::Herv), "herv");
    EXPECT_STREQ(TeFamilyTag(TeFamily::Unknown), "te");
}

}  // namespace
