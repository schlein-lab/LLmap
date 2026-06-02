// LLmap — SpliceSiteDb tests.
//
// Validates the textbook PWM scoring against synthesised canonical and
// scrambled junctions; verifies U2/U12 discrimination; checks back-splice
// detection; runs a sanity check on the polypyrimidine + branch-point
// scanners.

#include "annot/splice_site_db.h"

#include <gtest/gtest.h>

#include <string>

using namespace llmap::annot;

namespace {

SpliceSiteDb MakeDb() {
    SpliceSiteDb d;
    d.LoadDefaults();
    return d;
}

}  // namespace

TEST(SpliceSiteDb, CanonicalU2DonorScoresHigh) {
    auto db = MakeDb();
    // donor 'GT' + intron 5p 'AAGT...' — canonical-ish
    auto r = db.ScoreJunction(/*donor*/ "GT",
                                /*acceptor*/ "AG",
                                /*intron_3p*/ "CTTCTCCTTCAG",
                                /*intron_5p*/ "AAGTAAGT");
    EXPECT_GT(r.donor_score, 0.5f)
        << "canonical GT|RAGT donor should score >0.5 (got "
        << r.donor_score << ")";
    EXPECT_GT(r.acceptor_score, 0.4f);
    EXPECT_EQ(r.spliceosome_class, 0u);  // U2 major
}

TEST(SpliceSiteDb, ScrambledJunctionScoresLow) {
    auto db = MakeDb();
    auto r = db.ScoreJunction("CC",   // not GT or AT
                                "CC",   // not AG or AC
                                "AAAAAAAAAACC",
                                "CCCCCCCC");
    EXPECT_LT(r.donor_score,    0.3f);
    EXPECT_LT(r.acceptor_score, 0.3f);
    EXPECT_EQ(r.spliceosome_class, 2u);  // non-canonical
}

TEST(SpliceSiteDb, U12MinorClassDetected) {
    auto db = MakeDb();
    // AT|ATCCTT donor, TCCTTRAC|AC acceptor (R=A) → strong U12
    auto r = db.ScoreJunction("AT", "AC",
                                "ATCCTTAC",   // acceptor consensus
                                "ATCCTTAA");  // donor consensus
    EXPECT_EQ(r.spliceosome_class, 1u)
        << "expected U12 minor class; donor=" << r.donor_score
        << " acc=" << r.acceptor_score;
}

TEST(SpliceSiteDb, PolypyrimidineDetection) {
    EXPECT_NEAR(SpliceSiteDb::PolypyrimidineScore("CCCCTTCCTC"), 1.0f, 1e-5f);
    EXPECT_NEAR(SpliceSiteDb::PolypyrimidineScore("AAAGGGAAA"),  0.0f, 1e-5f);
    EXPECT_NEAR(SpliceSiteDb::PolypyrimidineScore(""),           0.0f, 1e-5f);
    EXPECT_NEAR(SpliceSiteDb::PolypyrimidineScore("CCAA"),       0.5f, 1e-5f);
}

TEST(SpliceSiteDb, BranchPointReportedWhenPresent) {
    auto db = MakeDb();
    // YNYURAC at offset ~-30 before acceptor — should be detected
    std::string intron_3p = "GCATATCTTACAG";  // ends with ACAG (acceptor)
    // The PolyPyr/branch-point detector requires ≥7-mer; we pad in front
    std::string padded = "AAGAAAAAGCAT" + intron_3p;
    auto r = db.ScoreJunction("GT", "AG",
                                padded,
                                "GTAAGT");
    // Branch-point offset is set only when the scanner finds ≥0.5 match;
    // synthetic Y-rich region above should trip it.
    // We don't pin the exact offset (PWM-dependent), but it should be ≤0
    // (i.e. valid offset reported).
    EXPECT_LE(r.branch_point_offset, 0);
}

TEST(SpliceSiteDb, IsBackSpliceConsistentRequiresAcceptorBeforeDonor) {
    auto db = MakeDb();
    // canonical: acceptor>donor — NOT a back-splice
    EXPECT_FALSE(db.IsBackSpliceConsistent(
        /*donor*/ 1000, /*acceptor*/ 5000,
        "GT", "AG"));
    // back-splice: acceptor<donor with canonical motifs
    EXPECT_TRUE(db.IsBackSpliceConsistent(
        /*donor*/ 5000, /*acceptor*/ 1000,
        "GT", "AG"));
    // back-splice with bogus motifs
    EXPECT_FALSE(db.IsBackSpliceConsistent(
        5000, 1000, "TT", "TT"));
    // GC/AC variant — still canonical-class
    EXPECT_TRUE(db.IsBackSpliceConsistent(
        5000, 1000, "GC", "AC"));
}

TEST(SpliceSiteDb, MissingDefaultsReturnsZeros) {
    SpliceSiteDb d;  // LoadDefaults NOT called
    auto r = d.ScoreJunction("GT", "AG", "CCCCC", "GTAAGT");
    EXPECT_EQ(r.donor_score, 0.0f);
    EXPECT_EQ(r.acceptor_score, 0.0f);
}
