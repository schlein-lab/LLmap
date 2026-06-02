// LLmap — AidFootprintDetector tests.
//
// Synthesise a Sγ-region anchor + a read that carries AID-mediated
// C→T edits in WRC hotspot motifs. Verify the detector picks them up
// and that confidence reflects hotspot fraction.

#include "rnamod/aid_footprint.h"

#include <gtest/gtest.h>

#include <string>

using namespace llmap::rnamod;

namespace {

// Build a synthetic Sγ4-like sequence: GC-rich, tandem WRC motifs.
//   anchor:  AGCT  AGCT  AGCT  AGCT  AGCT  AGCT  AGCT  AGCT (32 nt)
//   read   : same, but C→T at every WRC-hotspot C position
//
// In "AGCT" the C at offset 2 is preceded by AG (R=A is purine, but we
// want the leading dinucleotide WR — A and G — for AID hotspot).
// Actually AID hotspot is WRC: W=A/T, R=A/G, C is the cytosine being
// deaminated. So context is preceding 2bp = WR + the C itself. For
// AGCT we have A(W?yes) G(R?yes) C → AID hotspot.
constexpr const char* kAnchor =
    "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";  // 32 nt
constexpr const char* kReadAllEdits =
    "AGTTAGTTAGTTAGTTAGTTAGTTAGTTAGTT";  // all C→T at hotspot positions
constexpr const char* kReadNoEdits =
    "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";  // identical to anchor

}  // namespace

TEST(AidFootprint, DetectsAllHotspotEdits) {
    AidFootprintDetector det;
    auto r = det.Detect(kReadAllEdits, kAnchor, "S_gamma4");
    EXPECT_TRUE(r.detected);
    EXPECT_GE(r.n_c_to_u_events, 3u);
    EXPECT_GT(r.confidence, 0.9f) << "all edits in hotspot → confidence ~1.0";
    EXPECT_EQ(r.switch_region_id, "S_gamma4");
}

TEST(AidFootprint, NoEditsReturnsNotDetected) {
    AidFootprintDetector det;
    auto r = det.Detect(kReadNoEdits, kAnchor, "S_gamma4");
    EXPECT_FALSE(r.detected);
    EXPECT_EQ(r.n_c_to_u_events, 0u);
}

TEST(AidFootprint, FewerThanThreeEditsNotDetected) {
    AidFootprintDetector det;
    // Only 2 edits — below threshold even though both are in hotspot
    std::string r2 =
        "AGTTAGTTAGCTAGCTAGCTAGCTAGCTAGCT";  // 2 C→T at positions 2,6
    auto r = det.Detect(r2, kAnchor, "S_gamma4");
    EXPECT_FALSE(r.detected);
    EXPECT_EQ(r.n_c_to_u_events, 2u)
        << "events still counted, but below the 3-event threshold";
}

TEST(AidFootprint, EmptySwitchRegionIdReturnsImmediately) {
    AidFootprintDetector det;
    auto r = det.Detect(kReadAllEdits, kAnchor, /*switch_id=*/"");
    EXPECT_FALSE(r.detected);
    EXPECT_EQ(r.n_c_to_u_events, 0u);
}

TEST(AidFootprint, LengthMismatchReturnsNotDetected) {
    AidFootprintDetector det;
    std::string short_read = "AGTT";
    auto r = det.Detect(short_read, kAnchor, "S_gamma4");
    EXPECT_FALSE(r.detected);
    EXPECT_EQ(r.n_c_to_u_events, 0u);
}
