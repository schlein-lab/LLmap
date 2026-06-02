// LLmap — unit tests for the IGH-mode upstream-anchor classifier.
//
// Tests cover the contract documented in `src/classify/upstream_anchor.h`:
//
//   (a) Read with CH1 + canonical upstream promoter   → canonical
//   (b) Read with CH1 + chimdup upstream promoter     → chimdup
//   (c) Read with CH1 only (no flank)                 → none, mismatches=0
//
// The in-memory mock catalog plays the role of the real SegDup loader
// (still in flight in `src/catalog/`). It implements `ISegDupCatalog`
// directly; when the real loader lands the same tests should keep passing
// if the loader is wired with the same `IGHG4_chimdup_tandem` entry.

#include "classify/segdup_catalog_iface.h"
#include "classify/upstream_anchor.h"

#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

#include <gtest/gtest.h>

using llmap::classify::igh::CatalogEntryView;
using llmap::classify::igh::ClassificationResult;
using llmap::classify::igh::Classify;
using llmap::classify::igh::ClassifyOptions;
using llmap::classify::igh::FlankingKmerMismatches;
using llmap::classify::igh::ISegDupCatalog;
using llmap::classify::igh::MotifMatchesWindow;
using llmap::classify::igh::ParseRelativeWindow;
using llmap::classify::igh::PromoterMotifMatch;
using llmap::classify::igh::PromoterSignature;
using llmap::classify::igh::ReadInput;
using llmap::classify::igh::ToString;

namespace {

// ---------------------------------------------------------------------------
// Test fixtures.
// ---------------------------------------------------------------------------

// Verbatim motifs from catalog/curated/IGHG4_chimdup_tandem.json. Hard-coded
// *only here* in the test so the implementation has zero motif knowledge.
constexpr std::string_view kCanonicalMotif =
    "CGGTTCTT...GTCTATCTGCGATGG";
constexpr std::string_view kDuplicateMotif =
    "GGGTTCTT...AACTGTCCGCGAGG";

// Synthetic CH1 (uppercase, A/C/G/T only). 50 bp; identity of bases is
// irrelevant for these tests.
constexpr std::string_view kSyntheticCh1 =
    "ACACACGTGTACACGTACGTACACGTGTACACGTACGTACACGTGTACAC";

// Build a synthetic 100 bp upstream window that contains the canonical
// motif's prefix and suffix in order, with random filler in between.
std::string CanonicalUpstream() {
    // prefix = CGGTTCTT (8 bp), suffix = GTCTATCTGCGATGG (15 bp).
    std::string out;
    out += "AAAAAAAAAAAA";        // 12 bp padding
    out += "CGGTTCTT";             // canonical prefix
    out += "TTTTAAAACCCCGGGGAAATC"; // 21 bp spacer
    out += "GTCTATCTGCGATGG";      // canonical suffix
    out += "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    out.resize(100, 'A');
    return out;
}

// Build a synthetic 100 bp upstream window containing the duplicate motif.
std::string ChimdupUpstream() {
    std::string out;
    out += "AAAAAAAAAAAA";
    out += "GGGTTCTT";             // chimdup prefix (G→C swap vs canonical)
    out += "TTTTAAAACCCCGGGGAAATC";
    out += "AACTGTCCGCGAGG";       // chimdup suffix
    out += "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    out.resize(100, 'A');
    return out;
}

/// In-memory mock catalog. Looks up by haplotype-class tag or by
/// (target_id, pos). Tests should populate `entries_` directly in the body.
class MockCatalog final : public ISegDupCatalog {
public:
    struct Region {
        std::string target_id;
        std::uint64_t start{0};
        std::uint64_t end{0};
        std::string haplotype_class;
    };

    void AddEntry(CatalogEntryView v) {
        entries_[v.haplotype_class] = std::move(v);
    }

    void AddRegion(Region r) { regions_.push_back(std::move(r)); }

    std::optional<CatalogEntryView> LookupByPosition(
        std::string_view target_id, std::uint64_t pos) const override {
        for (const auto& r : regions_) {
            if (r.target_id == target_id && pos >= r.start && pos < r.end) {
                auto it = entries_.find(r.haplotype_class);
                if (it != entries_.end()) return it->second;
            }
        }
        return std::nullopt;
    }

    std::optional<CatalogEntryView> LookupByHaplotypeClass(
        std::string_view tag) const override {
        auto it = entries_.find(std::string{tag});
        if (it == entries_.end()) return std::nullopt;
        return it->second;
    }

private:
    std::unordered_map<std::string, CatalogEntryView> entries_;
    std::vector<Region> regions_;
};

/// Populate a chimdup-tandem catalog entry mirroring the curated JSON.
CatalogEntryView MakeChimdupEntry() {
    CatalogEntryView v;
    v.locus_id = "IGHG4_chimdup_tandem";
    v.haplotype_class = "IGHG4_chimdup_homozygous";
    v.mapping_primary.kmer_size = 25;
    v.mapping_primary.max_mismatch = 2;
    v.mapping_primary.include_flanking_bp = 200;
    v.mapping_primary.include_flanking_anchor = "CH1_start";

    PromoterSignature sig;
    sig.window_relative_to_anchor = "-100..0";
    sig.anchor = "IGHG4_CH1_start";
    sig.canonical_motif = std::string{kCanonicalMotif};
    sig.duplicate_motif = std::string{kDuplicateMotif};
    v.promoter_signature = std::move(sig);

    return v;
}

/// Compose a read: <upstream(100 bp)> + <CH1(50 bp)>. The anchor (CH1 start)
/// sits at offset 100. Returns the read seq and the anchor offset.
struct Read {
    std::string seq;
    std::uint32_t anchor_offset;
};
Read MakeRead(std::string_view upstream) {
    Read r;
    r.seq.assign(upstream);
    r.seq.append(kSyntheticCh1);
    r.anchor_offset = static_cast<std::uint32_t>(upstream.size());
    return r;
}

}  // namespace

// ---------------------------------------------------------------------------
// Low-level utility tests.
// ---------------------------------------------------------------------------

TEST(UpstreamAnchorUtil, ParseRelativeWindowNominal) {
    auto p = ParseRelativeWindow("-100..0");
    ASSERT_TRUE(p.has_value());
    EXPECT_EQ(p->first, -100);
    EXPECT_EQ(p->second, 0);
}

TEST(UpstreamAnchorUtil, ParseRelativeWindowSwapsOrder) {
    auto p = ParseRelativeWindow("50..-50");  // misordered → swapped
    ASSERT_TRUE(p.has_value());
    EXPECT_LE(p->first, p->second);
}

TEST(UpstreamAnchorUtil, ParseRelativeWindowRejectsGarbage) {
    EXPECT_FALSE(ParseRelativeWindow("").has_value());
    EXPECT_FALSE(ParseRelativeWindow("foo..bar").has_value());
    EXPECT_FALSE(ParseRelativeWindow("100").has_value());
}

TEST(UpstreamAnchorUtil, MotifEllipsisOrderedSubseq) {
    // Canonical motif against canonical window → match.
    EXPECT_TRUE(
        MotifMatchesWindow(kCanonicalMotif, CanonicalUpstream()));
    // Canonical motif against chimdup window → no match (prefix CGGTTCTT
    // is absent; chimdup has GGGTTCTT).
    EXPECT_FALSE(
        MotifMatchesWindow(kCanonicalMotif, ChimdupUpstream()));
    // Duplicate motif against chimdup window → match.
    EXPECT_TRUE(
        MotifMatchesWindow(kDuplicateMotif, ChimdupUpstream()));
}

TEST(UpstreamAnchorUtil, MotifMatchesNoEllipsis) {
    EXPECT_TRUE(MotifMatchesWindow("ACGTACGT", "TTACGTACGTAA"));
    EXPECT_FALSE(MotifMatchesWindow("ACGTACGT", "TTACGCACGTAA"));
}

TEST(UpstreamAnchorUtil, MotifEmptyInputs) {
    EXPECT_FALSE(MotifMatchesWindow("", "ACGT"));
    EXPECT_FALSE(MotifMatchesWindow("ACGT", ""));
}

TEST(UpstreamAnchorUtil, FlankingKmerMismatchesIdentical) {
    std::string s = "ACGTACGTACGT";
    EXPECT_EQ(0u, FlankingKmerMismatches(s, s, /*k=*/5, /*mm_cap=*/2));
}

TEST(UpstreamAnchorUtil, FlankingKmerMismatchesDivergent) {
    std::string a = "ACGTACGTACGT";
    std::string b = "TGCATGCATGCA";  // fully different
    EXPECT_GT(FlankingKmerMismatches(a, b, 5, 2), 0u);
}

TEST(UpstreamAnchorUtil, FlankingKmerMismatchesCaseInsensitive) {
    std::string a = "acgtacgtACGT";
    std::string b = "ACGTACGTACGT";
    EXPECT_EQ(0u, FlankingKmerMismatches(a, b, 5, 2));
}

TEST(UpstreamAnchorUtil, ToStringTokens) {
    EXPECT_EQ(ToString(PromoterMotifMatch::None),      "none");
    EXPECT_EQ(ToString(PromoterMotifMatch::Canonical), "canonical");
    EXPECT_EQ(ToString(PromoterMotifMatch::Chimdup),   "chimdup");
    EXPECT_EQ(ToString(PromoterMotifMatch::Mixed),     "mixed");
}

// ---------------------------------------------------------------------------
// End-to-end classification.
// ---------------------------------------------------------------------------

class UpstreamAnchorE2E : public ::testing::Test {
protected:
    MockCatalog catalog_;
    void SetUp() override {
        catalog_.AddEntry(MakeChimdupEntry());
        catalog_.AddRegion({"chr14", 105'625'000ULL, 105'650'000ULL,
                            "IGHG4_chimdup_homozygous"});
    }
};

TEST_F(UpstreamAnchorE2E, CanonicalReadCallsCanonicalMotif) {
    auto rd = MakeRead(CanonicalUpstream());
    ReadInput in;
    in.read_seq = rd.seq;
    in.target_id = "chr14";
    in.mapping_start = 105'626'000ULL;
    in.upstream_offset_in_read = rd.anchor_offset;

    auto res = Classify(catalog_, in);

    EXPECT_EQ(res.haplotype_class_call, "IGHG4_chimdup_homozygous");
    EXPECT_TRUE(res.promoter_test_ran);
    EXPECT_TRUE(res.flank_window_used);
    EXPECT_EQ(res.promoter_motif_match, PromoterMotifMatch::Canonical);
}

TEST_F(UpstreamAnchorE2E, ChimdupReadCallsChimdupMotif) {
    auto rd = MakeRead(ChimdupUpstream());
    ReadInput in;
    in.read_seq = rd.seq;
    in.target_id = "chr14";
    in.mapping_start = 105'626'000ULL;
    in.upstream_offset_in_read = rd.anchor_offset;

    auto res = Classify(catalog_, in);

    EXPECT_EQ(res.haplotype_class_call, "IGHG4_chimdup_homozygous");
    EXPECT_TRUE(res.promoter_test_ran);
    EXPECT_EQ(res.promoter_motif_match, PromoterMotifMatch::Chimdup);
    // Canonical reference vs chimdup window → expect non-zero mismatches.
    EXPECT_GT(res.flanking_kmer_mismatches, 0u);
}

TEST_F(UpstreamAnchorE2E, NoFlankReadHasZeroMismatchesAndNoMotif) {
    // Read = CH1 only, anchor at offset 0 → upstream window is empty.
    Read rd;
    rd.seq.assign(kSyntheticCh1);
    rd.anchor_offset = 0;

    ReadInput in;
    in.read_seq = rd.seq;
    in.target_id = "chr14";
    in.mapping_start = 105'626'000ULL;
    in.upstream_offset_in_read = rd.anchor_offset;

    auto res = Classify(catalog_, in);

    EXPECT_EQ(res.haplotype_class_call, "IGHG4_chimdup_homozygous");
    EXPECT_TRUE(res.promoter_test_ran);
    EXPECT_FALSE(res.flank_window_used);  // no upstream bases available
    EXPECT_EQ(res.promoter_motif_match, PromoterMotifMatch::None);
    EXPECT_EQ(res.flanking_kmer_mismatches, 0u);
}

TEST_F(UpstreamAnchorE2E, UnknownPositionReturnsEmptyCall) {
    // Outside the catalog's region.
    ReadInput in;
    in.read_seq = kSyntheticCh1;
    in.target_id = "chr1";
    in.mapping_start = 100'000ULL;
    in.upstream_offset_in_read = 0;

    auto res = Classify(catalog_, in);

    EXPECT_TRUE(res.haplotype_class_call.empty());
    EXPECT_EQ(res.CallOrUnknown(), "unknown");
    EXPECT_FALSE(res.promoter_test_ran);
}

TEST_F(UpstreamAnchorE2E, HypothesisOverridesPositionLookup) {
    auto rd = MakeRead(ChimdupUpstream());
    ReadInput in;
    in.read_seq = rd.seq;
    in.target_id = "chr1";        // wrong target
    in.mapping_start = 100'000ULL; // wrong position
    in.upstream_offset_in_read = rd.anchor_offset;

    ClassifyOptions opts;
    opts.hypothesis = "IGHG4_chimdup_homozygous";

    auto res = Classify(catalog_, in, opts);

    EXPECT_EQ(res.haplotype_class_call, "IGHG4_chimdup_homozygous");
    EXPECT_EQ(res.promoter_motif_match, PromoterMotifMatch::Chimdup);
}

TEST_F(UpstreamAnchorE2E, MixedReadCallsMixed) {
    // Construct an upstream containing BOTH motifs (e.g. chimeric junction
    // read that spans canonical and chimdup signatures). Both motifs must
    // fit in the catalog's "-100..0" window: keep the combined block under
    // 100 bp and place the anchor right after it.
    std::string mixed;
    mixed += "CGGTTCTT";          // 8  bp — canonical prefix
    mixed += "AAAAAAAA";          // 8  bp — spacer
    mixed += "GTCTATCTGCGATGG";   // 15 bp — canonical suffix
    mixed += "TTTT";              // 4  bp
    mixed += "GGGTTCTT";          // 8  bp — chimdup prefix
    mixed += "AAAAAAAA";          // 8  bp — spacer
    mixed += "AACTGTCCGCGAGG";    // 14 bp — chimdup suffix
    // total = 65 bp; well within 100 bp window
    ASSERT_LE(mixed.size(), 100u);

    auto rd = MakeRead(mixed);
    ReadInput in;
    in.read_seq = rd.seq;
    in.target_id = "chr14";
    in.mapping_start = 105'626'000ULL;
    in.upstream_offset_in_read = rd.anchor_offset;

    auto res = Classify(catalog_, in);
    EXPECT_EQ(res.promoter_motif_match, PromoterMotifMatch::Mixed);
}
