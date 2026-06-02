// LLmap — TranscriptKindClassifier implementation.
//
// Hand-coded heuristics, deliberately simple and fast. Each detector
// returns bool; Classify() composes them in an ordered cascade with
// the most specific rules winning. Open-ended fallback assigns a
// CustomKindTag so the discovery channel stays observable.

#include "annot/transcript_kind_classifier.h"

#include "annot/splice_site_db.h"

#include <algorithm>
#include <cctype>
#include <string>

namespace llmap::annot {

namespace {

bool HasTag(const anchor::AnchorRecord& a, std::string_view tag) {
    return std::any_of(a.tags.begin(), a.tags.end(),
                        [&](const std::string& t) { return t == tag; });
}

bool HasTagPrefix(const anchor::AnchorRecord& a, std::string_view prefix) {
    return std::any_of(a.tags.begin(), a.tags.end(),
                        [&](const std::string& t) {
                            return t.size() >= prefix.size()
                                && std::equal(prefix.begin(), prefix.end(), t.begin());
                        });
}

bool StartsWith(std::string_view s, std::string_view prefix) {
    return s.size() >= prefix.size()
        && std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::string ToUpper(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) out.push_back(std::toupper(static_cast<unsigned char>(c)));
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Specialised detectors
// ---------------------------------------------------------------------------

bool TranscriptKindClassifier::LooksLikeSterileGermlineTranscript(
    const anchor::AnchorRecord& a) const {
    // Sterile germline transcripts are unique to IGH/IGK/IGL/TR loci:
    // a single transcript spans I-promoter → S-region → C-gene with
    // NO V-D-J upstream.
    const bool ig_locus = StartsWith(a.host_gene_id, "IGH")
                          || StartsWith(a.host_gene_id, "IGK")
                          || StartsWith(a.host_gene_id, "IGL")
                          || StartsWith(a.host_gene_id, "TR");
    if (!ig_locus) return false;

    // Switch-region tag is the strong signal.
    const bool has_switch = HasTagPrefix(a, "switch_region")
                             || HasTagPrefix(a, "Sgamma")
                             || HasTagPrefix(a, "Smu")
                             || HasTagPrefix(a, "Salpha")
                             || HasTagPrefix(a, "Sepsilon");
    if (!has_switch) return false;

    // Make sure we're not on a V-D-J segment.
    const bool has_vdj = HasTagPrefix(a, "V_gene")
                          || HasTagPrefix(a, "D_gene")
                          || HasTagPrefix(a, "J_gene")
                          || HasTagPrefix(a, "biotype:IG_V")
                          || HasTagPrefix(a, "biotype:IG_D")
                          || HasTagPrefix(a, "biotype:IG_J");
    return !has_vdj;
}

bool TranscriptKindClassifier::LooksLikeCircularRna(
    const anchor::AnchorRecord& a,
    const SpliceSiteDb& /*splice*/) const {
    // Back-splice junction is the structural marker. The
    // SpliceSiteDb-aware version (caller passes splice in here) could
    // additionally validate motif consistency, but is_circular() on the
    // anchor already checks for spliceosome_class == 3.
    if (a.is_circular()) return true;
    // Tag-based fallback: callers may have already marked it.
    return HasTag(a, "circular_rna") || HasTagPrefix(a, "biotype:circRNA");
}

bool TranscriptKindClassifier::LooksLikePreMrna(
    const anchor::AnchorRecord& a) const {
    return HasTag(a, "pre_mrna")
        || HasTag(a, "intron_retained")
        || HasTagPrefix(a, "biotype:retained_intron")
        || HasTagPrefix(a, "biotype:nonsense_mediated_decay");
}

bool TranscriptKindClassifier::LooksLikePirna(
    const anchor::AnchorRecord& a) const {
    // No exon boundaries (small RNAs are unsplicedish) AND length in
    // piRNA range. ref_end - ref_start used when sequence is empty.
    if (a.is_spliced()) return false;
    std::size_t len = a.sequence.size();
    if (len == 0 && a.has_genomic_coords()) {
        len = static_cast<std::size_t>(*a.ref_end - *a.ref_start);
    }
    return len >= 24 && len <= 33;
}

bool TranscriptKindClassifier::LooksLikeSirna(
    const anchor::AnchorRecord& a) const {
    if (a.is_spliced()) return false;
    std::size_t len = a.sequence.size();
    if (len == 0 && a.has_genomic_coords()) {
        len = static_cast<std::size_t>(*a.ref_end - *a.ref_start);
    }
    return len >= 20 && len <= 24;
}

core::TranscriptKind TranscriptKindClassifier::RefineSnornaSubclass(
    const anchor::AnchorRecord& a) const {
    if (a.sequence.empty()) return core::TranscriptKind::Unknown;
    const std::string up = ToUpper(a.sequence);

    // C/D box hallmark: 5' RUGAUGA  ... 3' CUGA terminator.
    const bool has_cdbox_5p = (up.find("AUGAUGA") != std::string::npos
                                || up.find("GUGAUGA") != std::string::npos
                                || up.find("ATGATGA") != std::string::npos
                                || up.find("GTGATGA") != std::string::npos);
    const bool has_cdbox_3p = (up.find("CUGA") != std::string::npos
                                || up.find("CTGA") != std::string::npos);
    if (has_cdbox_5p && has_cdbox_3p && up.size() >= 50 && up.size() <= 200) {
        return core::TranscriptKind::Snorna_CDbox;
    }

    // H/ACA box hallmark: ANANNA midway + ACA near 3' end.
    const bool has_anana = (up.find("ANANNA") != std::string::npos
                             || up.find("AAANNA") != std::string::npos
                             || up.find("ACANNA") != std::string::npos);
    const bool has_aca_3p =
        up.size() >= 3
        && (up.substr(up.size() - 3) == "ACA"
            || up.substr(up.size() - 4, 3) == "ACA"
            || up.substr(up.size() - 5, 3) == "ACA");
    if (has_anana && has_aca_3p && up.size() >= 100 && up.size() <= 220) {
        return core::TranscriptKind::Snorna_HacaBox;
    }

    return core::TranscriptKind::Unknown;
}

// ---------------------------------------------------------------------------
// Open-ended fallback
// ---------------------------------------------------------------------------

std::pair<core::TranscriptKind, core::CustomKindTag>
TranscriptKindClassifier::ClassifyOpenEnded(
    const anchor::AnchorRecord& a) const {
    core::CustomKindTag tag;
    std::size_t len = a.sequence.size();
    if (len == 0 && a.has_genomic_coords()) {
        len = static_cast<std::size_t>(*a.ref_end - *a.ref_start);
    }

    // Label encodes length-band so frequent unknowns are grouped.
    std::string band;
    if      (len < 50)    band = "short_lt50";
    else if (len < 200)   band = "short_50_200";
    else if (len < 500)   band = "medium_200_500";
    else if (len < 2000)  band = "medium_500_2k";
    else if (len < 10000) band = "long_2k_10k";
    else                  band = "very_long_10k+";

    tag.label = "novel_" + band;
    tag.reason_signature = "no_known_pattern;length_band=" + band
        + ";spliced=" + (a.is_spliced() ? "yes" : "no")
        + ";has_coords=" + (a.has_genomic_coords() ? "yes" : "no");

    return {core::TranscriptKind::NovelUnclassified, std::move(tag)};
}

// ---------------------------------------------------------------------------
// Top-level Classify — ordered cascade.
// ---------------------------------------------------------------------------

std::pair<core::TranscriptKind, std::optional<core::CustomKindTag>>
TranscriptKindClassifier::Classify(
    const anchor::AnchorRecord& anchor,
    const SpliceSiteDb& splice,
    const TissueContext& /*tissue*/) const {
    using K = core::TranscriptKind;

    // 1. Most specific structural class wins first.
    if (LooksLikeCircularRna(anchor, splice)) {
        return {K::CircularRna, std::nullopt};
    }
    if (LooksLikeSterileGermlineTranscript(anchor)) {
        return {K::SterileGermline, std::nullopt};
    }
    if (LooksLikePreMrna(anchor)) {
        return {K::PreMrna, std::nullopt};
    }

    // 2. Short-RNA structural fingerprint.
    if (LooksLikeSirna(anchor)) return {K::Sirna,  std::nullopt};
    if (LooksLikePirna(anchor)) return {K::Pirna,  std::nullopt};

    // 3. snoRNA subclass refinement when the loader said Snorna_CDbox
    //    by default.
    if (anchor.kind == K::Snorna_CDbox || anchor.kind == K::Snorna_HacaBox) {
        auto refined = RefineSnornaSubclass(anchor);
        if (refined != K::Unknown) return {refined, std::nullopt};
    }

    // 4. Loader's coarse kind, if any.
    if (anchor.kind != K::Unknown) {
        return {anchor.kind, std::nullopt};
    }

    // 5. Open-ended discovery.
    auto [kind, tag] = ClassifyOpenEnded(anchor);
    return {kind, std::move(tag)};
}

}  // namespace llmap::annot
