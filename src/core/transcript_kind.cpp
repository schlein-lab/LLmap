// LLmap — TranscriptKind name <-> enum mapping.
//
// Kept in a separate .cpp so the header (transcript_kind.h) stays
// lightweight enough to include from many places. The name table is
// stable: existing labels are part of the on-disk Parquet / BAM tag
// contract and must never be renamed.

#include "core/transcript_kind.h"

#include <array>
#include <string_view>

namespace llmap::core {

namespace {

// One entry per enum value. Position == enum integer.
// New kinds: append, never reorder.
struct KindLabel {
    TranscriptKind kind;
    std::string_view name;
};

constexpr std::array<KindLabel, 30> kLabels = {{
    {TranscriptKind::Unknown,           "unknown"},
    {TranscriptKind::MatureMrna,        "mature_mrna"},
    {TranscriptKind::PreMrna,           "pre_mrna"},
    {TranscriptKind::SterileGermline,   "sterile_germline"},
    {TranscriptKind::CircularRna,       "circular_rna"},
    {TranscriptKind::Lncrna,            "lncrna"},
    {TranscriptKind::Snorna_CDbox,      "snorna_cdbox"},
    {TranscriptKind::Snorna_HacaBox,    "snorna_hacabox"},
    {TranscriptKind::Scarna,            "scarna"},
    {TranscriptKind::Mirna,             "mirna"},
    {TranscriptKind::Pirna,             "pirna"},
    {TranscriptKind::Sirna,             "sirna"},
    {TranscriptKind::Snrna_Major,       "snrna_major"},
    {TranscriptKind::Snrna_Minor,       "snrna_minor"},
    {TranscriptKind::Trna,              "trna"},
    {TranscriptKind::Rrna,              "rrna"},
    {TranscriptKind::Yrna,              "yrna"},
    {TranscriptKind::Vaultrna,          "vaultrna"},
    {TranscriptKind::Srprna_7sl,        "srprna_7sl"},
    {TranscriptKind::Rna_7sk,           "rna_7sk"},
    {TranscriptKind::TercTerra,         "terc_terra"},
    {TranscriptKind::Rmrp_Rnasep,       "rmrp_rnasep"},
    {TranscriptKind::Erna,              "erna"},
    {TranscriptKind::TirnaParna,        "tirna_parna"},
    {TranscriptKind::Antisense,         "antisense"},
    {TranscriptKind::Intergenic,        "intergenic"},
    {TranscriptKind::Mitochondrial,     "mitochondrial"},
    {TranscriptKind::RepeatDerived,     "repeat_derived"},
    {TranscriptKind::Viral,             "viral"},
    {TranscriptKind::Fusion,            "fusion"},
}};
// NovelUnclassified intentionally not in the lookup table — we want
// it to round-trip via its label explicitly so callers don't accidentally
// stringify a NovelUnclassified-with-no-CustomKindTag.

}  // namespace

const char* TranscriptKindName(TranscriptKind k) noexcept {
    for (const auto& [kind, name] : kLabels) {
        if (kind == k) return name.data();
    }
    if (k == TranscriptKind::NovelUnclassified) return "novel_unclassified";
    return "unknown";
}

std::optional<TranscriptKind>
ParseTranscriptKind(std::string_view s) noexcept {
    for (const auto& [kind, name] : kLabels) {
        if (name == s) return kind;
    }
    if (s == "novel_unclassified") return TranscriptKind::NovelUnclassified;
    return std::nullopt;
}

}  // namespace llmap::core
