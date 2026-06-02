// LLmap — AnchorSource name/parse helpers.
//
// Stable wire-level strings; never rename existing entries (downstream
// tools depend on the exact spelling).

#include "anchor/anchor_record.h"

#include <array>
#include <string_view>

namespace llmap::anchor {

namespace {

struct SourceLabel {
    AnchorSource src;
    std::string_view name;
};

constexpr std::array<SourceLabel, 12> kLabels = {{
    {AnchorSource::Unknown,           "unknown"},
    {AnchorSource::Gencode,           "gencode"},
    {AnchorSource::Mane,              "mane"},
    {AnchorSource::Refseq,            "refseq"},
    {AnchorSource::Imgt_GeneDb,       "imgt_genedb"},
    {AnchorSource::Imgt_Hla,          "imgt_hla"},
    {AnchorSource::Fantom5_Cat,       "fantom5_cat"},
    {AnchorSource::ChessDb,           "chessdb"},
    {AnchorSource::Pangenome_PerHap,  "pangenome_perhap"},
    {AnchorSource::Branch_Bubble,     "branch_bubble"},
    {AnchorSource::Computed_Cluster,  "computed_cluster"},
    {AnchorSource::Custom,            "custom"},
}};

}  // namespace

const char* AnchorSourceName(AnchorSource s) noexcept {
    for (const auto& [src, name] : kLabels) {
        if (src == s) return name.data();
    }
    return "unknown";
}

std::optional<AnchorSource> ParseAnchorSource(std::string_view s) noexcept {
    for (const auto& [src, name] : kLabels) {
        if (name == s) return src;
    }
    return std::nullopt;
}

}  // namespace llmap::anchor
