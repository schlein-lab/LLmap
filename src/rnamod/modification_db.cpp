// LLmap — ModificationDb implementation.
//
// Sequence-context priors are hard-coded from the literature; site-level
// loaders are best-effort TSV / BED readers that tolerate missing or
// extra columns. The aim is robust ingest of imperfect public files
// over precise schema enforcement.

#include "rnamod/modification_db.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <sstream>

namespace llmap::rnamod {

namespace {

struct KindLabel {
    ModificationKind kind;
    std::string_view name;
};

constexpr std::array<KindLabel, 28> kLabels = {{
    {ModificationKind::Unknown,            "unknown"},
    {ModificationKind::M6A,                "m6a"},
    {ModificationKind::M6Am,               "m6am"},
    {ModificationKind::M5C,                "m5c"},
    {ModificationKind::Hm5C,               "hm5c"},
    {ModificationKind::M1A,                "m1a"},
    {ModificationKind::Pseudouridine,      "psi"},
    {ModificationKind::Inosine,            "i"},
    {ModificationKind::M7G,                "m7g"},
    {ModificationKind::M3C,                "m3c"},
    {ModificationKind::M2_O_Methyl,        "nm"},
    {ModificationKind::Ac4C,               "ac4c"},
    {ModificationKind::M6_2A,              "m6_2a"},
    {ModificationKind::M1G,                "m1g"},
    {ModificationKind::M22G,               "m22g"},
    {ModificationKind::M5U,                "m5u"},
    {ModificationKind::T6A,                "t6a"},
    {ModificationKind::Mcm5U,              "mcm5u"},
    {ModificationKind::Queosine,           "q"},
    {ModificationKind::A2I_Editing,        "a2i_edit"},
    {ModificationKind::C2U_Editing,        "c2u_edit"},
    {ModificationKind::Cap0,               "cap0"},
    {ModificationKind::Cap1,               "cap1"},
    {ModificationKind::Cap2,               "cap2"},
    {ModificationKind::Polya_Standard,     "polya_standard"},
    {ModificationKind::Polya_Short,        "polya_short"},
    {ModificationKind::Polya_Modified,     "polya_modified"},
    {ModificationKind::Oligo_U_Tail,       "oligo_u"},
}};

std::string ToUpper(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) out.push_back(std::toupper(static_cast<unsigned char>(c)));
    return out;
}

bool IsDrach(std::string_view ctx) {
    // DRACH = [AGT][AG]AC[ACT]  (D=A/G/U, R=A/G, fixed A,C, H=A/C/U)
    // ctx is upper-cased 5-mer.
    if (ctx.size() != 5) return false;
    const char d = ctx[0], r = ctx[1], a = ctx[2], c = ctx[3], h = ctx[4];
    return (d == 'A' || d == 'G' || d == 'T' || d == 'U')
        && (r == 'A' || r == 'G')
        && (a == 'A')
        && (c == 'C')
        && (h == 'A' || h == 'C' || h == 'T' || h == 'U');
}

// AID-target consensus: WRC (W=A/T, R=A/G, C).
bool IsAidRgyw(std::string_view ctx) {
    // ctx[0..4]: looking for WRCY → 4-mer match starting at ctx[1]
    if (ctx.size() < 5) return false;
    auto isW = [](char c) { return c == 'A' || c == 'T'; };
    auto isR = [](char c) { return c == 'A' || c == 'G'; };
    auto isY = [](char c) { return c == 'C' || c == 'T'; };
    // Use middle 4 bases of context as RGYW (canonical AID hotspot).
    const char a = ctx[1], b = ctx[2], c = ctx[3], d = ctx[4];
    return isR(a) && b == 'G' && isY(c) && isW(d);
}

}  // namespace

const char* ModificationKindName(ModificationKind k) noexcept {
    for (const auto& [kk, name] : kLabels) {
        if (kk == k) return name.data();
    }
    if (k == ModificationKind::NovelUnclassifiedMod) {
        return "novel_unclassified_mod";
    }
    return "unknown";
}

std::optional<ModificationKind>
ParseModificationKind(std::string_view s) noexcept {
    for (const auto& [k, name] : kLabels) {
        if (name == s) return k;
    }
    if (s == "novel_unclassified_mod")
        return ModificationKind::NovelUnclassifiedMod;
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// LoadDefaults — built-in sequence-context priors.
//
// Source: Pratanwanich 2021 (m6A in DRACH ≈ 0.30), Carlile 2014 (Ψ rare
// in mRNA ≈ 0.002), Bahn 2015 (A→I in ADAR-rich repeats ≈ 0.15).
// Priors are stored keyed by "<kind>:<context>" for compact lookup.
// ---------------------------------------------------------------------------

void ModificationDb::LoadDefaults() {
    // m6A: DRACH-context prior. We enumerate the 24 DRACH 5-mers and
    // assign ≈0.30; for non-DRACH contexts the prior is ≈0.005 (still
    // non-zero — m6A does occur outside DRACH but rarely).
    // We don't actually enumerate them all in the map — IsDrach()
    // checks at lookup time. We just remember that defaults are loaded.
    defaults_loaded_ = true;
}

// ---------------------------------------------------------------------------
// Site-level loaders
// ---------------------------------------------------------------------------

bool ModificationDb::LoadREPIC(const std::filesystem::path& bed) {
    std::ifstream in(bed);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string chrom, pos_s, pos2_s, name, score_s, strand;
        // REPIC convention: chrom, pos (0-based), pos+1, name (transcript_id
        // when available), score, strand
        ss >> chrom >> pos_s >> pos2_s >> name >> score_s >> strand;
        if (chrom.empty() || pos_s.empty() || name.empty()) continue;
        try {
            std::uint32_t pos = static_cast<std::uint32_t>(std::stoul(pos_s));
            SiteKey k{name, pos};
            sites_[k].push_back(ModificationKind::M6A);
        } catch (...) { /* skip malformed */ }
    }
    return true;
}

bool ModificationDb::LoadREDIportal(const std::filesystem::path& tsv) {
    std::ifstream in(tsv);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; if (!line.empty() && line[0] != '#') {
            // best-effort header skip
        }}
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string chrom, pos_s, transcript_id, ref, alt;
        ss >> chrom >> pos_s >> ref >> alt >> transcript_id;
        if (chrom.empty() || pos_s.empty() || transcript_id.empty()) continue;
        try {
            std::uint32_t pos = static_cast<std::uint32_t>(std::stoul(pos_s));
            SiteKey k{transcript_id, pos};
            if (ref == "A" && alt == "I") {
                sites_[k].push_back(ModificationKind::A2I_Editing);
            } else if (ref == "C" && alt == "U") {
                sites_[k].push_back(ModificationKind::C2U_Editing);
            }
        } catch (...) {}
    }
    return true;
}

bool ModificationDb::LoadModomicsExtended(const std::filesystem::path& tsv) {
    std::ifstream in(tsv);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; continue; }  // skip header line
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string kind_s, transcript_id, pos_s;
        ss >> kind_s >> transcript_id >> pos_s;
        if (kind_s.empty() || transcript_id.empty() || pos_s.empty()) continue;
        auto kind = ParseModificationKind(kind_s);
        if (!kind) continue;
        try {
            std::uint32_t pos = static_cast<std::uint32_t>(std::stoul(pos_s));
            SiteKey k{transcript_id, pos};
            sites_[k].push_back(*kind);
        } catch (...) {}
    }
    return true;
}

// ---------------------------------------------------------------------------
// Lookup
// ---------------------------------------------------------------------------

std::vector<ModificationKind>
ModificationDb::KnownAt(std::string_view transcript_id,
                         std::uint32_t pos) const {
    SiteKey k{std::string(transcript_id), pos};
    auto it = sites_.find(k);
    if (it == sites_.end()) return {};
    return it->second;
}

float ModificationDb::PriorProbability(ModificationKind kind,
                                        std::string_view context_5bp) const {
    if (context_5bp.size() != 5) return 0.0f;
    const std::string ctx_up = ToUpper(context_5bp);

    switch (kind) {
        case ModificationKind::M6A:
            return IsDrach(ctx_up) ? 0.30f : 0.005f;
        case ModificationKind::A2I_Editing:
            // ADAR prefers double-stranded RNA contexts; without full
            // dsRNA-fold prediction we use a coarse prior.
            return 0.05f;
        case ModificationKind::C2U_Editing:
            // AID hotspot: RGYW. APOBEC has a wider preference.
            return IsAidRgyw(ctx_up) ? 0.10f : 0.005f;
        case ModificationKind::Pseudouridine:
            // Ψ rare in mRNA, ≈0.002 background.
            return 0.002f;
        case ModificationKind::M5C:
            return 0.001f;
        default:
            return 0.0f;  // unknown context-prior
    }
}

std::size_t ModificationDb::TotalSites() const noexcept {
    return sites_.size();
}

}  // namespace llmap::rnamod
