// LLmap — ModificationDb: orthogonal RNA-modification annotation layer.
//
// >170 distinct RNA modifications are catalogued in MODOMICS 2025;
// of these, ~18 are routinely detectable from long-read sequencing
// data (Nanopore RNA004 + Dorado, PacBio kineticsTools). LLmap does
// NOT call modifications from raw signal — that's modkit / m6Anet /
// EpiNano / dorado-rna territory. Instead, LLmap **imports + preserves**
// their output and uses the modification pattern as an extra
// disambiguating signal in the Multi-Signal Fusion engine (Block 4.5).
//
// Two database axes are loaded here:
//
//   1. Sequence-context priors. P(m6A | DRACH) ≈ 0.30 in mRNA per
//      Pratanwanich 2021. Used by the Fusion engine to soften
//      L_sequence at known mod-prone positions (a mismatch in DRACH
//      context could be modification, not error).
//
//   2. Site-level annotations from REPIC (m6A), REDIportal (A→I),
//      MODOMICS (Ψ, m5C, m1A, etc.). Per-transcript, per-position
//      site lists. Used to disambiguate sequence-identical paralogs
//      via their tissue-specific mod patterns:
//        canonical IGHG4 has m6A at pos N (REPIC: 7/8 B-cell samples)
//        ChimDup   IGHG4 has m6A at pos N (REPIC: 0/8 — never reported)
//        → a read with m6A-call at pos N biases the Wave-Collapse
//          posterior toward canonical
//
// Threading: read-only after load; multiple consumers may Lookup concurrently.

#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::rnamod {

// ===========================================================================
// ModificationKind — uint16_t for extensibility, ~30 documented entries.
// ===========================================================================

enum class ModificationKind : std::uint16_t {
    Unknown = 0,
    // mRNA-dominant
    M6A,             ///< N6-methyladenosine (most abundant mRNA mod)
    M6Am,            ///< N6,2'-O-dimethyladenosine (5'-cap +1)
    M5C,             ///< 5-methylcytosine
    Hm5C,            ///< 5-hydroxymethyl-C
    M1A,             ///< 1-methyladenosine
    Pseudouridine,   ///< Ψ
    Inosine,         ///< A→I editing via ADAR
    M7G,             ///< 7-methylguanosine (5'-cap)
    M3C,             ///< 3-methylcytosine
    M2_O_Methyl,     ///< 2'-O-methylation (Nm)
    Ac4C,            ///< 4-acetylcytidine
    M6_2A,           ///< N6,N6-dimethyladenosine
    // tRNA / snoRNA / rRNA
    M1G,
    M22G,
    M5U,
    T6A,
    Mcm5U,
    Queosine,
    // Editing events (base-changes, not stable mods)
    A2I_Editing,
    C2U_Editing,     ///< APOBEC / AID — central to IGH class-switch
    // Cap variants
    Cap0,
    Cap1,
    Cap2,
    // Poly-A variants
    Polya_Standard,
    Polya_Short,
    Polya_Modified,  ///< contains m6A in poly-A tail
    Oligo_U_Tail,
    // Discovery channel
    NovelUnclassifiedMod,
};

const char* ModificationKindName(ModificationKind k) noexcept;
std::optional<ModificationKind> ParseModificationKind(std::string_view s) noexcept;

// ===========================================================================
// ModificationCall — one observed modification on one read.
//
// These come from imported tool output (modkit BED, m6Anet TSV,
// PacBio kineticsTools GFF, etc.) and ride along on the AlignmentRecord
// of the read that carries them.
// ===========================================================================

struct ModificationCall {
    ModificationKind kind{ModificationKind::Unknown};
    std::uint32_t pos_in_read{0};        ///< 0-based offset in raw read
    std::optional<std::uint32_t>
        pos_in_transcript;                ///< if mappable to transcript coord
    float confidence{0.0f};               ///< tool-reported confidence [0,1]
    /// Source tool:
    ///   'M' = modkit (Nanopore RNA004)
    ///   'D' = dorado-rna integrated
    ///   'E' = Eligos2
    ///   'A' = m6Anet
    ///   'P' = PacBio kineticsTools
    ///   'O' = other / unknown
    char source_tool{'O'};
    /// For NovelUnclassifiedMod: free-form label.
    std::optional<std::string> custom_label;
};

// ===========================================================================
// ModificationDb — site-level annotations from REPIC / REDIportal / MODOMICS.
// ===========================================================================

class ModificationDb {
public:
    /// Built-in sequence-context priors. Always succeeds.
    void LoadDefaults();

    /// REPIC m6A site database — BED5+: chrom, pos, pos+1, ref, alt, sample_count
    bool LoadREPIC(const std::filesystem::path& bed);

    /// REDIportal A→I editing sites — TSV with chrom + pos + transcript_id.
    bool LoadREDIportal(const std::filesystem::path& tsv);

    /// MODOMICS extended catalog (~170 modifications). Best-effort TSV
    /// parser; unknown columns ignored.
    bool LoadModomicsExtended(const std::filesystem::path& tsv);

    /// Sites known at (transcript_id, pos).
    [[nodiscard]] std::vector<ModificationKind>
    KnownAt(std::string_view transcript_id, std::uint32_t pos) const;

    /// Sequence-context prior. Returns P(mod | flanking 5-mer context)
    /// per published frequencies. 0 ⇒ unknown / context not recognised.
    /// Example: PriorProbability(M6A, "GGACT") returns ≈0.30 (DRACH match).
    [[nodiscard]] float PriorProbability(ModificationKind kind,
                                          std::string_view context_5bp) const;

    /// Telemetry — number of distinct sites loaded.
    [[nodiscard]] std::size_t TotalSites() const noexcept;

private:
    // (transcript_id, pos) → list of known mods
    struct SiteKey {
        std::string transcript_id;
        std::uint32_t pos;
        bool operator==(const SiteKey& o) const noexcept {
            return transcript_id == o.transcript_id && pos == o.pos;
        }
    };
    struct SiteKeyHash {
        std::size_t operator()(const SiteKey& k) const noexcept {
            auto h1 = std::hash<std::string>{}(k.transcript_id);
            return h1 ^ (static_cast<std::size_t>(k.pos) * 0x9E3779B97F4A7C15ULL);
        }
    };

    std::unordered_map<SiteKey, std::vector<ModificationKind>, SiteKeyHash>
        sites_;

    // (kind, context_5bp) → prior probability
    std::unordered_map<std::string, float> priors_;

    bool defaults_loaded_{false};
};

}  // namespace llmap::rnamod
