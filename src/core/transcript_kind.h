// LLmap — TranscriptKind: comprehensive RNA-class taxonomy.
//
// Standalone header so alignment_record.h doesn't carry the full
// classification noise. Used by:
//   - AlignmentRecord (this kind goes into the lossless output)
//   - anchor::AnchorRecord (annotates each catalog/computed anchor)
//   - annot::TranscriptKindClassifier (the heuristic-based classifier
//     that fills the kind field at anchor-load time)
//
// Design notes:
//
//   * uint16_t under-the-hood so we can grow the enum without ABI
//     breakage. The bestiary already has ~28 entries; new RNA classes
//     discovered post-2026 plug in here.
//
//   * The enum is deliberately fine-grained (Snorna_CDbox vs
//     Snorna_HacaBox, Snrna_Major vs Snrna_Minor, etc.) because the
//     mapping likelihood + expression priors differ across these
//     biological sub-classes. Lumping them as one "small_rna" bucket
//     would lose information that the Multi-Signal Fusion engine
//     actually consumes.
//
//   * `NovelUnclassified` plus `CustomKindTag` form the open-ended
//     discovery channel. When the heuristic classifier sees a
//     transcript that doesn't match any known kind, it falls through
//     to NovelUnclassified and attaches a free-form CustomKindTag with
//     the heuristic signature that triggered. These tags are tracked
//     in <out>.lossless.summary.json so an unrecognised RNA class that
//     keeps recurring is visible — candidate for enum promotion in a
//     later release.
//
// Memory crossrefs:
//   - [[ighg4_sgamma4_identical_in_tandem_dup]] — SterileGermline is
//     mandatory for IGH class-switch detection (Iγ4→Sγ4→Cγ4 reads).
//   - [[hprc_ighg4_tandem_dup_pervasive]] — circular RNAs of IGH origin
//     are documented in B-cell sub-populations; MappedCircular catches
//     them as a first-class status.

#pragma once

#include <cstdint>
#include <optional>
#include <string>

namespace llmap::core {

/// Transcript / RNA-class taxonomy. uint16_t for headroom; current
/// occupants fit in <40 slots but more keep being discovered.
enum class TranscriptKind : std::uint16_t {
    Unknown = 0,

    // --- protein-coding family ------------------------------------------
    MatureMrna,             ///< spliced, polyadenylated, 5'-capped mRNA
    PreMrna,                ///< unspliced or partially spliced precursor
    SterileGermline,        ///< I-promoter → S-region → C-gene; no V-D-J;
                            ///< prerequisite for IGH class-switch detection

    // --- circular ---------------------------------------------------------
    CircularRna,            ///< back-spliced circRNA (BSJ-derived)

    // --- long ncRNA -------------------------------------------------------
    Lncrna,                 ///< long noncoding (incl. lincRNA, intronic-lnc)

    // --- small RNAs (size-class driven, all matter for LLmap's k=21 path) -
    Snorna_CDbox,           ///< C/D box snoRNA — 2'-O-Me guide
    Snorna_HacaBox,         ///< H/ACA box snoRNA — pseudouridylation guide
    Scarna,                 ///< small Cajal body RNA
    Mirna,                  ///< pri/pre/mature microRNA
    Pirna,                  ///< PIWI-interacting RNA (~26-31 nt)
    Sirna,                  ///< small interfering RNA (~21-23 nt)
    Snrna_Major,            ///< U1, U2, U4, U5, U6 (major spliceosome)
    Snrna_Minor,            ///< U11, U12, U4atac, U6atac (minor spliceosome)

    // --- structural / housekeeping ---------------------------------------
    Trna,                   ///< transfer RNA (cyto + mt)
    Rrna,                   ///< ribosomal (28S, 18S, 5.8S, 5S + pre-rRNA)
    Yrna,                   ///< Y-RNA (Ro60-associated)
    Vaultrna,               ///< vault RNA (vtRNA1-1 etc.)
    Srprna_7sl,             ///< 7SL — signal recognition particle
    Rna_7sk,                ///< 7SK — RNA pol II elongation regulator
    TercTerra,              ///< TERC / TERRA (telomere RNA component)
    Rmrp_Rnasep,            ///< RMRP, RNase-P, RNase-MRP

    // --- regulatory / pervasive ------------------------------------------
    Erna,                   ///< enhancer RNA
    TirnaParna,             ///< transcription initiation / promoter-assoc
    Antisense,              ///< natural antisense transcript (NAT)
    Intergenic,             ///< intergenic novel transcribed units

    // --- derived / repetitive / exogenous --------------------------------
    Mitochondrial,          ///< mt-mRNA, mt-rRNA, mt-tRNA, D-loop transcripts
    RepeatDerived,          ///< LINE/SINE/Alu/HERV-derived
    Viral,                  ///< exogenous viral RNA or viRNA

    // --- composite -------------------------------------------------------
    Fusion,                 ///< fusion transcript across two gene loci

    // --- open-ended discovery channel ------------------------------------
    NovelUnclassified,      ///< detected, doesn't match any known kind
};

/// Free-form extension tag for RNA classes the enum does not yet cover.
///
/// When `kind == NovelUnclassified`, the classifier attaches one of these
/// so downstream tooling can group anonymous discoveries by the heuristic
/// signature that triggered them. Frequent labels are promotion candidates
/// for a future enum slot.
struct CustomKindTag {
    /// Short stable label, e.g. "novel_short_120nt", "tipi_rna_like".
    /// Convention: snake_case ASCII, no spaces. Length cap soft-set at 64.
    std::string label;

    /// Which heuristic produced this label (helps when re-running the
    /// classifier with a different rule-set and comparing distributions).
    /// Example: "length_band:100-300 + no_polyA + low_gc".
    std::string reason_signature;
};

/// Convenience predicate — true for any of the "Mapped*" status families
/// in alignment_record.h; the kind matters but the status carries the
/// authoritative class membership.
constexpr bool IsCodingFamily(TranscriptKind k) noexcept {
    return k == TranscriptKind::MatureMrna
        || k == TranscriptKind::PreMrna
        || k == TranscriptKind::SterileGermline;
}

constexpr bool IsSmallRna(TranscriptKind k) noexcept {
    return k == TranscriptKind::Mirna
        || k == TranscriptKind::Pirna
        || k == TranscriptKind::Sirna
        || k == TranscriptKind::Snorna_CDbox
        || k == TranscriptKind::Snorna_HacaBox
        || k == TranscriptKind::Scarna
        || k == TranscriptKind::Snrna_Major
        || k == TranscriptKind::Snrna_Minor;
}

constexpr bool IsStructuralRna(TranscriptKind k) noexcept {
    return k == TranscriptKind::Trna
        || k == TranscriptKind::Rrna
        || k == TranscriptKind::Yrna
        || k == TranscriptKind::Vaultrna
        || k == TranscriptKind::Srprna_7sl
        || k == TranscriptKind::Rna_7sk
        || k == TranscriptKind::TercTerra
        || k == TranscriptKind::Rmrp_Rnasep;
}

/// Format helper for telemetry / output. Stable strings; never rename
/// existing entries (downstream tools depend on the exact spelling).
const char* TranscriptKindName(TranscriptKind k) noexcept;

/// Parse the string form back to the enum. Returns std::nullopt on
/// unknown labels (caller may then construct a CustomKindTag).
std::optional<TranscriptKind> ParseTranscriptKind(std::string_view s) noexcept;

}  // namespace llmap::core
