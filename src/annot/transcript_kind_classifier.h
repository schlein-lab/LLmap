// LLmap — TranscriptKindClassifier: heuristic refinement of anchor kind.
//
// The GENCODE loader assigns a coarse TranscriptKind from each entry's
// biotype string. That works for ~95 % of cases but misses several
// biologically important classes that the biotype field doesn't
// distinguish:
//
//   - Sterile germline transcripts (IGH Iγ→Sγ→Cγ, "transcript_type"
//     in GENCODE is just "protein_coding") — needs structural detection
//     based on the absence of V-D-J upstream + presence of an I-promoter
//     anchor + Switch-region.
//
//   - C/D-box vs H/ACA-box snoRNA — GENCODE labels both as "snoRNA".
//     The C/D box motif (RUGAUGA + CUGA terminator) and H/ACA box
//     (ANANNA + ACA) are sequence-detectable; we run the check at
//     load time.
//
//   - circRNA — only catchable when ExonBoundary records have
//     spliceosome_class == 3 (back-splice). Done as a structural
//     post-pass.
//
//   - pre-mRNA detection — when an anchor carries retained introns
//     (visible from coverage in the read, not from the catalog
//     directly), the classifier tags the kind. This happens at
//     read-time, not load-time, but the classifier owns the logic.
//
// The classifier is deliberately separated from the loader so different
// catalog sources (GENCODE, MANE, IMGT) can share the same refinement
// pipeline.

#pragma once

#include "anchor/anchor_record.h"
#include "core/transcript_kind.h"

#include <optional>
#include <string_view>
#include <utility>

namespace llmap::annot {

class SpliceSiteDb;  // fwd-decl to avoid include cycle

/// Tissue context — used as a sanity-check input when the classifier
/// has to decide between rare kinds (e.g. sterile-germline transcripts
/// in B-cells vs spurious IGH coverage in other tissues).
struct TissueContext {
    std::string label;        ///< "lymph", "pbmc", "b_cell", "lcl", ""
    std::string sample_id;    ///< optional sample identifier for logs
};

class TranscriptKindClassifier {
public:
    TranscriptKindClassifier() = default;

    /// Classify an anchor. May refine an existing kind set by the loader,
    /// or assign one when the loader left it Unknown.
    ///
    /// Returns the refined kind PLUS an optional CustomKindTag for
    /// open-ended cases (kind == NovelUnclassified).
    [[nodiscard]] std::pair<core::TranscriptKind,
                            std::optional<core::CustomKindTag>>
    Classify(const anchor::AnchorRecord& anchor,
              const SpliceSiteDb& splice,
              const TissueContext& tissue) const;

    // ----- Specialised detectors (public for direct use + tests) ---------

    /// True if the anchor looks like a sterile germline transcript. We
    /// require:
    ///   - tags contain "switch_region" OR host_gene_id starts with "IGH"/"IGK"/...
    ///   - no V/D/J segment tag present
    ///   - kind is currently Unknown or MatureMrna
    [[nodiscard]] bool LooksLikeSterileGermlineTranscript(
        const anchor::AnchorRecord& anchor) const;

    /// True if any ExonBoundary has spliceosome_class == 3 (back-splice).
    [[nodiscard]] bool LooksLikeCircularRna(
        const anchor::AnchorRecord& anchor,
        const SpliceSiteDb& splice) const;

    /// True if the anchor carries retained intron tags or is annotated
    /// as a pre-mRNA biotype.
    [[nodiscard]] bool LooksLikePreMrna(
        const anchor::AnchorRecord& anchor) const;

    /// piRNA detection from sequence: length 26-31 nt, no exon boundaries,
    /// often germline / somatic. We don't reproduce piRBase here; this is
    /// a fast structural check.
    [[nodiscard]] bool LooksLikePirna(
        const anchor::AnchorRecord& anchor) const;

    /// siRNA detection: length 21-23 nt, no exon boundaries, often
    /// generated from dsRNA precursors.
    [[nodiscard]] bool LooksLikeSirna(
        const anchor::AnchorRecord& anchor) const;

    /// Heuristic snoRNA subclass refinement:
    ///   C/D box: hallmark RUGAUGA … CUGA, length 60-120 nt
    ///   H/ACA box: hallmark ANANNA … ACA, length 120-180 nt
    /// Returns the refined kind, or Unknown if neither matches.
    [[nodiscard]] core::TranscriptKind RefineSnornaSubclass(
        const anchor::AnchorRecord& anchor) const;

    /// Open-ended fallback. Returns NovelUnclassified plus a label
    /// describing the heuristic signature that triggered.
    [[nodiscard]] std::pair<core::TranscriptKind, core::CustomKindTag>
    ClassifyOpenEnded(const anchor::AnchorRecord& anchor) const;
};

}  // namespace llmap::annot
