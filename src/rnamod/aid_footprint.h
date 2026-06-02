// LLmap — AID-mediated C→U editing footprint detector.
//
// AID (Activation-Induced Deaminase, AICDA) deaminates cytosine to
// uracil inside IGH switch regions during class-switch recombination
// (CSR) and somatic hypermutation (SHM). The C→U edits cluster in
// short repetitive Sγ/Sα/Sε/Sμ tracts and are detectable from long-read
// data as a C→T mismatch density that exceeds the platform baseline.
//
// Why a dedicated module:
//
//   - The signal is small (typically 5-20 edits across a few kb), so
//     calling it from generic read-vs-reference mismatch counting is
//     drowned by sequencing noise. We constrain detection to anchors
//     tagged as Switch_region (i.e. Sγ4 / Sγ3 / Sμ / Sα / Sε from
//     IMGT germline) and look for the WRC AID hotspot motif before
//     accepting an edit as bona fide AID.
//
//   - The footprint is THE marker for active class-switching, the
//     biology that drives the IGHG4-canonical-vs-dup story in
//     [[ighg4_sgamma4_identical_in_tandem_dup]]. Lossless detection
//     here = directly observable IGH class-switch dynamics.
//
//   - AID footprint is *orthogonal* to AlignmentStatus, TranscriptKind,
//     and SplicingState — a successfully mapped IGHG4 mRNA read can
//     also carry an AID signature on its switch-region tail. So we
//     expose AidFootprint as an opt-in payload attached to any
//     AlignmentRecord via its existing aid_footprint field.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::rnamod {

/// Per-read AID footprint summary.
struct AidFootprintResult {
    bool detected{false};
    std::uint32_t n_c_to_u_events{0};
    std::vector<std::uint32_t> c_to_u_positions_in_read;
    std::string switch_region_id;        ///< e.g. "S_gamma4"
    float confidence{0.0f};               ///< [0,1]; AID-hotspot fraction
};

class AidFootprintDetector {
public:
    /// Detect AID footprint on a single read aligned to an anchor.
    ///
    ///   read_seq          — the read sequence (5'→3', forward strand)
    ///   anchor_seq        — the corresponding reference window (same length)
    ///   anchor_switch_id  — the switch region id from the anchor tags
    ///                        (empty ⇒ not a switch-region anchor, returns
    ///                         detected=false)
    ///
    /// Minimum-evidence threshold: at least 3 C→U edits, of which at
    /// least 50 % must fall in an AID hotspot motif (WRC). Below the
    /// threshold, detected stays false but n_c_to_u_events is still
    /// populated so callers can inspect borderline cases.
    [[nodiscard]] AidFootprintResult Detect(
        std::string_view read_seq,
        std::string_view anchor_seq,
        std::string_view anchor_switch_id) const;
};

}  // namespace llmap::rnamod
