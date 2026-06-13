// LLmap — consensus-contig placement propagation (assemble-then-map, line A).
//
// After the transient probe of a ConsensusContig is mapped to a locus, the
// placement propagates to every anchored member read — LOSSLESSLY: each read keeps
// its OWN sequence; only the per-read reference coordinate + effective strand are
// derived from the contig-frame layout. The probe's bases are never propagated
// (the probe is a throwaway region-finder; the truth is the reads).
//
// Coordinate transform (contig frame → reference):
//   * probe mapped FORWARD: contig column c → ref (probe.ref_start + c).
//     member at offset o → ref_start = probe.ref_start + o;
//     effective read strand = m.reverse.
//   * probe mapped REVERSE: the contig is reverse-complemented onto the ref, so
//     contig column c → ref (probe.ref_start + length - 1 - c). A member spanning
//     contig columns [o, o+L) lands at ref_start = probe.ref_start + length - (o+L);
//     effective read strand = !m.reverse.
//   where L = the member read's own length (reads[read_idx].size()).
//
// Unanchored members (anchored == false) are NOT placed here — the caller maps them
// individually (lossless: no read is dropped, it just falls back to solo mapping).

#pragma once

#include "mapper/cluster_consensus.h"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace llmap::mapper {

// Where the transient probe landed (probe[0] is at ref_start in the probe-forward
// frame; reverse == the probe mapped to the reverse strand).
struct ProbePlacement {
    std::string  ref_name;
    std::int64_t ref_start{0};   // 0-based reference coordinate of probe column 0
    bool         reverse{false};
};

// A propagated per-read placement. The read's own sequence (reads[read_idx],
// reverse-complemented iff `reverse`) is what gets aligned/CIGAR'd against the
// region — never the probe.
struct MemberPlacement {
    std::size_t  read_idx{0};
    std::string  ref_name;
    std::int64_t ref_start{0};   // 0-based ref coordinate of the read's 5'-on-ref base
    bool         reverse{false}; // effective read strand on the reference
    bool         anchored{false};// false ⇒ not placed by the contig; map this read solo
};

// Propagate a probe placement to all members. Unanchored members are returned with
// anchored == false (and ref_start 0) so the caller can route them to solo mapping.
[[nodiscard]] std::vector<MemberPlacement> PropagatePlacement(
    const ConsensusContig& contig,
    std::span<const std::string> reads,
    const ProbePlacement& probe);

}  // namespace llmap::mapper
