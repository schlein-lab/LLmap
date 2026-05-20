// LLmap — post-hoc IGH locus re-sort.
//
// Runs AFTER standard mapping. For every read whose primary alignment lands in
// the IGH constant region (genomic) or on an IGH transcript/copy target
// (transcriptomic), it re-scans the read against paralog-specific exon anchors
// and, when the exact-anchor evidence contradicts the mapper's call, re-points
// the read to its true copy of origin and records the paralog assignment.
//
// The same entry point serves both the integrated pipeline (mutating in-memory
// AlignmentRecords before BAM output) and the standalone `igh-resort` command
// (re-sorting an existing minimap2/other BAM via the SAM stream).

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "core/alignment_record.h"
#include "igh_locus/igh_anchor_catalog.h"
#include "igh_locus/igh_region.h"

namespace llmap::igh_locus {

struct ResortOptions {
    // Relabel primary target to the anchor-implied copy on disagreement.
    bool relabel{true};
    // Require at least this many distinct exons of one copy before acting.
    int min_exons{2};
    // Populate record.paralog_assignment from the anchor evidence.
    bool set_paralog{true};
    bool verbose{false};
    // Genomic gating region (ignored for transcriptomic targets).
    IghRegion region{IghRegion::Default()};
};

struct ResortStats {
    std::size_t n_examined{0};   // records seen
    std::size_t n_in_locus{0};   // mapped IGH reads considered
    std::size_t n_matched{0};    // reads with a confident anchor call
    std::size_t n_resorted{0};   // reads whose target/coords were changed
    std::size_t n_paralog_set{0};// reads whose paralog_assignment was set
};

// Mutate `records` in place. `read_seqs` must be parallel to `records` (same
// index = same read); entries may be empty when a sequence is unavailable, in
// which case that read is left untouched.
ResortStats ApplyResort(const IghAnchorCatalog& catalog,
                        std::vector<AlignmentRecord>& records,
                        const std::vector<std::string>& read_seqs,
                        const ResortOptions& opts = {});

}  // namespace llmap::igh_locus
