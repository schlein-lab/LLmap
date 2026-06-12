// LLmap — Transcript-Mode record building for `llmap align`.
//
// When `--mode transcript` is active, a spliced read is reported by the
// classical chainer as several per-exon alignments. This helper joins them
// (mapping::ApplyTranscriptStage) into one spliced AlignmentRecord with an
// N-op CIGAR. The DNA/GenomeReads path does not call any of this.

#pragma once

#include <string>
#include <vector>

#include "annot/splice_site_db.h"
#include "classical/classical_pipeline.h"
#include "core/alignment_record.h"
#include "mapping/transcript_stage.h"

namespace llmap::cli::align_internal {

// Build a ref_id → sequence lookup from the parallel reference vectors loaded
// in cmd_align. The returned functor holds a private index; the string_views
// point into `seqs`, which must outlive the functor (it does — both live for
// the whole align run).
[[nodiscard]] mapping::RefSeqLookup MakeRefLookup(
    const std::vector<std::string>& names,
    const std::vector<std::string>& seqs);

// Turn one read's classical alignments into a spliced AlignmentRecord. Falls
// back to the primary alignment when nothing can be joined. Caller guarantees
// res.HasAlignment().
[[nodiscard]] AlignmentRecord BuildTranscriptRecord(
    const classical::ReadAlignmentResult& res,
    std::uint32_t read_len,
    const mapping::RefSeqLookup& ref_lookup,
    const annot::SpliceSiteDb& splice_db);

}  // namespace llmap::cli::align_internal
