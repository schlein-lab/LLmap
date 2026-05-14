// LLmap — ClassicalPipeline: chain extension with WFA2 alignment.

#include "classical/classical_pipeline.h"

namespace llmap::classical {

std::optional<WFA2Result> ClassicalPipeline::AlignGap(
    std::string_view query_seq,
    uint32_t ref_id,
    uint32_t query_start,
    uint32_t query_end,
    uint32_t ref_start,
    uint32_t ref_end) const {

    if (ref_id >= ref_seqs_.size()) {
        return std::nullopt;
    }

    const auto& ref_seq = ref_seqs_[ref_id];
    if (ref_end > ref_seq.size()) {
        ref_end = static_cast<uint32_t>(ref_seq.size());
    }
    if (ref_start >= ref_end || query_start >= query_end) {
        return std::nullopt;
    }

    auto query_substr = query_seq.substr(query_start, query_end - query_start);
    auto ref_substr = std::string_view(ref_seq).substr(ref_start, ref_end - ref_start);

    return aligner_.Align(query_substr, ref_substr);
}

std::optional<ClassicalAlignment> ClassicalPipeline::ExtendChain(
    std::string_view query_seq,
    const Chain& chain,
    const std::vector<Anchor>& anchors) const {

    if (chain.anchors.empty()) {
        return std::nullopt;
    }

    ClassicalAlignment aln;
    aln.ref_start = chain.ref_start;
    aln.ref_end = chain.ref_end;
    aln.query_start = chain.query_start;
    aln.query_end = chain.query_end;
    aln.score = chain.score;

    uint8_t k = config_.minimizer_config.k;
    uint32_t prev_query = chain.query_start;
    uint32_t prev_ref = chain.ref_start;
    size_t matches = 0;
    size_t mismatches = 0;
    size_t insertions = 0;
    size_t deletions = 0;

    bool have_ref_seqs = !ref_seqs_.empty() && chain.ref_id < ref_seqs_.size();

    for (uint32_t anchor_idx : chain.anchors) {
        if (anchor_idx >= anchors.size()) continue;
        const auto& anchor = anchors[anchor_idx];

        int32_t query_gap = static_cast<int32_t>(anchor.query_pos) - static_cast<int32_t>(prev_query);
        int32_t ref_gap = static_cast<int32_t>(anchor.ref_pos) - static_cast<int32_t>(prev_ref);

        if (query_gap > 0 && ref_gap > 0) {
            // Try WFA2 alignment for the gap between anchors
            if (have_ref_seqs && query_gap >= 1 && ref_gap >= 1) {
                auto gap_result = AlignGap(
                    query_seq, chain.ref_id,
                    prev_query, anchor.query_pos,
                    prev_ref, anchor.ref_pos);

                if (gap_result) {
                    // Use WFA2 CIGAR
                    for (const auto& elem : gap_result->cigar) {
                        aln.cigar.push_back(elem);
                    }
                    matches += gap_result->num_matches;
                    mismatches += gap_result->num_mismatches;
                    insertions += gap_result->num_insertions;
                    deletions += gap_result->num_deletions;
                } else {
                    // Fallback to interpolation
                    uint32_t aligned = static_cast<uint32_t>(std::min(query_gap, ref_gap));
                    if (aligned > 0) {
                        aln.cigar.push_back({CigarOp::Match, aligned});
                        matches += aligned;
                    }
                    if (query_gap > ref_gap) {
                        uint32_t ins = static_cast<uint32_t>(query_gap - ref_gap);
                        aln.cigar.push_back({CigarOp::Insertion, ins});
                        insertions += ins;
                    } else if (ref_gap > query_gap) {
                        uint32_t del = static_cast<uint32_t>(ref_gap - query_gap);
                        aln.cigar.push_back({CigarOp::Deletion, del});
                        deletions += del;
                    }
                }
            } else {
                // Original interpolation fallback
                uint32_t aligned = static_cast<uint32_t>(std::min(query_gap, ref_gap));
                if (aligned > 0) {
                    aln.cigar.push_back({CigarOp::Match, aligned});
                    matches += aligned;
                }
                if (query_gap > ref_gap) {
                    uint32_t ins = static_cast<uint32_t>(query_gap - ref_gap);
                    aln.cigar.push_back({CigarOp::Insertion, ins});
                    insertions += ins;
                } else if (ref_gap > query_gap) {
                    uint32_t del = static_cast<uint32_t>(ref_gap - query_gap);
                    aln.cigar.push_back({CigarOp::Deletion, del});
                    deletions += del;
                }
            }
        } else if (query_gap > 0) {
            aln.cigar.push_back({CigarOp::Insertion, static_cast<uint32_t>(query_gap)});
            insertions += query_gap;
        } else if (ref_gap > 0) {
            aln.cigar.push_back({CigarOp::Deletion, static_cast<uint32_t>(ref_gap)});
            deletions += ref_gap;
        }

        // k-mer match (the anchor itself)
        aln.cigar.push_back({CigarOp::Match, k});
        matches += k;

        prev_query = anchor.query_pos + k;
        prev_ref = anchor.ref_pos + k;
    }

    // Handle trailing gap after last anchor
    if (prev_query < chain.query_end) {
        uint32_t trail_q = chain.query_end - prev_query;
        uint32_t trail_r = (chain.ref_end > prev_ref) ? chain.ref_end - prev_ref : 0;

        if (have_ref_seqs && trail_q > 0 && trail_r > 0) {
            auto gap_result = AlignGap(
                query_seq, chain.ref_id,
                prev_query, chain.query_end,
                prev_ref, chain.ref_end);

            if (gap_result) {
                for (const auto& elem : gap_result->cigar) {
                    aln.cigar.push_back(elem);
                }
                matches += gap_result->num_matches;
                mismatches += gap_result->num_mismatches;
                insertions += gap_result->num_insertions;
                deletions += gap_result->num_deletions;
            } else {
                aln.cigar.push_back({CigarOp::Match, trail_q});
                matches += trail_q;
            }
        } else if (trail_q > 0) {
            aln.cigar.push_back({CigarOp::Match, trail_q});
            matches += trail_q;
        }
    }

    // Merge adjacent CIGAR operations
    std::vector<CigarElement> merged;
    for (const auto& elem : aln.cigar) {
        if (elem.length == 0) continue;
        if (!merged.empty() && merged.back().op == elem.op) {
            merged.back().length += elem.length;
        } else {
            merged.push_back(elem);
        }
    }
    aln.cigar = std::move(merged);

    // Statistics
    size_t total_aligned = matches + mismatches + insertions + deletions;
    aln.identity = total_aligned > 0
        ? static_cast<float>(matches) / static_cast<float>(total_aligned)
        : 0.0f;

    return aln;
}

}  // namespace llmap::classical
