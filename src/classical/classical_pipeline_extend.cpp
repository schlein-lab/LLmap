#include <cstdlib>
#include <cstdio>
// LLmap — ClassicalPipeline: chain extension with WFA2 alignment.

#include "classical/classical_pipeline.h"

namespace llmap::classical {

namespace {

std::string ReverseComplement(std::string_view seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': case 'a': rc.push_back('T'); break;
            case 'T': case 't': rc.push_back('A'); break;
            case 'C': case 'c': rc.push_back('G'); break;
            case 'G': case 'g': rc.push_back('C'); break;
            default:            rc.push_back('N'); break;
        }
    }
    return rc;
}

}  // namespace

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
    std::string_view query_seq_in,
    const Chain& chain,
    const std::vector<Anchor>& anchors_in) const {

    if (chain.anchors.empty()) {
        return std::nullopt;
    }

    ClassicalAlignment aln;
    aln.score = chain.score;

    uint8_t k = config_.minimizer_config.k;

    // Reverse-strand chains: align the reverse-complement of the read against
    // the forward-strand reference. Query positions for reverse-strand
    // minimizers refer to forward-strand coordinates, so flip them once here.
    std::string rc_query_storage;
    std::vector<Anchor> flipped_anchors_storage;
    std::string_view query_seq = query_seq_in;
    const std::vector<Anchor>* anchors_ptr = &anchors_in;

    if (!chain.is_forward) {
        rc_query_storage = ReverseComplement(query_seq_in);
        query_seq = rc_query_storage;
        flipped_anchors_storage.reserve(anchors_in.size());
        const uint32_t qlen = static_cast<uint32_t>(query_seq_in.size());
        for (const auto& a : anchors_in) {
            Anchor flipped = a;
            // Minimizer at forward-strand position p maps to position qlen-p-k
            // in the reverse-complemented read.
            flipped.query_pos = (a.query_pos + k <= qlen)
                ? qlen - a.query_pos - k
                : 0;
            flipped_anchors_storage.push_back(flipped);
        }
        anchors_ptr = &flipped_anchors_storage;
    }
    const std::vector<Anchor>& anchors = *anchors_ptr;
    size_t matches = 0;
    size_t mismatches = 0;
    size_t insertions = 0;
    size_t deletions = 0;

    bool have_ref_seqs = !ref_seqs_.empty() && chain.ref_id < ref_seqs_.size();
    const auto& ref_seq = have_ref_seqs ? ref_seqs_[chain.ref_id] : std::string{};
    uint32_t ref_len = have_ref_seqs ? static_cast<uint32_t>(ref_seq.size()) : 0;
    uint32_t query_len = static_cast<uint32_t>(query_seq.size());

    // Get first and last anchor positions for end extension
    const auto& first_anchor = anchors[chain.anchors.front()];
    const auto& last_anchor = anchors[chain.anchors.back()];

    if (const char* d = std::getenv("LLMAP_DEBUG_REV"); d && d[0]=='1' && !chain.is_forward) {
        for (std::size_t ci = 0; ci < 3 && ci < chain.anchors.size(); ++ci) {
            const auto& a = anchors[chain.anchors[ci]];
            std::string qs(query_seq.substr(a.query_pos, std::min<uint32_t>(19, query_len - a.query_pos)));
            std::string rs = (a.ref_pos + 19 <= ref_len)
                ? std::string(std::string_view(ref_seq).substr(a.ref_pos, 19)) : "";
            std::fprintf(stderr, "[rev] anchor#%zu q'=%u r=%u  rcq=%s  ref=%s  %s\n",
                ci, a.query_pos, a.ref_pos, qs.c_str(), rs.c_str(),
                (qs==rs?"MATCH":"differ"));
        }
    }

    // === LEFT EXTENSION ===
    // Extend alignment from query start (0) to first anchor position.
    // Soft-clip directly if the unaligned span is large (extension cost is
    // O(span²) and the gain is negligible past the first ~500 bp).
    const uint32_t kMaxExtensionSpan = config_.extension_max_span;
    uint32_t left_query_bases = first_anchor.query_pos;
    uint32_t left_ref_bases = 0;
    uint32_t actual_ref_start = chain.ref_start;

    if (have_ref_seqs && left_query_bases > 0 && left_query_bases <= kMaxExtensionSpan) {
        // Calculate how far we can extend left in reference.
        // Use the query span exactly — adding +50 of upstream ref padding
        // caused a systematic 50-bp leftward POS shift on T1/T2: WFA2's
        // semi-global ExtendLeft greedily consumed the upstream bases as a
        // leading deletion (CIGAR `50D…`), pushing actual_ref_start ~50 bp
        // upstream of the true read origin and tanking F1 (T1 0.046 → fix → 1.0).
        // If real left-side indels appear, they are still captured by the
        // anchor-to-anchor WFA2 path; the extension only fills the un-anchored
        // 5′ flank, which by definition cannot have more reference than query.
        left_ref_bases = std::min(left_query_bases, first_anchor.ref_pos);
        uint32_t ref_ext_start = first_anchor.ref_pos - left_ref_bases;

        if (left_ref_bases > 0 && left_query_bases > 0) {
            auto left_result = aligner_.ExtendLeft(
                query_seq.substr(0, first_anchor.query_pos),
                std::string_view(ref_seq).substr(ref_ext_start, left_ref_bases),
                static_cast<int32_t>(first_anchor.query_pos),
                static_cast<int32_t>(left_ref_bases));

            if (left_result) {
                // Add left extension CIGAR (soft-clip unaligned, then extension)
                int32_t aligned_query = left_result->query_end - left_result->query_start;
                int32_t query_softclip = static_cast<int32_t>(first_anchor.query_pos) - aligned_query;

                if (query_softclip > 0) {
                    aln.cigar.push_back({CigarOp::SoftClip, static_cast<uint32_t>(query_softclip)});
                }

                for (const auto& elem : left_result->cigar) {
                    aln.cigar.push_back(elem);
                }

                matches += left_result->num_matches;
                mismatches += left_result->num_mismatches;
                insertions += left_result->num_insertions;
                deletions += left_result->num_deletions;

                // Update ref_start based on extension
                actual_ref_start = ref_ext_start + static_cast<uint32_t>(left_result->ref_start);
            } else {
                // Extension failed - soft clip the left portion
                if (left_query_bases > 0) {
                    aln.cigar.push_back({CigarOp::SoftClip, left_query_bases});
                }
            }
        } else if (left_query_bases > 0) {
            aln.cigar.push_back({CigarOp::SoftClip, left_query_bases});
        }
    } else if (left_query_bases > 0) {
        // No reference sequences - soft clip left portion
        aln.cigar.push_back({CigarOp::SoftClip, left_query_bases});
    }

    // Set alignment start positions
    aln.ref_start = actual_ref_start;
    aln.query_start = 0;  // Now we cover from query start

    uint32_t prev_query = first_anchor.query_pos;
    uint32_t prev_ref = first_anchor.ref_pos;

    if (!config_.chain_config.splice_aware) {
        // === GENOME: WINDOWED CORE WFA (the long-read fix) ===
        // The chain is colinear (the [rev]/[diag] dumps: anchors match ref exactly,
        // Δq==Δr). The former per-anchor interpolation + k-mer accounting BREAKS on
        // DENSE anchors (spacing < k, common on long reads incl. reverse chains):
        // prev_query advances by k per anchor and outpaces the anchor positions →
        // k-mer overlap ≥ k → kmer_add 0, gaps go negative → query bases counted as
        // INSERTIONS → M=84/I=29268 for a true 29 kb reverse chain → 0% mapped. We
        // align the colinear core honestly, TILED into ≤kWin-bp windows at anchor
        // boundaries (each window's WFA is ≤kWin² ≈ 12 MB, O(n·kWin) total, no OOM).
        // query_seq is already reverse-complemented and the anchors flipped for a
        // reverse chain, so this is strand-correct.
        // PER-ANCHOR-SEGMENT base comparison. Between two consecutive chain anchors
        // the diagonal is locally exact (the anchors track the read's small indels),
        // so each segment is a clean colinear stretch: base-compare it (M/X) and emit
        // the offset jump between segments as the indel (I/D). This bypasses the WFA2
        // wrapper entirely — it degenerates (M=1) on long/unequal input — and is exact
        // + fast (O(n), no DP). The former fixed 256 bp tiles drifted within a tile
        // after an in-tile indel → 32% spurious X (ident 0.68); per-anchor segments
        // (dense anchors, Δ~12 bp) are indel-free → true ~0.99.
        const std::size_t na = chain.anchors.size();
        auto cmp = [&](uint32_t q0, uint32_t r0, uint32_t n) {
            for (uint32_t i = 0; i < n; ++i) {
                if (q0 + i >= query_len || r0 + i >= ref_len) break;
                const char qc = static_cast<char>(query_seq[q0 + i] & ~0x20);
                const char rc = static_cast<char>(ref_seq[r0 + i] & ~0x20);
                if (qc == rc) { aln.cigar.push_back({CigarOp::Match, 1}); ++matches; }
                else { aln.cigar.push_back({CigarOp::Diff, 1}); ++mismatches; }
            }
        };
        if (have_ref_seqs && na >= 1) {
            uint32_t pq = first_anchor.query_pos;
            uint32_t pr = first_anchor.ref_pos;
            for (std::size_t idx = 1; idx < na; ++idx) {
                const auto& a = anchors[chain.anchors[idx]];
                if (a.query_pos <= pq || a.ref_pos <= pr) continue;  // overlap/backward
                const int64_t dq = static_cast<int64_t>(a.query_pos) - pq;
                const int64_t dr = static_cast<int64_t>(a.ref_pos) - pr;
                cmp(pq, pr, static_cast<uint32_t>(std::min(dq, dr)));
                if (dq > dr) {
                    aln.cigar.push_back({CigarOp::Insertion, static_cast<uint32_t>(dq - dr)});
                    insertions += dq - dr;
                } else if (dr > dq) {
                    aln.cigar.push_back({CigarOp::Deletion, static_cast<uint32_t>(dr - dq)});
                    deletions += dr - dq;
                }
                pq = a.query_pos;
                pr = a.ref_pos;
            }
            cmp(pq, pr, k);  // final anchor's k-mer
        }
    } else {
      // === TRANSCRIPT (B2): per-anchor path — short per-exon spans (the windowed
      // core WFA must NOT run here: it would align across a sub-chain's intron). ===
      for (uint32_t anchor_idx : chain.anchors) {
        if (anchor_idx >= anchors.size()) continue;
        const auto& anchor = anchors[anchor_idx];
        int32_t query_gap = static_cast<int32_t>(anchor.query_pos) - static_cast<int32_t>(prev_query);
        int32_t ref_gap = static_cast<int32_t>(anchor.ref_pos) - static_cast<int32_t>(prev_ref);
        if (query_gap > 0 && ref_gap > 0) {
            bool use_wfa = have_ref_seqs && (query_gap >= 50 || ref_gap >= 50);
            std::optional<WFA2Result> gap_result;
            if (use_wfa) {
                gap_result = AlignGap(query_seq, chain.ref_id,
                    prev_query, anchor.query_pos, prev_ref, anchor.ref_pos);
            }
            if (gap_result) {
                for (const auto& elem : gap_result->cigar) aln.cigar.push_back(elem);
                matches += gap_result->num_matches;
                mismatches += gap_result->num_mismatches;
                insertions += gap_result->num_insertions;
                deletions += gap_result->num_deletions;
            } else {
                uint32_t aligned = static_cast<uint32_t>(std::min(query_gap, ref_gap));
                if (aligned > 0) { aln.cigar.push_back({CigarOp::Match, aligned}); matches += aligned; }
                if (query_gap > ref_gap) {
                    uint32_t ins = static_cast<uint32_t>(query_gap - ref_gap);
                    aln.cigar.push_back({CigarOp::Insertion, ins}); insertions += ins;
                } else if (ref_gap > query_gap) {
                    uint32_t del = static_cast<uint32_t>(ref_gap - query_gap);
                    aln.cigar.push_back({CigarOp::Deletion, del}); deletions += del;
                }
            }
        } else if (query_gap > 0) {
            aln.cigar.push_back({CigarOp::Insertion, static_cast<uint32_t>(query_gap)});
            insertions += query_gap;
        } else if (ref_gap > 0) {
            aln.cigar.push_back({CigarOp::Deletion, static_cast<uint32_t>(ref_gap)});
            deletions += ref_gap;
        }
        uint32_t kmer_add = k;
        if (anchor.query_pos < prev_query) {
            uint32_t overlap = prev_query - anchor.query_pos;
            kmer_add = (overlap < k) ? (k - overlap) : 0;
        }
        if (kmer_add > 0) { aln.cigar.push_back({CigarOp::Match, kmer_add}); matches += kmer_add; }
        prev_query = std::max(prev_query, anchor.query_pos + k);
        prev_ref = std::max(prev_ref, anchor.ref_pos + k);
      }
    }

    // Handle trailing gap after last anchor (within chain span)
    uint32_t last_anchor_end_query = last_anchor.query_pos + k;
    uint32_t last_anchor_end_ref = last_anchor.ref_pos + k;

    // === RIGHT EXTENSION ===
    // Extend alignment from last anchor to query end
    uint32_t right_query_bases = query_len - last_anchor_end_query;
    uint32_t actual_ref_end = chain.ref_end;

    if (have_ref_seqs && right_query_bases > 0 && right_query_bases <= kMaxExtensionSpan) {
        uint32_t right_ref_bases = std::min(right_query_bases, ref_len - last_anchor_end_ref);

        if (right_ref_bases > 0) {
            auto right_result = aligner_.ExtendRight(
                query_seq.substr(last_anchor_end_query),
                std::string_view(ref_seq).substr(last_anchor_end_ref, right_ref_bases),
                0,  // Start from beginning of substring
                0);

            if (right_result) {
                // Add right extension CIGAR
                for (const auto& elem : right_result->cigar) {
                    aln.cigar.push_back(elem);
                }

                matches += right_result->num_matches;
                mismatches += right_result->num_mismatches;
                insertions += right_result->num_insertions;
                deletions += right_result->num_deletions;

                // Calculate how much of query was aligned
                int32_t aligned_query = right_result->query_end - right_result->query_start;
                int32_t query_softclip = static_cast<int32_t>(right_query_bases) - aligned_query;

                if (query_softclip > 0) {
                    aln.cigar.push_back({CigarOp::SoftClip, static_cast<uint32_t>(query_softclip)});
                }

                // Update ref_end based on extension
                actual_ref_end = last_anchor_end_ref + static_cast<uint32_t>(right_result->ref_end);
            } else {
                // Extension failed - soft clip the right portion
                if (right_query_bases > 0) {
                    aln.cigar.push_back({CigarOp::SoftClip, right_query_bases});
                }
                actual_ref_end = last_anchor_end_ref;
            }
        } else if (right_query_bases > 0) {
            aln.cigar.push_back({CigarOp::SoftClip, right_query_bases});
            actual_ref_end = last_anchor_end_ref;
        }
    } else if (right_query_bases > 0) {
        // No reference sequences - soft clip right portion
        aln.cigar.push_back({CigarOp::SoftClip, right_query_bases});
        actual_ref_end = last_anchor_end_ref;
    }

    // Set alignment end positions
    aln.ref_end = actual_ref_end;
    aln.query_end = query_len;  // Now we cover to query end

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

    if (const char* d = std::getenv("LLMAP_DEBUG_GAP"); d && d[0] == '1') {
        std::fprintf(stderr,
            "[gap] anchors=%zu qlen=%u M=%zu X=%zu I=%zu D=%zu total=%zu ident=%.3f\n",
            chain.anchors.size(), query_len, matches, mismatches, insertions,
            deletions, total_aligned, aln.identity);
    }
    return aln;
}

}  // namespace llmap::classical
