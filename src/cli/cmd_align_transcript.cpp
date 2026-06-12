// LLmap — Transcript-Mode record building implementation.

#include "cli/cmd_align_transcript.h"

#include <algorithm>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::cli::align_internal {

namespace {

// Junction probability from the splice-site PWM scores (each already in
// [0,1]). Both ends must look like a real splice site, so we take the
// product, with a lossless floor so a sub-chain is never silently dropped.
float JunctionProbFromSplice(const annot::SpliceScoreResult& s) {
    const float p = s.donor_score * s.acceptor_score;
    return std::max(0.05f, std::min(1.0f, p));
}

AlignmentHit SplicedToHit(const mapping::SplicedAlignment& sa) {
    AlignmentHit hit;
    hit.target_id = sa.ref_id;
    hit.start = sa.ref_start;
    hit.end = sa.ref_end;
    hit.cigar = CigarString{sa.cigar};
    hit.score = sa.score;
    hit.nm = 0;  // recomputed downstream if needed; spliced gap is not a mismatch
    hit.is_reverse = (sa.strand == '-');
    return hit;
}

AlignmentHit ClassicalToHit(const classical::ClassicalAlignment& aln) {
    AlignmentHit hit;
    hit.target_id = aln.ref_name;
    hit.start = static_cast<std::uint64_t>(aln.ref_start);
    hit.end = static_cast<std::uint64_t>(aln.ref_end);
    hit.cigar = CigarString{aln.CigarString()};
    hit.score = aln.score;
    hit.nm =
        static_cast<std::uint32_t>((1.0f - aln.identity) * aln.AlignedBases());
    hit.is_reverse = !aln.is_forward;
    return hit;
}

// Ranking key: prefer spliced (multi-exon) alignments, then higher score,
// then longer query span.
bool Better(const mapping::SplicedAlignment& a,
            const mapping::SplicedAlignment& b) {
    if (a.is_spliced != b.is_spliced) return a.is_spliced;
    if (a.score != b.score) return a.score > b.score;
    return (a.query_end - a.query_start) > (b.query_end - b.query_start);
}

}  // namespace

mapping::RefSeqLookup MakeRefLookup(const std::vector<std::string>& names,
                                    const std::vector<std::string>& seqs) {
    auto index = std::make_shared<std::unordered_map<std::string, std::string_view>>();
    const std::size_t n = std::min(names.size(), seqs.size());
    index->reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        (*index)[names[i]] = std::string_view(seqs[i]);
    }
    return [index](std::string_view id) -> std::string_view {
        const auto it = index->find(std::string(id));
        return it == index->end() ? std::string_view{} : it->second;
    };
}

AlignmentRecord BuildTranscriptRecord(const classical::ReadAlignmentResult& res,
                                      std::uint32_t read_len,
                                      const mapping::RefSeqLookup& ref_lookup,
                                      const annot::SpliceSiteDb& splice_db) {
    // All per-exon alignments become linear sub-chains.
    std::vector<mapping::LinearSubChain> subs;
    subs.reserve(res.alignments.size());
    for (const auto& aln : res.alignments) {
        mapping::LinearSubChain s;
        s.ref_id = aln.ref_name;
        s.ref_start = static_cast<std::uint64_t>(aln.ref_start);
        s.ref_end = static_cast<std::uint64_t>(aln.ref_end);
        s.query_start = static_cast<std::uint32_t>(aln.query_start);
        s.query_end = static_cast<std::uint32_t>(aln.query_end);
        s.score = aln.score;
        s.strand = aln.is_forward ? '+' : '-';
        s.cigar = aln.CigarString();
        subs.push_back(std::move(s));
    }

    mapping::JunctionScorer scorer =
        [&splice_db](std::string_view d, std::string_view a, std::string_view i3,
                     std::string_view i5) -> float {
        return JunctionProbFromSplice(splice_db.ScoreJunction(d, a, i3, i5));
    };

    // R-A: in Transcript-Mode the joiner threshold is the lossless floor, not a
    // canonicality gate. An intron-sized, query-colinear gap (the geometry
    // checks in GapLooksLikeIntron) is merged even when the splice motif is
    // weak/non-canonical — exactly the novel-alt-exon / NMD-escape case the
    // reference annotation would miss. The motif strength is preserved as the
    // jM junction confidence; it never blocks the merge or fragments the read.
    mapping::TranscriptStageConfig cfg;
    cfg.joiner.min_junction_probability = 0.05f;  // == the scorer's floor

    auto spliced = mapping::ApplyTranscriptStage(subs, ref_lookup, scorer, cfg);

    if (spliced.empty()) {
        // Should not happen (res.HasAlignment()), but stay safe: primary only.
        const auto* primary = res.PrimaryAlignment();
        return make_mapped(res.query_name, read_len, ClassicalToHit(*primary));
    }

    std::size_t best = 0;
    for (std::size_t i = 1; i < spliced.size(); ++i) {
        if (Better(spliced[i], spliced[best])) best = i;
    }

    AlignmentHit primary = SplicedToHit(spliced[best]);
    std::vector<AlignmentHit> alts;
    alts.reserve(spliced.size() - 1);
    for (std::size_t i = 0; i < spliced.size(); ++i) {
        if (i != best) alts.push_back(SplicedToHit(spliced[i]));
    }
    return make_mapped(res.query_name, read_len, std::move(primary),
                       std::move(alts));
}

}  // namespace llmap::cli::align_internal
