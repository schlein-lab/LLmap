// LLmap — Transcript-Mode spliced stage implementation.

#include "mapping/transcript_stage.h"

#include <algorithm>
#include <cctype>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace llmap::mapping {

namespace {

char Upper(char c) noexcept {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

char Complement(char c) noexcept {
    switch (Upper(c)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

// Uppercase copy of a reference slice [start, end).
std::string Slice(std::string_view ref, std::uint64_t start, std::uint64_t end) {
    std::string s;
    if (start >= end || end > ref.size()) return s;
    s.reserve(end - start);
    for (std::uint64_t i = start; i < end; ++i) s.push_back(Upper(ref[i]));
    return s;
}

// Reverse-complement (uppercased).
std::string RevComp(std::string_view ref, std::uint64_t start, std::uint64_t end) {
    std::string s;
    if (start >= end || end > ref.size()) return s;
    s.reserve(end - start);
    for (std::uint64_t i = end; i > start; --i) s.push_back(Complement(ref[i - 1]));
    return s;
}

std::uint64_t Min64(std::uint64_t a, std::uint64_t b) { return a < b ? a : b; }
std::uint64_t Max64(std::uint64_t a, std::uint64_t b) { return a > b ? a : b; }

}  // namespace

// ---------------------------------------------------------------------------
// ExtractJunctionMotifs
// ---------------------------------------------------------------------------

bool ExtractJunctionMotifs(const LinearSubChain& a, const LinearSubChain& b,
                           std::string_view ref_seq, std::string& donor2,
                           std::string& acceptor2, std::string& intron3p,
                           std::string& intron5p) {
    // Intron occupies the reference gap [a.ref_end, b.ref_start) (0-based,
    // half-open). Need at least 2 bp at each end and in-bounds positions.
    const std::uint64_t i_lo = a.ref_end;       // first intron base
    const std::uint64_t i_hi = b.ref_start;     // one past last intron base
    if (i_hi < i_lo + 4) return false;          // need room for both 2-mers
    if (i_hi > ref_seq.size()) return false;

    if (a.strand != '-') {
        // Sense == reference forward. Donor at intron 5' end, acceptor at 3'.
        donor2    = Slice(ref_seq, i_lo, i_lo + 2);
        acceptor2 = Slice(ref_seq, i_hi - 2, i_hi);
        intron5p  = Slice(ref_seq, i_lo, Min64(i_lo + 8, i_hi));
        intron3p  = Slice(ref_seq, Max64(i_hi >= 50 ? i_hi - 50 : 0, i_lo), i_hi);
    } else {
        // Sense == reverse complement. The transcript's donor is at the high
        // reference coordinate (b side), acceptor at the low side; both 2-mers
        // are reverse-complemented so the scorer sees sense GT/AG-class motifs.
        donor2    = RevComp(ref_seq, i_hi - 2, i_hi);
        acceptor2 = RevComp(ref_seq, i_lo, i_lo + 2);
        intron5p  = RevComp(ref_seq, Max64(i_hi >= 8 ? i_hi - 8 : 0, i_lo), i_hi);
        intron3p  = RevComp(ref_seq, i_lo, Min64(i_lo + 50, i_hi));
    }
    return !donor2.empty() && !acceptor2.empty();
}

// ---------------------------------------------------------------------------
// ApplyTranscriptStage
// ---------------------------------------------------------------------------

std::vector<SplicedAlignment> ApplyTranscriptStage(
    std::span<const LinearSubChain> sub_chains, const RefSeqLookup& ref_lookup,
    const JunctionScorer& score_junction, const TranscriptStageConfig& cfg) {

    std::vector<SplicedAlignment> out;
    if (sub_chains.empty()) return out;

    // Group by (ref_id, strand), preserving first-seen order for determinism.
    using Key = std::pair<std::string, char>;
    std::vector<std::pair<Key, std::vector<LinearSubChain>>> groups;
    auto group_for = [&](const std::string& ref, char strand)
        -> std::vector<LinearSubChain>& {
        for (auto& g : groups) {
            if (g.first.first == ref && g.first.second == strand) return g.second;
        }
        groups.push_back({Key{ref, strand}, {}});
        return groups.back().second;
    };
    for (const auto& sc : sub_chains) group_for(sc.ref_id, sc.strand).push_back(sc);

    for (auto& [key, group] : groups) {
        // Order along the reference (the order JoinSplicedChains expects).
        std::sort(group.begin(), group.end(),
                  [](const LinearSubChain& x, const LinearSubChain& y) {
                      if (x.ref_start != y.ref_start) return x.ref_start < y.ref_start;
                      return x.query_start < y.query_start;
                  });

        // Per-pair junction probabilities from reference splice motifs.
        std::vector<float> probs;
        if (group.size() > 1) {
            probs.reserve(group.size() - 1);
            const std::string_view ref =
                ref_lookup ? ref_lookup(key.first) : std::string_view{};
            for (std::size_t i = 0; i + 1 < group.size(); ++i) {
                float p = 0.5f;  // default when motifs / scorer unavailable
                std::string d, acc, i3, i5;
                if (!ref.empty() && score_junction &&
                    ExtractJunctionMotifs(group[i], group[i + 1], ref, d, acc, i3,
                                          i5)) {
                    p = score_junction(d, acc, i3, i5);
                    p = std::clamp(p, 0.0f, 1.0f);
                }
                probs.push_back(p);
            }
        }

        const SplicedChainResult res = JoinSplicedChains(group, probs, cfg.joiner);

        for (const auto& chain : res.chains) {
            const auto& subs = chain.sub_chains;
            if (subs.empty()) continue;
            SplicedAlignment a;
            a.ref_id = chain.ref_id;
            a.strand = chain.strand;
            a.score = chain.total_score;
            a.ref_start = subs.front().ref_start;
            a.ref_end = subs.back().ref_end;
            std::uint32_t qs = subs.front().query_start;
            std::uint32_t qe = subs.front().query_end;
            for (const auto& s : subs) {
                qs = std::min(qs, s.query_start);
                qe = std::max(qe, s.query_end);
            }
            a.query_start = qs;
            a.query_end = qe;
            a.cigar = EmitSplicedCigar(chain);
            a.is_spliced = subs.size() > 1;
            for (const auto& j : chain.junctions) {
                a.junctions.emplace_back(j.donor_ref_pos, j.acceptor_ref_pos);
                a.junction_conf.push_back(j.probability);
            }
            out.push_back(std::move(a));
        }
    }
    return out;
}

}  // namespace llmap::mapping
