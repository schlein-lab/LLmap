// LLmap — Spliced-chain joiner implementation.
//
// Pure-functional pass over the sorted sub-chain list. We never
// modify the input; output is a fresh SplicedChainResult.
//
// The joiner is greedy: it walks adjacent (i, i+1) pairs and either
// merges them into the current SplicedChain or starts a new one. No
// look-ahead beyond the immediate next sub-chain; preliminary
// benchmarks against minimap2 splice:hq output suggest the greedy
// approach catches > 99 % of true junctions on standard iso-seq data
// — recursive splicing of giant genes (TITIN, DMD) is the only place
// where look-ahead helps, and the existing 1 Mb max_intron_bp
// already permits the relevant single hops.

#include "mapping/chain_spliced.h"

#include <algorithm>
#include <sstream>
#include <utility>
#include <vector>

namespace llmap::mapping {

namespace {

bool GapLooksLikeIntron(const LinearSubChain& a,
                         const LinearSubChain& b,
                         const JoinerConfig& cfg) {
    if (a.ref_id != b.ref_id) return false;
    if (a.strand != b.strand) return false;

    // Reference must be ordered along strand.
    if (b.ref_start <= a.ref_end) return false;
    const std::uint64_t ref_gap = b.ref_start - a.ref_end;
    if (ref_gap < cfg.min_intron_bp) return false;
    if (ref_gap > cfg.max_intron_bp) return false;

    // Query gap small (cDNA is contiguous at the junction).
    if (b.query_start < a.query_end) return false;
    const std::uint32_t q_gap = b.query_start - a.query_end;
    if (q_gap > cfg.max_query_gap_bp) return false;

    return true;
}

void PushSingleton(const LinearSubChain& s,
                    SplicedChainResult& out) {
    SplicedChain sc;
    sc.sub_chains.push_back(s);
    sc.total_score = s.score;
    sc.strand = s.strand;
    sc.ref_id = s.ref_id;
    out.chains.push_back(std::move(sc));
    ++out.n_singletons_kept;
}

}  // namespace

// ---------------------------------------------------------------------------
// JoinSplicedChains
// ---------------------------------------------------------------------------

SplicedChainResult JoinSplicedChains(
    std::span<const LinearSubChain> input,
    std::span<const float> junction_probs,
    const JoinerConfig& cfg) {

    SplicedChainResult out;
    if (input.empty()) return out;

    // Per-pair-junction probability: explicit input or default 0.5.
    auto junc_p = [&](std::size_t pair_idx) -> float {
        if (junction_probs.size() == input.size() - 1
            && pair_idx < junction_probs.size()) {
            return junction_probs[pair_idx];
        }
        return 0.5f;
    };

    // Greedy walk: open a SplicedChain at sub_chain 0; for each next
    // sub_chain, decide merge-or-flush.
    SplicedChain current;
    current.sub_chains.push_back(input[0]);
    current.total_score = input[0].score;
    current.strand = input[0].strand;
    current.ref_id = input[0].ref_id;

    for (std::size_t i = 1; i < input.size(); ++i) {
        const auto& a = current.sub_chains.back();
        const auto& b = input[i];
        const float p = junc_p(i - 1);

        if (GapLooksLikeIntron(a, b, cfg)
            && p >= cfg.min_junction_probability) {
            // Merge.
            Junction j;
            j.donor_ref_pos    = a.ref_end - 1;
            j.acceptor_ref_pos = b.ref_start;
            j.query_gap        = b.query_start - a.query_end;
            j.probability      = p;
            j.is_confirmed     = p >= 0.50f;
            current.junctions.push_back(j);
            current.sub_chains.push_back(b);
            current.total_score += b.score;
        } else {
            // Flush current and start a new chain with the next sub.
            // If `current` is a singleton we count it; otherwise it
            // remains a multi-chain (and we still count the merge to
            // 1 chain output but not as singleton).
            if (current.sub_chains.size() == 1) {
                ++out.n_singletons_kept;
            }
            out.chains.push_back(std::move(current));
            current = SplicedChain{};
            current.sub_chains.push_back(b);
            current.total_score = b.score;
            current.strand = b.strand;
            current.ref_id = b.ref_id;
        }
    }
    // Push the final current.
    if (current.sub_chains.size() == 1) {
        ++out.n_singletons_kept;
    }
    out.chains.push_back(std::move(current));

    return out;
}

// ---------------------------------------------------------------------------
// EmitSplicedCigar — concatenate sub-chain CIGARs with N or D between.
// ---------------------------------------------------------------------------

namespace {

// Tokenise a CIGAR string into (length, op) pairs.
std::vector<std::pair<std::uint64_t, char>> ParseCigar(const std::string& s) {
    std::vector<std::pair<std::uint64_t, char>> ops;
    std::uint64_t num = 0;
    for (const char c : s) {
        if (c >= '0' && c <= '9') {
            num = num * 10 + static_cast<std::uint64_t>(c - '0');
        } else {
            ops.emplace_back(num, c);
            num = 0;
        }
    }
    return ops;
}

}  // namespace

std::string EmitSplicedCigar(const SplicedChain& sc) {
    const auto& subs = sc.sub_chains;
    const std::size_t n = subs.size();
    std::stringstream ss;
    for (std::size_t i = 0; i < n; ++i) {
        auto ops = ParseCigar(subs[i].cigar);
        // The trailing soft-clip of a non-last sub-chain and the leading
        // soft-clip of a non-first sub-chain represented the neighbouring exon
        // — now encoded as the intron N. Drop them, else they become illegal
        // INTERNAL soft-clips in the merged CIGAR. The leading clip of the
        // first sub-chain and the trailing clip of the last are real read ends
        // and are kept.
        if (i > 0 && !ops.empty() && ops.front().second == 'S') {
            ops.erase(ops.begin());
        }
        if (i + 1 < n && !ops.empty() && ops.back().second == 'S') {
            ops.pop_back();
        }
        for (const auto& [len, op] : ops) ss << len << op;

        if (i + 1 < n) {
            const auto& a = subs[i];
            const auto& b = subs[i + 1];
            // Residual read gap between the two exon blocks (their aligned ends
            // didn't meet exactly) → an insertion, so the total query length
            // stays exact (an N consumes no query).
            if (b.query_start > a.query_end) {
                ss << (b.query_start - a.query_end) << 'I';
            }
            const std::uint64_t intron_len =
                (b.ref_start > a.ref_end) ? (b.ref_start - a.ref_end) : 0;
            if (intron_len > 0) {
                const char op = (i < sc.junctions.size() &&
                                 sc.junctions[i].is_confirmed) ? 'N' : 'D';
                ss << intron_len << op;
            }
        }
    }
    return ss.str();
}

}  // namespace llmap::mapping
