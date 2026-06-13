// LLmap — Splice-site (GT/AG) boundary snapping implementation.

#include "mapping/splice_snap.h"

#include <cctype>
#include <cstdint>
#include <cstdlib>

namespace llmap::mapping {

namespace {

inline char Up(char c) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

inline bool Motif2(std::string_view ref, std::int64_t pos, char m0, char m1) {
    if (pos < 0 || static_cast<std::uint64_t>(pos) + 2 > ref.size()) return false;
    return Up(ref[static_cast<std::size_t>(pos)]) == m0 &&
           Up(ref[static_cast<std::size_t>(pos) + 1]) == m1;
}

}  // namespace

SpliceSnap SnapJunction(const LinearSubChain& a, const LinearSubChain& b,
                        std::string_view ref_seq, char strand,
                        std::uint32_t window) {
    SpliceSnap r;
    r.donor_ref_pos = a.ref_end;        // defaults = unsnapped boundary
    r.acceptor_ref_pos = b.ref_start;
    r.query_split = a.query_end;

    if (ref_seq.empty() || b.ref_start < a.ref_end) return r;

    // Canonical intron 2-mers on the FORWARD reference:
    //   '+' transcript : GT (donor, intron 5') ... AG (acceptor, intron 3')
    //   '-' transcript : CT ............................ AC   (revcomp)
    const bool fwd = (strand != '-');
    const char d0 = fwd ? 'G' : 'C', d1 = 'T';
    const char c0 = 'A', c1 = fwd ? 'G' : 'C';

    // For each candidate read split q, the upstream exon ends at q (its ref end,
    // the donor, follows by the same delta) and the downstream exon starts at q
    // (its ref start, the acceptor, follows). So the read stays contiguous —
    // q_gap collapses to zero. Pick the canonical site of smallest total
    // boundary displacement; in-gap placements (q in [a.query_end,b.query_start])
    // all share the minimal displacement (= the slop), so they win over reaching
    // outside the gap.
    const auto qe_a = static_cast<std::int64_t>(a.query_end);
    const auto qs_b = static_cast<std::int64_t>(b.query_start);
    std::int64_t qlo = qe_a - static_cast<std::int64_t>(window);
    std::int64_t qhi = qs_b + static_cast<std::int64_t>(window);
    if (qlo < 0) qlo = 0;

    std::int64_t best_disp = -1;
    for (std::int64_t q = qlo; q <= qhi; ++q) {
        const std::int64_t da = q - qe_a;                       // exon-a grows by da
        const std::int64_t db = qs_b - q;                       // exon-b grows by db
        const std::int64_t donor = static_cast<std::int64_t>(a.ref_end) + da;
        const std::int64_t accpt = static_cast<std::int64_t>(b.ref_start) - db;
        // Keep both exons non-empty and the intron positive.
        if (donor <= static_cast<std::int64_t>(a.ref_start)) continue;
        if (accpt >= static_cast<std::int64_t>(b.ref_end)) continue;
        if (accpt - donor < 2) continue;
        if (Motif2(ref_seq, donor, d0, d1) && Motif2(ref_seq, accpt - 2, c0, c1)) {
            const std::int64_t disp = std::llabs(da) + std::llabs(db);
            if (best_disp < 0 || disp < best_disp) {
                best_disp = disp;
                r.snapped = true;
                r.donor_ref_pos = static_cast<std::uint64_t>(donor);
                r.acceptor_ref_pos = static_cast<std::uint64_t>(accpt);
                r.query_split = static_cast<std::uint32_t>(q);
                r.motif_score = 1.0f;
            }
        }
    }
    return r;
}

}  // namespace llmap::mapping
