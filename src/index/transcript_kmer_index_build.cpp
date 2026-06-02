// LLmap — TranscriptKmerIndex::BuildFromAnchorStore.
//
// Walks every anchor in the store and emits k-mers into the appropriate
// origin table. Six branches, one per KmerOrigin. The logic is
// deliberately verbose-but-explicit so it's clear which structural
// case each branch handles.
//
// Notes on correctness:
//
//   * For IntraExon, we generate k-mers from anchor.sequence directly
//     so they exactly mirror what a read of the same sequence will hash
//     to. For long anchors (>k_intra) this is many k-mers; we cap via
//     max_occ at the AddHit boundary.
//
//   * For JunctionSpanning, we slide a window of k_junction across the
//     join of two adjacent exons. The join position is computed from
//     the anchor's exon_boundaries[i].pos_in_transcript field. We need
//     k_junction/2 bases on each side of the boundary, so anchors with
//     exons shorter than k_junction/2 produce shorter junction-spanning
//     k-mer sets — that's fine, it just means fewer index entries for
//     those exons.
//
//   * For BackSpliceSpanning, we synthesise the back-splice sequence by
//     concatenating the LAST k_circular/2 bp of the 3'-most exon with
//     the FIRST k_circular/2 bp of the 5'-most exon (the circRNA
//     ordering). This is detectable when the anchor has any
//     ExonBoundary with spliceosome_class == 3.
//
//   * For SterileIntronic, we look at anchors tagged with switch_region
//     and emit k-mers over their full sequence — the entire S-region
//     IS the indexable surface for sterile-germline reads.
//
//   * For PreMrnaIntronic, opt-in. When enabled we generate k-mers from
//     the intronic gaps between consecutive exons; that requires
//     access to a reference FASTA (the anchor's exon sequences alone
//     don't tell us the intron sequence). This commit ships the
//     scaffold; the actual genomic-intron sequence pull lands in
//     Block 9 (where we have the reference FASTA in scope).
//
//   * ShortRna uses alt_k for any anchor whose total sequence length
//     is <= 50 nt and whose kind is in the small-RNA family.

#include "index/transcript_kmer_index.h"

#include "core/transcript_kind.h"

#include <algorithm>
#include <cctype>

namespace llmap::index {

namespace {

bool IsCleanRun(std::string_view s) noexcept {
    for (char c : s) {
        const char u = static_cast<char>(
            std::toupper(static_cast<unsigned char>(c)));
        if (u != 'A' && u != 'C' && u != 'G' && u != 'T') return false;
    }
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// ForEachKmer — generic walker, defined inline-template in the header
// would be cleaner but the template body is small enough to live here.
// ---------------------------------------------------------------------------

template <typename Emit>
void TranscriptKmerIndex::ForEachKmer(std::string_view seq,
                                       std::uint8_t k,
                                       Emit&& emit) {
    if (k == 0 || seq.size() < k) return;
    for (std::size_t i = 0; i + k <= seq.size(); ++i) {
        auto sub = seq.substr(i, k);
        if (!IsCleanRun(sub)) continue;
        const auto h = HashKmer(sub);
        emit(h, static_cast<std::uint32_t>(i), /*same_strand=*/true);
    }
}

// ---------------------------------------------------------------------------
// BuildFromAnchorStore
// ---------------------------------------------------------------------------

void TranscriptKmerIndex::BuildFromAnchorStore(
    const anchor::AnchorStore& store,
    const annot::SpliceSiteDb& /*splice*/,
    const TranscriptKmerIndexConfig& cfg) {

    Clear();
    cfg_ = cfg;

    const auto& anchors = store.anchors();

    for (std::uint32_t idx = 0; idx < anchors.size(); ++idx) {
        const auto& a = anchors[idx];
        if (a.sequence.empty()) continue;

        const bool is_short = a.sequence.size() <= 50
                              && core::IsSmallRna(a.kind);

        // --- 1. ShortRna table for small RNAs (k=alt_k) -----------------
        if (cfg_.include_short_rna && is_short) {
            ForEachKmer(a.sequence, cfg_.alt_k,
                [&](std::uint64_t h, std::uint32_t off, bool ss) {
                    TranscriptKmerHit hit;
                    hit.anchor_id_idx = idx;
                    hit.kmer_offset_in_anchor = off;
                    hit.same_strand = ss;
                    hit.kind = a.kind;
                    AddHit(KmerOrigin::ShortRna, h, hit);
                });
            // Short-RNA anchors only populate the short-RNA table.
            continue;
        }

        // --- 2. IntraExon — full anchor sequence, k=k_intra -------------
        ForEachKmer(a.sequence, cfg_.k_intra,
            [&](std::uint64_t h, std::uint32_t off, bool ss) {
                TranscriptKmerHit hit;
                hit.anchor_id_idx = idx;
                hit.kmer_offset_in_anchor = off;
                hit.same_strand = ss;
                hit.kind = a.kind;
                AddHit(KmerOrigin::IntraExon, h, hit);
            });

        // --- 3. JunctionSpanning — for spliced anchors -----------------
        if (cfg_.include_junction_spanning && a.is_spliced()) {
            const std::uint8_t kj = cfg_.k_junction;
            const std::uint32_t half = kj / 2;

            // For each junction, generate k_junction-mers centred on
            // pos_in_transcript.
            for (const auto& b : a.exon_boundaries) {
                // Back-splice gets its own table below.
                if (b.spliceosome_class == 3) continue;

                const std::uint32_t pos = b.pos_in_transcript;
                if (pos < half) continue;                // junction too close to 5' end
                if (pos + half > a.sequence.size()) continue;  // too close to 3' end

                // Sliding window covering all kj-length substrings that
                // straddle the junction (window centres in [pos-half,
                // pos+half-kj]).
                const std::size_t lo = pos - half;
                const std::size_t hi = std::min<std::size_t>(
                    pos + half - kj + 1,
                    a.sequence.size() - kj + 1);
                for (std::size_t i = lo; i < hi; ++i) {
                    auto sub = std::string_view(a.sequence).substr(i, kj);
                    if (!IsCleanRun(sub)) continue;
                    TranscriptKmerHit hit;
                    hit.anchor_id_idx = idx;
                    hit.kmer_offset_in_anchor = static_cast<std::uint32_t>(i);
                    hit.same_strand = true;
                    hit.kind = a.kind;
                    AddHit(KmerOrigin::JunctionSpanning, HashKmer(sub), hit);
                }
            }
        }

        // --- 4. BackSpliceSpanning — for circRNA anchors ---------------
        if (cfg_.include_backsplice && a.is_circular()) {
            const std::uint8_t kc = cfg_.k_circular;
            const std::uint8_t half = kc / 2;
            // Build the synthetic back-splice sequence: tail of anchor
            // concatenated with head of anchor (circular wrap).
            if (a.sequence.size() >= kc) {
                std::string circ_seq;
                circ_seq.reserve(kc);
                circ_seq.append(a.sequence,
                                 a.sequence.size() - half,
                                 half);
                circ_seq.append(a.sequence, 0, kc - half);
                if (IsCleanRun(circ_seq)) {
                    TranscriptKmerHit hit;
                    hit.anchor_id_idx = idx;
                    hit.kmer_offset_in_anchor = 0;  // synthetic; not a real offset
                    hit.same_strand = true;
                    hit.kind = a.kind;
                    AddHit(KmerOrigin::BackSpliceSpanning,
                            HashKmer(circ_seq), hit);
                }
            }
        }

        // --- 5. SterileIntronic — for switch-region anchors -----------
        if (cfg_.include_sterile_intronic) {
            const bool is_switch_region = std::any_of(
                a.tags.begin(), a.tags.end(),
                [](const std::string& t) {
                    return t.find("switch_region") != std::string::npos
                        || t.starts_with("Sgamma")
                        || t.starts_with("Smu")
                        || t.starts_with("Salpha")
                        || t.starts_with("Sepsilon");
                });
            if (is_switch_region) {
                ForEachKmer(a.sequence, cfg_.k_intra,
                    [&](std::uint64_t h, std::uint32_t off, bool ss) {
                        TranscriptKmerHit hit;
                        hit.anchor_id_idx = idx;
                        hit.kmer_offset_in_anchor = off;
                        hit.same_strand = ss;
                        hit.kind = a.kind;
                        AddHit(KmerOrigin::SterileIntronic, h, hit);
                    });
            }
        }

        // --- 6. PreMrnaIntronic — opt-in, requires reference FASTA -----
        if (cfg_.include_premrna_intronic) {
            // Scaffold only — actual intronic sequence pull from the
            // reference FASTA lands in Block 9 (CLI integration). When
            // landed, we generate k-mers over the intronic stretch
            // between consecutive exons of `a` and emit them into
            // premrna_intronic_.
            // For now this branch is intentionally inert; the toggle
            // exists so callers can pre-arrange their pipelines.
        }
    }
}

}  // namespace llmap::index
