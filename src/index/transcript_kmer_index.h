// LLmap — TranscriptKmerIndex: exon-boundary-aware reverse k-mer index.
//
// Standard k-mer indices treat the genome as a flat sequence; that's
// fine for DNA mapping but it produces nonsense junction matches for
// transcript data. A k-mer that spans the exon-exon junction in
// cDNA has no contiguous occurrence in genomic DNA — yet a naïve
// genome-built index will happily return false hits at intronic
// positions where the bases just happen to match.
//
// TranscriptKmerIndex flips the direction (transcripts → genome, not
// genome → reads) AND segregates k-mers by their structural origin:
//
//   IntraExon          — entirely within a single exon. Safe to match
//                        against both genomic windows and reads.
//   JunctionSpanning   — overlaps an exon-exon junction in cDNA. Only
//                        valid against READS (which are also spliced);
//                        never queried against genomic windows.
//   BackSpliceSpanning — circRNA back-splice junction. Only valid
//                        against circRNA reads; ordered acceptor→donor.
//   SterileIntronic    — intronic w.r.t. mature mRNA but present in
//                        sterile germline transcripts (Iγ→Sγ→Cγ).
//                        Same shape as a normal k-mer but flagged so
//                        it doesn't match general genomic positions.
//   PreMrnaIntronic    — intronic in mature mRNA but expected in
//                        pre-mRNA reads. OPT-IN only (catalog can grow
//                        large).
//   ShortRna           — separate table for k≤22 (miRNA / piRNA /
//                        siRNA-band). Lookup uses the alt_k size.
//
// Six separate unordered_maps means the lookup APIs can restrict to the
// origins they trust:
//
//   QueryGenomeWindow  — uses ONLY {IntraExon, PreMrnaIntronic}
//   QueryReadKmer      — uses ALL six tables
//
// Caller picks. Default Config has PreMrnaIntronic off (off-by-default
// to avoid index bloat) — turn it on with include_premrna_intronic.

#pragma once

#include "anchor/anchor_record.h"
#include "anchor/anchor_store.h"
#include "annot/splice_site_db.h"
#include "core/transcript_kind.h"

#include <cstdint>
#include <filesystem>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::index {

// ===========================================================================
// Per-k-mer origin classification.
// ===========================================================================

enum class KmerOrigin : std::uint8_t {
    IntraExon = 0,
    JunctionSpanning,
    BackSpliceSpanning,
    SterileIntronic,
    PreMrnaIntronic,
    ShortRna,
};

const char* KmerOriginName(KmerOrigin o) noexcept;

// ===========================================================================
// Hit record.
// ===========================================================================

struct TranscriptKmerHit {
    /// Index into AnchorStore::anchors() — call store.anchors()[idx] to
    /// dereference. Stable while the store has not been mutated.
    std::uint32_t anchor_id_idx{0};
    /// Offset of the k-mer start INSIDE the anchor sequence (0-based).
    std::uint32_t kmer_offset_in_anchor{0};
    /// True if the k-mer matched on the same strand as the anchor.
    bool same_strand{true};
    KmerOrigin origin{KmerOrigin::IntraExon};
    /// Mirror of the anchor's TranscriptKind for fast filtering at
    /// lookup time (avoids dereferencing the AnchorRecord for cheap
    /// origin-based pre-filtering).
    core::TranscriptKind kind{core::TranscriptKind::Unknown};
};

// ===========================================================================
// Configuration.
// ===========================================================================

struct TranscriptKmerIndexConfig {
    /// k for intra-exon k-mers. Long enough to be specific in repeat-
    /// rich regions but not too long to over-fit short exons.
    /// Default 51 — see Block-3 design notes in the plan.
    std::uint8_t k_intra{51};

    /// k for junction-spanning k-mers. Shorter than k_intra because
    /// many exons are <50 bp (Hinge region is 33 bp); a longer k
    /// wouldn't fit across.
    std::uint8_t k_junction{31};

    /// k for back-splice (circRNA) junctions. Same as junction.
    std::uint8_t k_circular{31};

    /// Fallback k for small RNAs (miRNA, piRNA, siRNA). Must be < 24
    /// since mature small RNAs are ~21-31 nt total.
    std::uint8_t alt_k{21};

    /// Max occurrences per hash before we cap. Keeps repeat-rich IGH
    /// loci from blowing memory. 0 = no cap.
    std::uint32_t max_occ{200};

    /// Per-origin enable flags.
    bool include_junction_spanning{true};
    bool include_backsplice{true};
    bool include_sterile_intronic{true};
    bool include_premrna_intronic{false};   // off by default
    bool include_short_rna{true};
};

// ===========================================================================
// The index.
// ===========================================================================

class TranscriptKmerIndex {
public:
    TranscriptKmerIndex() = default;

    /// Build the index over an AnchorStore. SpliceSiteDb is used only
    /// when we need to verify motif consistency for borderline cases;
    /// the AnchorStore itself already carries exon_boundaries.
    ///
    /// Idempotent — calls Clear() first.
    void BuildFromAnchorStore(const anchor::AnchorStore& store,
                              const annot::SpliceSiteDb& splice,
                              const TranscriptKmerIndexConfig& cfg = {});

    /// Reverse direction: a genomic window (e.g. one chrom-region of
    /// the reference) is hashed and matched against the transcript
    /// k-mers we hold. Only IntraExon + PreMrnaIntronic tables are
    /// consulted — junction-spanning k-mers cannot validly match the
    /// genome.
    [[nodiscard]] std::vector<TranscriptKmerHit>
    QueryGenomeWindow(std::string_view window_seq) const;

    /// Forward direction: a read k-mer hash is looked up across ALL
    /// origin tables. Caller can post-filter by origin / kind if needed.
    [[nodiscard]] std::vector<TranscriptKmerHit>
    QueryReadKmer(std::uint64_t kmer_hash) const;

    /// Convenience: hash + query a string-form k-mer.
    [[nodiscard]] std::vector<TranscriptKmerHit>
    QueryReadKmer(std::string_view kmer_seq) const;

    /// Per-table counts, used by the lossless summary.
    [[nodiscard]] std::size_t TableSize(KmerOrigin o) const noexcept;

    /// Total k-mers across all tables.
    [[nodiscard]] std::size_t TotalKmers() const noexcept;

    void Clear();

    // ----- Serialisation (optional cache) ------------------------------
    bool Save(const std::filesystem::path& path) const;
    bool Load(const std::filesystem::path& path);

private:
    /// Canonical k-mer hashing — same primitive used by the existing
    /// MinimizerIndex (so a read's hash from either tool agrees).
    static std::uint64_t HashKmer(std::string_view k) noexcept;

    /// Generate every k-mer of length k from `seq` and call `emit(hash,
    /// offset, same_strand)` for each one. Skips windows with non-ACGT
    /// bases (we deliberately don't try to decode IUPAC codes — the
    /// caller's anchor sequences should already be clean).
    template <typename Emit>
    static void ForEachKmer(std::string_view seq,
                             std::uint8_t k,
                             Emit&& emit);

    /// Add a single (hash, hit) record into the named table, honouring
    /// max_occ cap.
    void AddHit(KmerOrigin origin,
                 std::uint64_t hash,
                 TranscriptKmerHit hit);

    // Six segregated hash tables — see header docs for why each exists.
    std::unordered_map<std::uint64_t,
                       std::vector<TranscriptKmerHit>> intra_exon_;
    std::unordered_map<std::uint64_t,
                       std::vector<TranscriptKmerHit>> junction_spanning_;
    std::unordered_map<std::uint64_t,
                       std::vector<TranscriptKmerHit>> backsplice_;
    std::unordered_map<std::uint64_t,
                       std::vector<TranscriptKmerHit>> sterile_intronic_;
    std::unordered_map<std::uint64_t,
                       std::vector<TranscriptKmerHit>> premrna_intronic_;
    std::unordered_map<std::uint64_t,
                       std::vector<TranscriptKmerHit>> short_rna_;

    TranscriptKmerIndexConfig cfg_;
};

}  // namespace llmap::index
