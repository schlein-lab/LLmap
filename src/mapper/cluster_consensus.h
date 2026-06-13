// LLmap — Cluster → consensus mini-contig (assemble-then-map, Block idea 2).
//
// Operator idea: Stage-1 maps reads against each other (self-interference) →
// clusters of reads from the same locus/transcript. Instead of mapping every
// read, assemble each cluster into a CONSENSUS MINI-CONTIG (longer, error-
// averaged) and map the contig — far fewer units (speed) and longer/more-unique
// units (recall: a full-transcript contig carries all exons → maps as one
// spliced unit, dissolving the molecule/0 91 bp-fragment problem). Reads inherit
// the placement via the existing member-propagation (Phase 3.5).
//
// This module is the missing assembly step; it is DECOUPLED from the clustering
// (Stage-1 / FAISS) — it takes an already-formed cluster's read sequences and
// returns the consensus + per-read offsets. So it builds/tests without FAISS;
// the clustering (FAISS or a CPU fallback) and the map→propagate wiring are the
// pipeline's job.
//
// V1 algorithm (forward strand): pick the longest read as the backbone, anchor
// each read to it by the modal shared-k-mer offset, and take a per-column
// majority vote → consensus. V2 (later): reverse-strand members, full POA,
// quality-weighted columns.

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::mapper {

// Placement of a cluster member inside the consensus contig (for propagation).
struct MemberOffset {
    std::size_t read_idx{0};   // index into the input read span
    std::int64_t offset{0};    // start of the read within the consensus (0-based)
    bool anchored{false};      // false ⇒ could not be anchored (kept as singleton)
};

struct ConsensusContig {
    std::string sequence;                  // the assembled consensus
    std::vector<MemberOffset> members;     // one per input read
    std::uint32_t n_anchored{0};           // reads that contributed to the consensus
};

struct ConsensusConfig {
    std::uint32_t k{15};                   // anchor k-mer length
    std::uint32_t min_shared_kmers{3};     // min shared anchors to place a read
};

// Assemble a cluster's reads into a consensus mini-contig. Reads are assumed to
// be from the same cluster (similar, same strand for V1). A single read returns
// itself as the contig. Reads that can't be anchored stay un-anchored (the
// caller maps them individually — lossless: nothing dropped).
[[nodiscard]] ConsensusContig AssembleConsensus(
    std::span<const std::string> reads, const ConsensusConfig& cfg = {});

}  // namespace llmap::mapper
