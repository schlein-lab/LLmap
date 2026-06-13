// LLmap — Cluster → mini-contig (assemble-then-map, Block idea 2), Line A.
//
// Operator idea: Stage-1 maps reads against each other (self-interference) →
// clusters of reads from the same locus/transcript. Instead of mapping every
// read, build ONE virtual construct per cluster and map it to the bucket-probable
// region — far fewer units (speed) and a longer/more-unique unit (recall): the
// construct spans the UNION of partial reads → longer than any single read →
// dissolves the molecule/0 91 bp-fragment problem. Reads inherit the placement
// via member-propagation (Phase 3.5), each with its own offset/strand.
//
// LINE A (operator-confirmed 2026-06-13) — the construct is a TRANSIENT SEARCH
// PROBE: a smoothed representative sequence used ONLY to find the region. It is
// thrown away after mapping and never persisted. The LOSSLESS TRUTH is the layout
// (members[] offset/strand + per-position depth) over the caller's real reads —
// no base is averaged away in the output, multi-allele columns stay visible
// (heteroplasmy / low-VAF somatic / A→I editing / paralog PSVs preserved).
//   * probe          → region-finding only, NOT truth, never written out.
//   * members/depth  → 100% sequence + coverage + strand, the actual output.
// If the probe maps but the fit is poor, the caller re-buckets (escalation);
// nothing is ever dropped.
//
// This module is DECOUPLED from clustering — it takes an already-formed cluster's
// reads (+ their strand from minimizer_cluster) and returns probe + layout. The
// map→propagate / re-bucket loop is the pipeline's job.

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::mapper {

// Placement of a cluster member inside the contig (for lossless propagation).
struct ConsensusMember {
    std::size_t  read_idx{0};   // index into the input read span
    std::int64_t offset{0};     // start of the read within the contig (0-based)
    bool         anchored{false}; // false ⇒ not anchored (caller maps it alone)
    bool         reverse{false};  // read is reverse-strand relative to the contig
};

struct ConsensusContig {
    // TRANSIENT search probe — region-finding ONLY. NOT lossless truth, never
    // persisted. The truth is members[] + depth[] over the caller's real reads.
    std::string probe;
    std::vector<ConsensusMember> members;  // one per input read (lossless layout)
    std::vector<std::uint32_t>   depth;    // per-position coverage over the contig
    std::size_t   length{0};               // contig length (== probe.size())
    std::uint32_t n_anchored{0};           // reads that contributed to the layout
};

struct ConsensusConfig {
    std::uint32_t k{15};                    // anchor k-mer length
    std::uint32_t min_shared_kmers{3};      // min shared anchors to place a read
};

// Assemble a cluster's reads into a probe + lossless layout. Reads are assumed
// from the same cluster. `reverse` (optional, from minimizer_cluster's
// ReadCluster.reverse; 1 ⇒ reverse strand) orients each read into the contig
// frame before layout; empty ⇒ all forward. A single read returns itself as the
// probe. Reads that can't be anchored stay un-anchored (caller maps them
// individually — lossless: nothing dropped).
[[nodiscard]] ConsensusContig AssembleConsensus(
    std::span<const std::string> reads, const ConsensusConfig& cfg = {},
    std::span<const std::uint8_t> reverse = {});

}  // namespace llmap::mapper
