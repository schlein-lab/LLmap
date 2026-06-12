// LLmap — Input sniffer: format + size-stats + Transcript-Mode resolution.
//
// Resolves `core::TranscriptMode::Auto` to a concrete mode from the
// input file, per docs/design/llmap_mode_architecture.md §3. The
// `core::TranscriptMode` enum (Block 9) promised this sniffer; this is
// its implementation.
//
// Deterministic, dependency-light (stdlib + llmap_core). No alignment
// pilot-pass in v1 (the doc's "≥10% long-N CIGAR" rule needs a real
// alignment pass — deferred; see TODO in ResolveMode). Detection uses:
//
//   1. format sniff via magic bytes / first line
//   2. for BAM/SAM: @PG / header program-token inspection
//   3. for FASTQ:   basename FLNC/isoseq/cdna/rna match
//   4. for FASTA:   sequence-length stats (median / N50 / n_seqs)
//
// Every resolution carries a human-readable `reason` string for the
// `[mode-detect]` log line.

#pragma once

#include "core/transcript_mode.h"

#include <cstdint>
#include <optional>
#include <string>

namespace llmap::io {

// ===========================================================================
// File format, sniffed from magic bytes / first line.
// ===========================================================================

enum class FileFormat : std::uint8_t {
    Unknown = 0,
    Fasta,
    Fastq,
    Sam,
    Bam,
};

[[nodiscard]] const char* FileFormatName(FileFormat f) noexcept;

// ===========================================================================
// FASTA sequence-length statistics (sampled).
// ===========================================================================

struct FastaStats {
    std::uint64_t n_seqs{0};              // sequences counted (whole file)
    std::uint64_t median_len{0};          // over the sampled lengths
    std::uint64_t n50{0};                 // over the sampled lengths
    std::uint64_t sampled{0};             // how many lengths went into stats
    bool          sampled_truncated{false};  // n_seqs > sample_limit
};

// ===========================================================================
// Full sniff result.
// ===========================================================================

struct SniffResult {
    FileFormat           format{FileFormat::Unknown};
    core::TranscriptMode mode{core::TranscriptMode::GenomeReads};
    std::string          reason;          // for the [mode-detect] log line
    std::optional<FastaStats> fasta_stats;
};

// ===========================================================================
// API
// ===========================================================================

/// Sniff the file format from the first bytes / first line.
/// Reads only a small prefix. Returns Unknown on read error / empty file.
[[nodiscard]] FileFormat SniffFormat(const std::string& path);

/// Compute FASTA length stats. `n_seqs` counts every `>` record in the
/// file (cheap full scan of header lines); median/N50 are over at most
/// `sample_limit` sequence lengths (the first ones encountered).
[[nodiscard]] FastaStats ComputeFastaStats(const std::string& path,
                                           std::uint64_t sample_limit = 1000);

/// Resolve a (possibly Auto) mode to a concrete one.
///
/// - `mode_override != Auto`            → returned verbatim (override wins).
/// - `has_reads && has_assembly`        → ReadsVsAssembly.
/// - otherwise sniff `primary_path` and apply the §3 heuristic.
///
/// `primary_path` is the reads path (or assembly path / first positional)
/// — whichever the CLI treats as the primary input.
[[nodiscard]] SniffResult ResolveMode(const std::string& primary_path,
                                      core::TranscriptMode mode_override,
                                      bool has_reads,
                                      bool has_assembly);

}  // namespace llmap::io
