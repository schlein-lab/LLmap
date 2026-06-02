// LLmap — Pangenome bridge: GAF (Graph Alignment Format) export.
//
// When `--pangenome-graph` is configured, LLmap writes a sidecar
// .gaf file in addition to its native BAM/Parquet output. GAF is the
// graph-aware alignment format used by vg / GraphAligner / Giraffe
// pipelines — exporting it lets downstream HPRC pangenome analyses
// consume LLmap output directly without lossy BAM->GAF transformation.
//
// GAF v1 column layout (we emit the 12 mandatory columns):
//
//   1.  query_name
//   2.  query_length
//   3.  query_start (0-based)
//   4.  query_end
//   5.  strand ('+' / '-')
//   6.  path                — comma-sep node ids in the pangenome graph;
//                              when a linear-genome anchor is exported
//                              we emit '>'+chrom+':'+start as a single
//                              "pseudo-node" so the GAF round-trips
//                              through vg without errors
//   7.  path_length
//   8.  path_start
//   9.  path_end
//   10. residue_matches
//   11. block_length
//   12. mapping_quality
//
// This module only declares the row + writer; the actual file I/O is
// driven by parquet_writer_transcript.cpp / bam_writer integration.

#pragma once

#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

namespace llmap::output {

struct GafRow {
    std::string query_name;
    std::uint32_t query_length{0};
    std::uint32_t query_start{0};
    std::uint32_t query_end{0};
    char strand{'+'};
    std::string path;                 // comma-sep node ids; see header
    std::uint64_t path_length{0};
    std::uint64_t path_start{0};
    std::uint64_t path_end{0};
    std::uint32_t residue_matches{0};
    std::uint32_t block_length{0};
    std::uint8_t mapping_quality{60};
};

/// Emit a GAF v1 row as a tab-separated string. No trailing newline.
[[nodiscard]] std::string EmitGafRow(const GafRow& row);

/// Write a vector of rows to disk, one per line. Returns true on
/// success. Existing file is overwritten.
[[nodiscard]] bool WriteGafFile(const std::filesystem::path& path,
                                  std::span<const GafRow> rows);

}  // namespace llmap::output
