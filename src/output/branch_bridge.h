// LLmap — BRANCH bridge: GAF-shape output compatible with BRANCH
// v0.9 Reconciliation engine.
//
// BRANCH ([[branch_tooling_state]]) consumes a slightly different GAF
// flavour than vg: it requires per-read provenance + cluster_id +
// a phasing tag at the start of the path column. This writer emits
// rows in that format so LLmap output can plug straight into BRANCH
// dual-phaser v0.9 pipelines.
//
// Format (tab-separated):
//
//   1. read_id
//   2. cluster_id (or '-' when no cluster)
//   3. haplotype  ('0', '1', or '.')
//   4. path       (BRANCH-style: starts with '@' + bubble id when
//                   read traverses a bubble, else '|' + chrom:start)
//   5. score
//   6. n_supporting_reads (set to 1 unless writer aggregates)
//   7. provenance ("LLMAP" by default; writer can override)

#pragma once

#include <cstdint>
#include <filesystem>
#include <span>
#include <string>

namespace llmap::output {

struct BranchGafRow {
    std::string read_id;
    std::uint32_t cluster_id{0};   // 0 ⇒ no cluster
    char haplotype{'.'};
    std::string path;
    std::int32_t score{0};
    std::uint32_t n_supporting_reads{1};
    std::string provenance{"LLMAP"};
};

[[nodiscard]] std::string EmitBranchGafRow(const BranchGafRow& row);

[[nodiscard]] bool WriteBranchGafFile(const std::filesystem::path& path,
                                       std::span<const BranchGafRow> rows);

}  // namespace llmap::output
