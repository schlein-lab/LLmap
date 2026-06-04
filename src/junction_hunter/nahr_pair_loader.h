// LLmap — junction_hunter: NAHR-pair definition loader.
//
// Loads pair definitions from the canonical TSV produced by Christian's
// genome-wide NAHR-linter (`genome_wide_nahr_2kb_pairs.tsv`). The TSV
// has a header row followed by one record per line; columns are TAB-
// separated. We tolerate extra trailing columns and unknown column
// orderings as long as the required fields are present.

#pragma once

#include "junction_hunter/junction_hunter_types.h"

#include <filesystem>
#include <vector>

namespace llmap::junction_hunter {

struct PairLoaderStatus {
    bool ok{false};
    std::size_t records_loaded{0};
    std::size_t records_skipped{0};
    std::string error;
};

/// Load NAHR-pair definitions from the linter TSV.
///
/// Required columns (any order, header-driven):
///   pair_id, chrom, lcr_up_start, lcr_up_end,
///   lcr_down_start, lcr_down_end, interior_start, interior_end,
///   identity
/// Optional column:
///   interior_len_bp / interior_kb
PairLoaderStatus LoadNahrPairsTsv(const std::filesystem::path& tsv,
                                  std::vector<NahrPair>& out);

}  // namespace llmap::junction_hunter
