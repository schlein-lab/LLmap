// LLmap — SegDup Catalog public API.
//
// Loads curated T1 (full-evidence) and bulk T2 (coords-only) segmental-
// duplication catalog entries from JSON / JSONL and provides coordinate
// and locus-id lookup.
//
// Layout follows catalog/schema/curated.schema.json (schema_version 0.1+,
// forward-compatible with 0.2 changes — unknown fields are tolerated).
//
// Namespace: llmap::catalog
// Module:    src/catalog/

#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::catalog {

/// Genomic interval keyed by assembly name and chromosome.
///
/// `start` is 0-based inclusive (BED-style), `end` is exclusive.
/// `strand` is '+', '-' or '.' for unknown.
struct GenomicCoords {
    std::string assembly;   ///< e.g. "GRCh38", "CHM13", "HPRC_pangenome_v1.1"
    std::string chrom;      ///< e.g. "chr14"
    std::int64_t start = 0; ///< 0-based inclusive
    std::int64_t end = 0;   ///< exclusive
    char strand = '.';      ///< '+', '-', or '.'
    std::string note;       ///< free-form annotation from JSON

    bool valid() const { return !assembly.empty() && !chrom.empty() && end > start; }
};

/// A diagnostic SNP from `diagnostic_features.discriminating_snps`.
///
/// Only the fields LLmap needs for mapping decisions are captured here;
/// remaining JSON content (provenance, frequencies in subpopulations,
/// notes) stays in the raw entry but is not surfaced through this struct.
struct DiagnosticSnp {
    std::string id;                ///< e.g. "IGHG4_CH1_pos69"
    std::string transcript;        ///< e.g. "IGHG4-201"
    std::int64_t cds_offset = -1;  ///< 1-based CDS offset, or -1 if absent
    std::int64_t grch38_pos = -1;  ///< 1-based GRCh38 chr position, or -1
    std::string chrom;             ///< e.g. "chr14" (derived from position key)
    std::string ref;               ///< reference allele
    std::string alt;               ///< alternate allele
    std::string consequence;       ///< e.g. "p.Ser7=", "synonymous"
    std::string snp_class;         ///< e.g. "ORIG_polymorph", "DUP_fixed"
    double freq_canonical = -1.0;  ///< -1 if absent in JSON
    double freq_dup = -1.0;        ///< -1 if absent in JSON
};

/// Primary mapping-strategy parameters (from `mapping_strategy.primary`).
struct MappingPrimary {
    std::int32_t kmer_size = 0;          ///< 0 if unspecified
    std::int32_t max_mismatch = 0;       ///< 0 if unspecified
    std::int32_t include_flanking_bp = 0;///< 0 if unspecified
    std::string include_flanking_anchor; ///< e.g. "CH1_start"
    bool require_unique_chain = false;
    std::string note;
};

/// One stage of `mapping_strategy.fallback_chain`.
struct MappingFallbackStage {
    std::int32_t stage = 0;
    std::string name;       ///< e.g. "relaxed_mismatch", "llm_checkpoint"
    std::string trigger;    ///< free-form description
    std::string rationale;
    // Optional fields kept as JSON-blob is overkill; we expose only the
    // common ones used by the mapper:
    std::optional<std::int32_t> kmer_size;
    std::optional<std::int32_t> max_mismatch;
    std::optional<std::int32_t> top_k;
    std::optional<bool> use_extension;
    std::optional<bool> emit_warning;
    std::optional<std::string> opt_in_flag;
};

/// Catalog `diagnostic_features.promoter_signature` block.
///
/// `canonical_motif` / `duplicate_motif` may contain a literal `...` meaning
/// "ordered subsequence within `window_relative_to_anchor`" (e.g.
/// `"CGGTTCTT...GTCTATCTGCGATGG"` requires prefix then suffix in 5'->3'
/// order, no spacer-length constraint). Match semantics live in the IGH-mode
/// classifier (`src/classify/`), the catalog only stores the raw strings.
struct PromoterSignature {
    std::string window_relative_to_anchor;     ///< e.g. "-100..0"
    std::string anchor;                        ///< e.g. "IGHG4_CH1_start"
    std::optional<std::string> canonical_motif;
    std::optional<std::string> duplicate_motif;
    std::string note;                          ///< free-form
};

/// A T1 (curated) or T2 (bulk) catalog entry.
///
/// T2 bulk records populate only `locus_id`, `structural_architecture`,
/// `coords_by_assembly` and `tier`. T1 curated records additionally fill
/// `haplotype_class`, `discriminating_snps`, `mapping_primary`, and
/// `fallback_chain`.
struct SegDupCatalogEntry {
    /// Catalog tier of the entry (T1 = curated, T2 = bulk coords-only).
    enum class Tier { T1_Curated, T2_Bulk };

    std::string locus_id;                ///< unique identifier
    std::string human_name;              ///< human-readable name (T1 only)
    std::string version;                 ///< e.g. "v2026.Q2"
    std::string schema_version;          ///< e.g. "0.1", "0.2"

    std::string structural_architecture; ///< enum-like, see schema
    std::string haplotype_class;         ///< locus-specific (T1 only)
    std::vector<std::string> mechanism;        ///< T1 only
    std::vector<std::string> clinical_function;///< T1 only
    std::string nahr_status;                   ///< T1 only

    /// Coordinates indexed by assembly name. Multiple assemblies allowed.
    std::unordered_map<std::string, GenomicCoords> coords_by_assembly;

    /// Diagnostic SNPs (T1 only; empty for T2).
    std::vector<DiagnosticSnp> discriminating_snps;

    MappingPrimary mapping_primary;             ///< T1 only
    std::vector<MappingFallbackStage> fallback_chain; ///< T1 only

    /// Promoter signature from diagnostic_features (T1 only, optional).
    /// Filled when the JSON has `diagnostic_features.promoter_signature`.
    std::optional<PromoterSignature> promoter_signature;

    Tier tier = Tier::T1_Curated;
    std::filesystem::path source_path;   ///< origin file (debug aid)

    /// Convenience accessors.
    bool is_curated() const { return tier == Tier::T1_Curated; }
    bool is_bulk()    const { return tier == Tier::T2_Bulk; }
    std::size_t n_snps() const { return discriminating_snps.size(); }

    /// True iff `pos` falls within this entry's coordinates for `assembly`.
    /// Half-open interval [start, end). Position is 0-based.
    bool contains(std::string_view assembly,
                  std::string_view chrom,
                  std::int64_t pos) const;

    /// True iff [`start`, `end`) overlaps this entry's coordinates.
    bool overlaps(std::string_view assembly,
                  std::string_view chrom,
                  std::int64_t start,
                  std::int64_t end) const;
};

/// Result of loading a single file.
struct LoadStatus {
    bool ok = false;
    std::size_t records_loaded = 0;
    std::size_t records_skipped = 0;
    std::string error;          ///< empty when ok == true
};

/// In-memory SegDup catalog with coordinate and id indices.
///
/// Coordinate lookup is O(n) over entries-per-(assembly,chrom). For the
/// expected catalog size (≈100 T1 + 10k T2 records) a sorted-by-start
/// scan is sufficient; an interval tree can be added later without
/// changing this API.
class SegDupCatalog {
public:
    SegDupCatalog() = default;

    /// Load a single curated T1 JSON file.
    LoadStatus load_curated_file(const std::filesystem::path& path);

    /// Load every `*.json` in `dir` as curated T1 entries.
    /// Non-recursive. Returns aggregated counts.
    LoadStatus load_curated_dir(const std::filesystem::path& dir);

    /// Load a bulk T2 JSONL file. One JSON object per line.
    LoadStatus load_bulk_jsonl(const std::filesystem::path& path);

    /// Add an externally-built entry. Useful for tests and synthetic data.
    void add_entry(SegDupCatalogEntry entry);

    /// Total entry count.
    std::size_t size() const { return entries_.size(); }

    /// Borrowed view of all entries (stable while catalog is unmodified).
    const std::vector<SegDupCatalogEntry>& entries() const { return entries_; }

    /// Find first entry whose coords contain (assembly, chrom, pos).
    /// Returns std::nullopt if nothing matches.
    std::optional<SegDupCatalogEntry>
    lookup_by_coords(std::string_view assembly,
                     std::string_view chrom,
                     std::int64_t pos) const;

    /// Find all entries whose coords overlap [start, end).
    std::vector<SegDupCatalogEntry>
    lookup_overlapping(std::string_view assembly,
                       std::string_view chrom,
                       std::int64_t start,
                       std::int64_t end) const;

    /// Find every entry whose `locus_id` equals `id`.
    /// (Returns vector to allow future versioned duplicates; currently
    /// `locus_id` is expected to be unique per catalog.)
    std::vector<SegDupCatalogEntry>
    lookup_by_locus_id(std::string_view id) const;

    /// Count of entries matching a structural_architecture value.
    std::size_t count_architecture(std::string_view architecture) const;

private:
    std::vector<SegDupCatalogEntry> entries_;
    /// locus_id → indices into entries_
    std::unordered_map<std::string, std::vector<std::size_t>> by_locus_id_;
    /// (assembly + '\0' + chrom) → indices into entries_, sorted by start
    std::unordered_map<std::string, std::vector<std::size_t>> by_chrom_;

    void reindex_entry(std::size_t idx);
    static std::string chrom_key(std::string_view assembly,
                                 std::string_view chrom);
};

/// Parse one curated-JSON document already in memory.
/// Returns the entry on success; `err` populated on failure.
std::optional<SegDupCatalogEntry>
parse_curated_json(std::string_view json_text,
                   const std::filesystem::path& source_path,
                   std::string& err);

/// Parse one bulk-T2 line (single JSON object).
std::optional<SegDupCatalogEntry>
parse_bulk_jsonl_line(std::string_view json_text,
                      const std::filesystem::path& source_path,
                      std::string& err);

}  // namespace llmap::catalog
