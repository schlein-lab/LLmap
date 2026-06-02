// LLmap — ExpressionDb: tissue × cell-type expression priors for the
// Multi-Signal Fusion engine.
//
// Sources (loaded one-by-one; partial sets are valid):
//
//   GTEx v8 bulk           — 54 tissues × ~17 k samples, median TPM
//                             per gene per tissue. Public via dbGaP /
//                             GTEx portal.
//   GTEx LR (Q3 2025)      — long-read PacBio Iso-Seq, ~88 samples,
//                             per-transcript abundance. Public.
//   Tabula Sapiens         — single-cell, ~480 k cells, organ × cell-type.
//                             Public via Zenodo.
//   HCA (Human Cell Atlas) — multi-tissue single-cell. Public via HCA portal.
//   Human Protein Atlas    — RNA-tissue + RNA-cell-type. Public via HPA.
//   Recount3               — 60 k+ harmonised SRA samples × tissue.
//                             Public via recount.bio.
//
// All loaders are best-effort: missing files / unknown columns are
// tolerated, partial loads are valid, and a missing tissue or
// cell-type in the lookup returns a neutral default rather than
// throwing.
//
// Plumbing into the Fusion engine:
//
//   - L_expression_prior(anchor, tissue) reads ExpectedTpm(tx_id, tissue)
//     and applies the sigmoid mapping spec'd in Plan-Block 4.5
//     (log10(TPM+1) → [0.10, 0.98]).
//   - L_depth_coverage(anchor) uses ExpectedTpm too — same lookup,
//     different downstream model (Negative-Binomial per-exon).
//   - L_barcode_context(read) uses ExpectedCellTypeFraction(tx_id,
//     tissue, cell_type) when cell_type is known.

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::fusion {

class ExpressionDb {
public:
    // ===== Loaders ========================================================

    /// GTEx v8 bulk TPM matrix in GCT format (gene_id × tissue → median TPM).
    /// First two lines are GCT header; data starts at line 3.
    bool LoadGtexBulkGct(const std::filesystem::path& gct);

    /// GTEx Long-Read transcript abundance (Parquet or TSV).
    /// Expected columns: transcript_id, tissue, tpm.
    bool LoadGtexLongReadTsv(const std::filesystem::path& tsv);

    /// Tabula Sapiens cell-type expression — TSV with columns
    /// (transcript_id, organ, cell_type, mean_expression).
    bool LoadTabulaSapiensTsv(const std::filesystem::path& tsv);

    /// Human Cell Atlas subset — same TSV shape as Tabula Sapiens.
    bool LoadHcaTsv(const std::filesystem::path& tsv);

    /// Human Protein Atlas RNA-tissue export (TSV).
    bool LoadHpaTsv(const std::filesystem::path& tsv);

    /// Recount3 harmonised tissue medians — TSV (gene_id × tissue → median TPM).
    bool LoadRecount3Tsv(const std::filesystem::path& tsv);

    // ===== Lookups =======================================================

    /// Expected TPM for a transcript in a given tissue. Returns -1.0
    /// when the (transcript, tissue) pair is unknown — callers should
    /// treat that as 'neutral' rather than 'zero expression'.
    [[nodiscard]] float ExpectedTpm(std::string_view transcript_id,
                                     std::string_view tissue) const;

    /// Expected cell-type fraction expressing the transcript above
    /// background, in the requested (tissue, cell_type) context.
    /// Returns -1.0 when unknown.
    [[nodiscard]] float ExpectedCellTypeFraction(
        std::string_view transcript_id,
        std::string_view tissue,
        std::string_view cell_type) const;

    /// Diversity of expression across all known tissues; higher entropy
    /// means the transcript is broadly expressed, lower entropy means
    /// tissue-specific. Returns 0.0 when no data.
    [[nodiscard]] float ExpressionEntropy(
        std::string_view transcript_id) const;

    /// Telemetry — number of (transcript, tissue) records loaded.
    [[nodiscard]] std::size_t TotalEntries() const noexcept;
    [[nodiscard]] std::size_t DistinctTranscripts() const noexcept;
    [[nodiscard]] std::size_t DistinctTissues() const noexcept;

    void Clear();

private:
    // (transcript_id, tissue) → TPM
    struct Key {
        std::string transcript_id;
        std::string tissue;
        bool operator==(const Key& o) const noexcept {
            return transcript_id == o.transcript_id && tissue == o.tissue;
        }
    };
    struct KeyHash {
        std::size_t operator()(const Key& k) const noexcept {
            return std::hash<std::string>{}(k.transcript_id)
                ^ (std::hash<std::string>{}(k.tissue) * 0x9E3779B97F4A7C15ULL);
        }
    };

    std::unordered_map<Key, float, KeyHash> tpm_;

    // (transcript_id, tissue, cell_type) → fraction expressing
    struct CtKey {
        std::string transcript_id;
        std::string tissue;
        std::string cell_type;
        bool operator==(const CtKey& o) const noexcept {
            return transcript_id == o.transcript_id
                && tissue == o.tissue
                && cell_type == o.cell_type;
        }
    };
    struct CtKeyHash {
        std::size_t operator()(const CtKey& k) const noexcept {
            return std::hash<std::string>{}(k.transcript_id)
                ^ (std::hash<std::string>{}(k.tissue) * 0x9E3779B97F4A7C15ULL)
                ^ (std::hash<std::string>{}(k.cell_type) * 0xc6a4a7935bd1e995ULL);
        }
    };
    std::unordered_map<CtKey, float, CtKeyHash> cell_type_frac_;
};

}  // namespace llmap::fusion
