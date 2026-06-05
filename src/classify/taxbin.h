// LLmap — Mode-6 Taxbin engine: deterministic taxonomic binding of un-human
// reads from per-species evidence + read-to-read cluster consensus.
//
// Mode-6 replaces fragile external classifiers (Kraken/BLAST) with an
// LLmap-native pipeline (LLmap-only mandate). The orchestration is:
//
//   1. Stage 1 (self_interference::AllpairPipeline): read-to-read similarity →
//      Leiden clusters. Reads that look alike are grouped BEFORE any reference
//      is consulted.
//   2. Cascade over a species panel: each reference (classical::ClassicalPipeline
//      or Stage2 WaveCollapse) yields a per-read confidence. Specific→broad
//      order; reads that bind strongly short-circuit, only the still-unbound
//      residual is carried to the next, broader reference ("keep mapping until
//      it finds something").
//   3. This engine: turn the per-read × per-species confidence matrix into an
//      explicit likelihood distribution (including a NOVEL/dark pseudo-species),
//      then COLLAPSE ambiguous reads onto their cluster's consensus species —
//      a read with weak individual signal inherits the species its read-to-read
//      cluster collectively supports.
//
// Everything here is deterministic and uncertainty-honest: the full likelihood
// vector is emitted (no hidden collapse), the hard call is argmax of a
// transparent posterior, and there are NO arbitrary point scores.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace llmap::classify {

// One reference in the cascade panel.
struct SpeciesPanelEntry {
    std::string label;       // species / taxon name, e.g. "EBV", "human", "mycoplasma"
    std::string index_path;  // LLmap reference index for this species
    int cascade_rank = 0;    // ascending: 0 = most specific/likely, consulted first
};

// Per-read evidence assembled from the cascade. per_read_species_conf[r][s] is
// the best alignment confidence (identity / posterior, in [0,1]) of read r
// against species s; 0 means "no hit / not consulted".
struct TaxbinInput {
    std::vector<std::string> read_ids;                       // length N
    std::vector<std::string> species_labels;                 // length S (panel order)
    std::vector<std::vector<float>> per_read_species_conf;   // [N][S]
    std::vector<std::uint32_t> cluster_ids;                  // [N] from Stage 1; empty = no clustering
};

struct TaxbinConfig {
    float bind_threshold = 0.80f;        // a read "binds" species s if conf >= this
    float cluster_weight = 0.50f;        // posterior = (1-w)*own + w*cluster_consensus
    float min_margin = 0.10f;            // need top1-top2 >= this for a confident call
    bool  enable_cluster_collapse = true;
};

// NOVEL is represented as the last entry of every likelihood vector
// (index == species_labels.size()).
struct TaxbinReadResult {
    std::string read_id;
    std::vector<float> likelihood;   // length S+1, normalized; last = NOVEL/dark
    int   top_species = -1;          // index into species_labels, or -1 == NOVEL
    float top_prob = 0.0f;
    float margin = 0.0f;             // top1 - top2 of the posterior
    bool  novel = false;             // bound to no reference species
    bool  by_cluster = false;        // call decided by cluster consensus, not own evidence
    std::uint32_t cluster_id = 0;
};

struct TaxbinClusterResult {
    std::uint32_t cluster_id = 0;
    std::size_t size = 0;
    int   consensus_species = -1;        // -1 == NOVEL
    float purity = 0.0f;                 // fraction of members whose own argmax == consensus
    std::vector<float> aggregate_likelihood;  // length S+1
};

struct TaxbinResult {
    std::vector<TaxbinReadResult> reads;
    std::vector<TaxbinClusterResult> clusters;
    std::size_t num_species = 0;
    std::size_t num_novel = 0;
    std::size_t num_collapsed_by_cluster = 0;
};

// Pure, deterministic engine. No I/O, no randomness — safe to unit-test.
class Taxbin {
public:
    explicit Taxbin(TaxbinConfig cfg = {}) : cfg_(cfg) {}

    TaxbinResult Run(const TaxbinInput& in) const;

    // Exposed for testing: per-read raw likelihood over (species..., NOVEL).
    // raw[s] = conf[s]; raw[NOVEL] = max(0, bind_threshold - max_s conf[s]).
    // Normalized to sum 1. If a read has no evidence, mass concentrates on NOVEL.
    static std::vector<float> RawLikelihood(const std::vector<float>& conf,
                                            float bind_threshold);

    const TaxbinConfig& Config() const { return cfg_; }

private:
    TaxbinConfig cfg_;
};

}  // namespace llmap::classify
