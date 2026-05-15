// LLmap -- region annotation types.
//
// The "region annotation" layer attaches per-base mapping-parameter
// overrides to a reference. Each interval can override any of the
// chain-DP parameters (k, max_occ, lambda_scale, identity threshold,
// etc.). Lookup is O(log N) per position via the interval tree.
//
// Layer priority (highest wins on overlap):
//   1. agent decision     (--llm=auto runtime hits)
//   2. specific_loci      (per-locus database)
//   3. taxonomy           (universal classifier output)
//   4. default            (taxonomy says "unique_single_copy")

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace llmap::annot {

// Reference position
struct Position {
    uint32_t ref_id = 0;
    uint32_t pos = 0;
};

// Parameter overrides for a position interval.
// Every field is optional. Unset = "inherit from default / lower layer".
struct ParamOverride {
    std::optional<uint8_t> k;                     // minimizer k
    std::optional<uint8_t> w;                     // minimizer w
    std::optional<uint32_t> max_occ;              // index max_occ cap
    std::optional<float> lambda_scale;            // 1.0 = default
    std::optional<float> identity_threshold;      // 0.0 - 1.0
    std::optional<float> anchor_weight_scale;     // scales w_a in chain score
    std::optional<bool> report_multi_position;    // emit top-N chains
    std::optional<bool> require_psv_disambig;     // need PSV catalog
    std::optional<bool> allow_high_mismatch;      // for hypervariable regions
    std::optional<bool> require_llm_at_runtime;   // force agent call

    // Compose self over base: keep our set fields, fall back to base for unset.
    void OverlayOver(const ParamOverride& base);
};

// Source layer that produced this annotation
enum class AnnotationLayer : uint8_t {
    Default       = 0,  // hard-coded baseline
    Taxonomy      = 1,  // classifier output
    SpecificLocus = 2,  // per-locus database
    AgentDecision = 3,  // live LLM call
    Stochastic    = 4,  // online empirical update
};

// One annotation interval: [start, end) on a given ref contig.
struct AnnotationInterval {
    uint32_t ref_id;
    uint32_t start;
    uint32_t end;
    std::string region_name;     // e.g. "tandem_repeat_alpha_satellite"
    std::string source;          // e.g. "taxonomy" or "specific_locus:chr1_centromere"
    AnnotationLayer layer = AnnotationLayer::Default;
    ParamOverride params;
};

// Local features extracted per window before classification.
struct WindowFeatures {
    uint32_t ref_id = 0;
    uint32_t start = 0;
    uint32_t end = 0;
    float shannon_5mer = 0.0f;
    float gc_content = 0.0f;
    int32_t tandem_period = -1;       // -1 = none detected
    uint32_t kmer_multiplicity_p95 = 1;
    float palindrome_density = 0.0f;
    float orf_density = 0.0f;
    std::string consensus_match;      // empty if none
};

}  // namespace llmap::annot
