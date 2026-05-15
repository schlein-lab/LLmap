// LLmap -- local feature extraction over reference windows.
//
// All features are deterministic functions of the reference sequence alone.
// They form the input to the rule-based classifier (or, equivalently, to a
// trained ML classifier if the user supplies one). No LLM involved at this
// stage.

#pragma once

#include "annot/annot_types.h"

#include <span>
#include <string_view>
#include <vector>

namespace llmap::annot {

struct FeatureExtractorConfig {
    uint32_t window_bp = 1000;     // 1 kb default; pass --window-bp 100 for viral
    uint32_t step_bp = 1000;       // non-overlapping windows by default
    uint32_t entropy_kmer_k = 5;   // k for Shannon-entropy computation
    uint32_t tandem_min_period = 2;
    uint32_t tandem_max_period = 1000;
};

// Compute features for every window of a single contig.
//
// seq.size() must equal contig length. ref_id is just propagated into
// WindowFeatures.ref_id for downstream bookkeeping.
std::vector<WindowFeatures> ExtractFeatures(
    uint32_t ref_id,
    std::string_view seq,
    const FeatureExtractorConfig& cfg = {});

// Compute features for every contig in a reference, given (name, seq) pairs.
// ref_id is the index into the input span.
std::vector<WindowFeatures> ExtractFeaturesAll(
    std::span<const std::string> contig_seqs,
    const FeatureExtractorConfig& cfg = {});

// --- low-level primitives, exposed for unit tests ---

float ShannonEntropyKmer(std::string_view window, uint32_t k);
float GcContent(std::string_view window);

// Detects strongest tandem period in [min_p, max_p]. Returns -1 if no
// significant signal. Uses autocorrelation against the window.
int32_t DetectTandemPeriod(std::string_view window, uint32_t min_p, uint32_t max_p);

// Density of palindromic k-mers (matching reverse-complement of themselves)
// over the window for a fixed small k (default 6). Useful for ITR/inverted
// repeat detection.
float PalindromeDensity(std::string_view window, uint32_t k = 6);

}  // namespace llmap::annot
