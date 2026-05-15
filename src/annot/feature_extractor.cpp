#include "annot/feature_extractor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace llmap::annot {

namespace {

inline int BaseIdx(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

inline char Complement(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

}  // namespace

float ShannonEntropyKmer(std::string_view window, uint32_t k) {
    if (window.size() < k || k == 0 || k > 8) return 0.0f;

    // Count canonical k-mer occurrences (encoded as base-4 integers).
    // We use a flat array up to k=8 (65536 entries) for cache efficiency.
    std::vector<uint32_t> counts(1u << (2 * k), 0);
    uint32_t total = 0;
    uint64_t encoded = 0;
    int valid = 0;
    uint64_t mask = (1ull << (2 * k)) - 1;

    for (size_t i = 0; i < window.size(); ++i) {
        int b = BaseIdx(window[i]);
        if (b < 0) {
            encoded = 0;
            valid = 0;
            continue;
        }
        encoded = ((encoded << 2) | static_cast<uint64_t>(b)) & mask;
        ++valid;
        if (valid >= static_cast<int>(k)) {
            ++counts[encoded];
            ++total;
        }
    }

    if (total == 0) return 0.0f;

    float entropy = 0.0f;
    float inv_total = 1.0f / static_cast<float>(total);
    for (uint32_t c : counts) {
        if (c == 0) continue;
        float p = c * inv_total;
        entropy -= p * std::log2(p);
    }
    return entropy;
}

float GcContent(std::string_view window) {
    if (window.empty()) return 0.0f;
    uint32_t gc = 0, valid = 0;
    for (char c : window) {
        int b = BaseIdx(c);
        if (b < 0) continue;
        ++valid;
        if (b == 1 || b == 2) ++gc;  // C or G
    }
    return valid > 0 ? static_cast<float>(gc) / static_cast<float>(valid) : 0.0f;
}

int32_t DetectTandemPeriod(std::string_view window, uint32_t min_p, uint32_t max_p) {
    if (window.size() < 2 * min_p) return -1;
    uint32_t cap = std::min<uint32_t>(max_p, static_cast<uint32_t>(window.size() / 2));
    if (cap < min_p) return -1;

    // Cheap autocorrelation: for each candidate period p, count matching
    // bases between window[0..N-p] and window[p..N]. Pick the p with the
    // highest match fraction provided it crosses a significance threshold.
    int32_t best_period = -1;
    float best_frac = 0.0f;
    const float kSignificance = 0.85f;

    for (uint32_t p = min_p; p <= cap; ++p) {
        uint32_t matches = 0;
        uint32_t compared = 0;
        for (size_t i = 0; i + p < window.size(); ++i) {
            int b1 = BaseIdx(window[i]);
            int b2 = BaseIdx(window[i + p]);
            if (b1 < 0 || b2 < 0) continue;
            ++compared;
            if (b1 == b2) ++matches;
        }
        if (compared == 0) continue;
        float frac = static_cast<float>(matches) / static_cast<float>(compared);
        if (frac > best_frac) {
            best_frac = frac;
            best_period = static_cast<int32_t>(p);
        }
    }

    return best_frac >= kSignificance ? best_period : -1;
}

float PalindromeDensity(std::string_view window, uint32_t k) {
    if (window.size() < k || k == 0 || k > 8) return 0.0f;

    uint32_t palindromes = 0;
    uint32_t total_kmers = 0;
    for (size_t i = 0; i + k <= window.size(); ++i) {
        bool valid = true;
        bool is_palindrome = true;
        for (uint32_t j = 0; j < k; ++j) {
            int b1 = BaseIdx(window[i + j]);
            int b2 = BaseIdx(window[i + k - 1 - j]);
            if (b1 < 0 || b2 < 0) {
                valid = false;
                break;
            }
            // Palindrome means window[i+j] is complement of window[i+k-1-j]
            int complement_idx = 3 - b2;  // A<->T (0<->3), C<->G (1<->2)
            if (b1 != complement_idx) {
                is_palindrome = false;
                break;
            }
        }
        if (!valid) continue;
        ++total_kmers;
        if (is_palindrome) ++palindromes;
    }
    return total_kmers > 0
        ? static_cast<float>(palindromes) / static_cast<float>(total_kmers)
        : 0.0f;
}

std::vector<WindowFeatures> ExtractFeatures(
    uint32_t ref_id,
    std::string_view seq,
    const FeatureExtractorConfig& cfg) {

    std::vector<WindowFeatures> out;
    if (seq.empty() || cfg.window_bp == 0) return out;

    out.reserve(seq.size() / cfg.step_bp + 1);
    for (uint32_t start = 0; start < seq.size(); start += cfg.step_bp) {
        uint32_t end = std::min<uint32_t>(start + cfg.window_bp,
                                          static_cast<uint32_t>(seq.size()));
        auto window = seq.substr(start, end - start);
        WindowFeatures f;
        f.ref_id = ref_id;
        f.start = start;
        f.end = end;
        f.shannon_5mer = ShannonEntropyKmer(window, cfg.entropy_kmer_k);
        f.gc_content = GcContent(window);
        f.tandem_period = DetectTandemPeriod(window, cfg.tandem_min_period,
                                             cfg.tandem_max_period);
        f.palindrome_density = PalindromeDensity(window, 6);
        // kmer_multiplicity_p95, orf_density, consensus_match require external
        // input (full-reference index, codon table, motif library) and are
        // filled in by ExtractFeaturesAll / the classifier stage.
        out.push_back(std::move(f));
    }
    return out;
}

std::vector<WindowFeatures> ExtractFeaturesAll(
    std::span<const std::string> contig_seqs,
    const FeatureExtractorConfig& cfg) {

    std::vector<WindowFeatures> out;
    for (size_t i = 0; i < contig_seqs.size(); ++i) {
        auto feats = ExtractFeatures(static_cast<uint32_t>(i),
                                     contig_seqs[i], cfg);
        out.insert(out.end(), feats.begin(), feats.end());
    }
    return out;
}

}  // namespace llmap::annot
