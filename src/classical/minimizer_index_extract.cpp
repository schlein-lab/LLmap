#include "classical/minimizer_index.h"
#include "classical/minimizer_index_impl.h"

#include <deque>

namespace llmap::classical {

using detail::EncodeBase;
using detail::IsValidBase;
using detail::Complement;

// Invertible hash function (from minimap2)
uint64_t HashKmer(uint64_t kmer, uint64_t seed) {
    kmer = (~kmer + (kmer << 21)) ^ seed;
    kmer = kmer ^ (kmer >> 24);
    kmer = (kmer + (kmer << 3)) + (kmer << 8);
    kmer = kmer ^ (kmer >> 14);
    kmer = (kmer + (kmer << 2)) + (kmer << 4);
    kmer = kmer ^ (kmer >> 28);
    kmer = kmer + (kmer << 31);
    return kmer;
}

uint64_t ReverseComplementKmer(uint64_t kmer, uint8_t k) {
    uint64_t rc = 0;
    for (uint8_t i = 0; i < k; ++i) {
        rc = (rc << 2) | Complement(kmer & 3);
        kmer >>= 2;
    }
    return rc;
}

uint64_t CanonicalKmer(uint64_t kmer, uint8_t k) {
    uint64_t rc = ReverseComplementKmer(kmer, k);
    return std::min(kmer, rc);
}

// Standalone minimizer extraction
std::vector<Minimizer> ExtractMinimizers(
    std::string_view sequence,
    const MinimizerConfig& config)
{
    std::vector<Minimizer> result;
    const size_t k = config.k;
    const size_t w = config.w;

    if (sequence.size() < k) return result;

    const size_t kmer_mask = (1ULL << (2 * k)) - 1;

    // Sliding window of (hash, position, is_reverse) for minimizer selection
    struct WindowEntry {
        uint64_t hash;
        uint32_t pos;
        bool is_reverse;
    };
    std::deque<WindowEntry> window;

    // Build initial k-mer
    uint64_t fwd_kmer = 0;
    uint64_t rev_kmer = 0;
    size_t valid_bases = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        char c = sequence[i];

        if (!IsValidBase(c)) {
            valid_bases = 0;
            fwd_kmer = 0;
            rev_kmer = 0;
            continue;
        }

        uint8_t base = EncodeBase(c);
        fwd_kmer = ((fwd_kmer << 2) | base) & kmer_mask;
        rev_kmer = (rev_kmer >> 2) | (static_cast<uint64_t>(Complement(base)) << (2 * (k - 1)));
        ++valid_bases;

        if (valid_bases < k) continue;

        // Position of this k-mer
        uint32_t pos = static_cast<uint32_t>(i - k + 1);

        // Choose canonical k-mer
        uint64_t kmer_to_hash;
        bool is_reverse;
        if (config.canonical) {
            if (fwd_kmer <= rev_kmer) {
                kmer_to_hash = fwd_kmer;
                is_reverse = false;
            } else {
                kmer_to_hash = rev_kmer;
                is_reverse = true;
            }
        } else {
            kmer_to_hash = fwd_kmer;
            is_reverse = false;
        }

        uint64_t hash = HashKmer(kmer_to_hash, config.seed);

        // Add to window
        WindowEntry entry{hash, pos, is_reverse};

        // Remove entries that have fallen out of window
        while (!window.empty() && window.front().pos + w <= pos) {
            window.pop_front();
        }

        // Remove entries from back that are larger than new entry (monotonic deque)
        while (!window.empty() && window.back().hash > hash) {
            window.pop_back();
        }

        window.push_back(entry);

        // Output minimizer once we have a full window
        if (pos >= w - 1 || i == sequence.size() - 1) {
            // The front of the deque is the minimum in the current window
            const auto& min_entry = window.front();

            // Only output if this is the first time we're seeing this minimizer
            // (avoids duplicates for minimizers that span multiple windows)
            if (result.empty() || result.back().pos != min_entry.pos) {
                Minimizer m;
                m.hash = min_entry.hash;
                m.pos = min_entry.pos;
                m.is_reverse = min_entry.is_reverse;
                result.push_back(m);
            }
        }
    }

    return result;
}

}  // namespace llmap::classical
