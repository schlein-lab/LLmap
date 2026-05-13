#include "classical/minimizer_index.h"
#include "classical/minimizer_index_impl.h"
#include "core/arena.h"

#include <deque>

namespace llmap::classical {

namespace {

// Sliding window entry for monotonic deque
struct WindowEntry {
    uint64_t hash;
    uint32_t pos;
    bool is_reverse;
};

// Core extraction logic, templated on output type
template <typename OutputContainer>
void ExtractMinimizersImpl(
    std::string_view sequence,
    const MinimizerConfig& config,
    OutputContainer& result)
{
    using detail::EncodeBase;
    using detail::IsValidBase;
    using detail::Complement;

    const size_t k = config.k;
    const size_t w = config.w;

    if (sequence.size() < k) return;

    const size_t kmer_mask = (1ULL << (2 * k)) - 1;

    // Use ScratchBuffer for the deque to avoid per-read allocation
    thread_local core::ScratchBuffer<WindowEntry> window_buf;
    window_buf.clear();

    size_t window_front = 0;  // Index of front in window_buf

    uint64_t fwd_kmer = 0;
    uint64_t rev_kmer = 0;
    size_t valid_bases = 0;

    uint32_t last_output_pos = UINT32_MAX;

    for (size_t i = 0; i < sequence.size(); ++i) {
        char c = sequence[i];

        if (!IsValidBase(c)) {
            valid_bases = 0;
            fwd_kmer = 0;
            rev_kmer = 0;
            window_buf.clear();
            window_front = 0;
            continue;
        }

        uint8_t base = EncodeBase(c);
        fwd_kmer = ((fwd_kmer << 2) | base) & kmer_mask;
        rev_kmer = (rev_kmer >> 2) | (static_cast<uint64_t>(Complement(base)) << (2 * (k - 1)));
        ++valid_bases;

        if (valid_bases < k) continue;

        uint32_t pos = static_cast<uint32_t>(i - k + 1);

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

        // Remove entries that have fallen out of window
        while (window_front < window_buf.size() &&
               window_buf[window_front].pos + w <= pos) {
            ++window_front;
        }

        // Remove entries from back that are larger (monotonic property)
        while (window_buf.size() > window_front &&
               window_buf.back().hash > hash) {
            window_buf.pop_back();
        }

        window_buf.push_back({hash, pos, is_reverse});

        // Output minimizer once we have a full window
        if (pos >= w - 1 || i == sequence.size() - 1) {
            const auto& min_entry = window_buf[window_front];

            if (min_entry.pos != last_output_pos) {
                Minimizer m;
                m.hash = min_entry.hash;
                m.pos = min_entry.pos;
                m.is_reverse = min_entry.is_reverse;
                result.push_back(m);
                last_output_pos = min_entry.pos;
            }
        }
    }
}

}  // namespace

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

// Standalone minimizer extraction (allocates result vector)
std::vector<Minimizer> ExtractMinimizers(
    std::string_view sequence,
    const MinimizerConfig& config)
{
    std::vector<Minimizer> result;
    ExtractMinimizersImpl(sequence, config, result);
    return result;
}

// Zero-allocation version using scratch buffer
void ExtractMinimizersInto(
    std::string_view sequence,
    const MinimizerConfig& config,
    core::ScratchBuffer<Minimizer>& out)
{
    out.clear();
    ExtractMinimizersImpl(sequence, config, out);
}

}  // namespace llmap::classical
