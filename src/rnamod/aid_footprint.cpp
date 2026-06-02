// LLmap — AidFootprintDetector implementation.
//
// Detection rule (Methot & Di Noia 2017, Casellas et al 2016):
//   * Walk the read/anchor in pairs (assumed pre-aligned, equal length).
//   * Count positions where anchor=C and read=T (the AID signature).
//   * For each such position, check whether the LEADING dinucleotide
//     on the anchor matches the AID WRC hotspot (W={A,T}, R={A,G}, C).
//   * Require ≥3 events AND ≥50 % in-hotspot for a positive call.
//
// We trust the caller to feed us aligned anchor + read of identical
// length. If lengths differ we return detected=false rather than
// guessing at the alignment; the upstream pipeline should resolve that.

#include "rnamod/aid_footprint.h"

#include <cctype>

namespace llmap::rnamod {

namespace {

inline char Upper(char c) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

bool IsW(char c) { return c == 'A' || c == 'T' || c == 'U'; }
bool IsR(char c) { return c == 'A' || c == 'G'; }

}  // namespace

AidFootprintResult AidFootprintDetector::Detect(
    std::string_view read_seq,
    std::string_view anchor_seq,
    std::string_view anchor_switch_id) const {

    AidFootprintResult r;
    if (anchor_switch_id.empty()) {
        // Not a switch-region anchor; nothing to detect.
        return r;
    }
    if (read_seq.size() != anchor_seq.size()) {
        // Length mismatch — caller should have aligned read to anchor first.
        return r;
    }

    std::uint32_t hotspot_hits = 0;

    // Walk positions; we need leading 2 bp of context for WRC hotspot
    // check, so start at i=2.
    for (std::size_t i = 2; i < anchor_seq.size(); ++i) {
        const char a = Upper(anchor_seq[i]);
        const char b = Upper(read_seq[i]);
        if (a != 'C') continue;
        if (b != 'T' && b != 'U') continue;

        // C→U / C→T event. Record position.
        ++r.n_c_to_u_events;
        r.c_to_u_positions_in_read.push_back(
            static_cast<std::uint32_t>(i));

        // Hotspot check: previous-two-bp WR + this-bp C.
        const char w = Upper(anchor_seq[i - 2]);
        const char rR = Upper(anchor_seq[i - 1]);
        if (IsW(w) && IsR(rR)) {
            ++hotspot_hits;
        }
    }

    r.switch_region_id = std::string(anchor_switch_id);

    if (r.n_c_to_u_events >= 3) {
        const float hotspot_frac =
            static_cast<float>(hotspot_hits)
            / static_cast<float>(r.n_c_to_u_events);
        if (hotspot_frac >= 0.5f) {
            r.detected = true;
            r.confidence = hotspot_frac;
        }
    }
    return r;
}

}  // namespace llmap::rnamod
