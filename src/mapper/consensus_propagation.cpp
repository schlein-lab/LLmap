// LLmap — consensus-contig placement propagation implementation.

#include "mapper/consensus_propagation.h"

namespace llmap::mapper {

std::vector<MemberPlacement> PropagatePlacement(
    const ConsensusContig& contig,
    std::span<const std::string> reads,
    const ProbePlacement& probe) {
    std::vector<MemberPlacement> out;
    out.reserve(contig.members.size());

    const std::int64_t L = static_cast<std::int64_t>(contig.length);
    for (const auto& m : contig.members) {
        MemberPlacement mp;
        mp.read_idx = m.read_idx;
        mp.ref_name = probe.ref_name;
        mp.anchored = m.anchored;
        if (!m.anchored) {
            // Caller maps this read solo — no contig-derived coordinate.
            out.push_back(mp);
            continue;
        }

        const std::int64_t read_len =
            (m.read_idx < reads.size())
                ? static_cast<std::int64_t>(reads[m.read_idx].size())
                : 0;
        const std::int64_t o = m.offset;

        if (!probe.reverse) {
            // contig column c → ref (ref_start + c).
            mp.ref_start = probe.ref_start + o;
            mp.reverse = m.reverse;
        } else {
            // contig reverse-complemented onto the ref: column c → ref_start+L-1-c.
            // member [o, o+read_len) → ref_start = ref_start + L - (o + read_len).
            mp.ref_start = probe.ref_start + L - (o + read_len);
            mp.reverse = !m.reverse;
        }
        out.push_back(mp);
    }
    return out;
}

}  // namespace llmap::mapper
