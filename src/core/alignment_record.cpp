#include "core/alignment_record.h"

namespace llmap {

bool AlignmentRecord::is_lossless_consistent() const noexcept {
    if (read_id.empty()) return false;
    switch (status) {
        case AlignmentStatus::Mapped:
            return primary.has_value() && !rejection_reason.has_value();
        case AlignmentStatus::Tentative:
            return !primary.has_value() && !tentative_targets.empty() && rejection_reason.has_value();
        case AlignmentStatus::Unmapped:
            return !primary.has_value() && rejection_reason.has_value();
    }
    return false;
}

AlignmentRecord make_mapped(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    std::vector<AlignmentHit> alternatives
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::Mapped;
    r.primary = std::move(primary);
    r.alternatives = std::move(alternatives);
    return r;
}

AlignmentRecord make_tentative(
    std::string read_id,
    std::uint32_t read_len,
    std::vector<TentativeTarget> targets,
    RejectionReason reason
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::Tentative;
    r.tentative_targets = std::move(targets);
    r.rejection_reason = reason;
    return r;
}

AlignmentRecord make_unmapped(
    std::string read_id,
    std::uint32_t read_len,
    RejectionReason reason
) {
    AlignmentRecord r;
    r.read_id = std::move(read_id);
    r.read_len = read_len;
    r.status = AlignmentStatus::Unmapped;
    r.rejection_reason = reason;
    return r;
}

}  // namespace llmap
