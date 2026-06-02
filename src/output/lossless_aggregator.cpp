// LLmap — Lossless aggregator implementation.
//
// Counter updates run under a single mutex; for the per-thread call
// volume we expect (one Observe() per emitted record, batched ~50k at
// a time) the contention is negligible. If it ever shows up in profiles
// we can switch to per-thread shards with periodic merge.

#include "output/lossless_aggregator.h"

#include "core/transcript_kind.h"

#include <fstream>
#include <sstream>

namespace llmap::output {

namespace {

// JSON-escape a string for emitting into the summary file. Minimal
// implementation — only handles the subset of characters that can show
// up in read IDs and enum names. Read IDs are guaranteed by FASTQ spec
// to be ASCII printable without quotes or backslashes, so we just emit
// them verbatim with a defensive replacement of `"` and `\`.
std::string JsonEscape(std::string_view s) {
    std::string out;
    out.reserve(s.size() + 4);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

}  // namespace

void LosslessAggregator::SetExpectedInputCount(std::uint64_t n) noexcept {
    std::lock_guard<std::mutex> g(mu_);
    c_.n_input_reads_declared = n;
}

bool LosslessAggregator::Observe(const AlignmentRecord& rec) {
    const bool ok = rec.is_lossless_consistent();

    std::lock_guard<std::mutex> g(mu_);
    ++c_.n_records_emitted;

    // Per-status (status enum guaranteed by IsAnyMapped switch to be 0..9)
    const auto status_idx = static_cast<std::size_t>(rec.status);
    if (status_idx < c_.by_status.size()) {
        ++c_.by_status[status_idx];
    }

    // Per-rejection-reason (only set when rejection_reason has a value)
    if (rec.rejection_reason.has_value()) {
        const auto rj_idx = static_cast<std::size_t>(*rec.rejection_reason);
        if (rj_idx < c_.by_rejection.size()) {
            ++c_.by_rejection[rj_idx];
        }
    }

    // Per-TranscriptKind
    ++c_.by_kind[static_cast<std::uint16_t>(rec.transcript_kind)];

    // Invariant tracking
    if (!ok) {
        c_.invariant_ok = false;
        if (c_.offending_read_ids.size() < kMaxOffenders) {
            c_.offending_read_ids.push_back(rec.read_id);
        }
    }

    return ok;
}

LosslessCounters LosslessAggregator::Snapshot() const {
    std::lock_guard<std::mutex> g(mu_);
    return c_;
}

bool LosslessAggregator::WriteSummary(
    const std::filesystem::path& path) const {

    LosslessCounters snap = Snapshot();

    // Cross-check: n_records vs declared input count
    const bool n_matches = (snap.n_input_reads_declared == 0)  // unknown ⇒ skip
        || snap.n_input_reads_declared == snap.n_records_emitted;

    std::ofstream out(path);
    if (!out) return false;

    out << "{\n";
    out << "  \"n_input_reads_declared\": " << snap.n_input_reads_declared << ",\n";
    out << "  \"n_records_emitted\":     " << snap.n_records_emitted << ",\n";
    out << "  \"counts_match\":           " << (n_matches ? "true" : "false") << ",\n";
    out << "  \"lossless_invariant_ok\": " << (snap.invariant_ok ? "true" : "false") << ",\n";

    // by_status block
    out << "  \"by_status\": {\n";
    bool first = true;
    for (std::size_t i = 0; i < snap.by_status.size(); ++i) {
        if (snap.by_status[i] == 0) continue;
        if (!first) out << ",\n";
        first = false;
        out << "    \"" << AlignmentStatusName(static_cast<AlignmentStatus>(i))
            << "\": " << snap.by_status[i];
    }
    out << "\n  },\n";

    // by_rejection block
    out << "  \"by_rejection\": {\n";
    first = true;
    for (std::size_t i = 0; i < snap.by_rejection.size(); ++i) {
        if (snap.by_rejection[i] == 0) continue;
        if (!first) out << ",\n";
        first = false;
        out << "    \"" << RejectionReasonName(static_cast<RejectionReason>(i))
            << "\": " << snap.by_rejection[i];
    }
    out << "\n  },\n";

    // by_kind block
    out << "  \"by_kind\": {\n";
    first = true;
    for (const auto& [k, n] : snap.by_kind) {
        if (n == 0) continue;
        if (!first) out << ",\n";
        first = false;
        out << "    \""
            << core::TranscriptKindName(static_cast<core::TranscriptKind>(k))
            << "\": " << n;
    }
    out << "\n  }";

    // Offenders (only when invariant broke)
    if (!snap.invariant_ok) {
        out << ",\n  \"offending_read_ids\": [\n";
        for (std::size_t i = 0; i < snap.offending_read_ids.size(); ++i) {
            out << "    \"" << JsonEscape(snap.offending_read_ids[i]) << "\"";
            if (i + 1 < snap.offending_read_ids.size()) out << ",";
            out << "\n";
        }
        out << "  ]";
    }

    out << "\n}\n";
    return static_cast<bool>(out);
}

void LosslessAggregator::Reset() {
    std::lock_guard<std::mutex> g(mu_);
    c_ = LosslessCounters{};
}

}  // namespace llmap::output
