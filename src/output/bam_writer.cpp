// LLmap — BAM/SAM output writer implementation.

#include "output/bam_writer.h"

#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace llmap::output {

namespace {

// SAM flags
constexpr std::uint16_t SAM_UNMAPPED = 0x4;
constexpr std::uint16_t SAM_SECONDARY = 0x100;
constexpr std::uint16_t SAM_SUPPLEMENTARY = 0x800;

// Map rejection reason to string for tags
std::string RejectionReasonTag(RejectionReason reason) {
    switch (reason) {
        case RejectionReason::NoSeeds: return "NO_SEEDS";
        case RejectionReason::LowSeedDensity: return "LOW_SEED_DENSITY";
        case RejectionReason::ChainScoreBelowThreshold: return "LOW_CHAIN_SCORE";
        case RejectionReason::AmbiguousNoAnchor: return "AMBIGUOUS";
        case RejectionReason::FailedWfa2Extension: return "WFA2_FAILED";
        case RejectionReason::DidNotConverge: return "NO_CONVERGE";
        case RejectionReason::HostContamination: return "HOST_CONTAM";
        case RejectionReason::LowComplexity: return "LOW_COMPLEXITY";
    }
    return "UNKNOWN";
}

// Calculate MAPQ from confidence and alternatives
std::uint8_t CalculateMapq(const AlignmentRecord& rec) {
    if (rec.status == AlignmentStatus::Unmapped) return 0;
    if (rec.status == AlignmentStatus::Tentative) return 0;

    // For mapped: higher confidence = higher MAPQ
    if (!rec.confidence_scores.empty()) {
        float conf = rec.confidence_scores[0];
        if (conf >= 0.99f) return 60;
        if (conf >= 0.95f) return 40;
        if (conf >= 0.90f) return 30;
        if (conf >= 0.80f) return 20;
        return 10;
    }

    // Fallback: penalize based on number of alternatives
    if (rec.alternatives.empty()) return 60;
    if (rec.alternatives.size() == 1) return 30;
    if (rec.alternatives.size() <= 3) return 20;
    return 10;
}

}  // namespace

// ========== Implementation class ==========

class BamWriterImpl {
public:
    std::ofstream file;
    std::string last_error;
    BamWriterStats stats;
    std::unordered_map<std::string, std::size_t> ref_name_to_idx;
    bool header_written = false;
};

// ========== BamWriter ==========

BamWriter::BamWriter(const std::filesystem::path& path,
                     std::span<const ReferenceSequence> references,
                     const BamWriterConfig& config)
    : path_(path), config_(config),
      references_(references.begin(), references.end()),
      impl_(std::make_unique<BamWriterImpl>()) {}

BamWriter::~BamWriter() {
    Close();
}

BamWriter::BamWriter(BamWriter&&) noexcept = default;
BamWriter& BamWriter::operator=(BamWriter&&) noexcept = default;

std::unique_ptr<BamWriter> BamWriter::Create(
    const std::filesystem::path& path,
    std::span<const ReferenceSequence> references,
    const BamWriterConfig& config) {

    // BAM requires htslib
    if (config.format == BamOutputFormat::BAM && !HtslibAvailable()) {
        return nullptr;
    }

    auto writer = std::unique_ptr<BamWriter>(
        new BamWriter(path, references, config));
    if (!writer->Initialize()) {
        return nullptr;
    }
    return writer;
}

bool BamWriter::Initialize() {
    impl_->file.open(path_, std::ios::out | std::ios::trunc);
    if (!impl_->file) {
        impl_->last_error = "Failed to open output file: " + path_.string();
        return false;
    }

    // Build reference name lookup
    for (std::size_t i = 0; i < references_.size(); ++i) {
        impl_->ref_name_to_idx[references_[i].name] = i;
    }

    // Write SAM header
    if (config_.write_header) {
        impl_->file << "@HD\tVN:1.6\tSO:unsorted\n";

        // Reference sequences
        for (const auto& ref : references_) {
            impl_->file << "@SQ\tSN:" << ref.name
                        << "\tLN:" << ref.length << "\n";
        }

        // Read group
        impl_->file << "@RG\tID:" << config_.read_group_id
                    << "\tSM:" << config_.sample_name
                    << "\tPL:" << config_.platform << "\n";

        // Program
        impl_->file << "@PG\tID:" << config_.program_name
                    << "\tPN:" << config_.program_name
                    << "\tVN:" << config_.program_version
                    << "\tCL:llmap align\n";

        impl_->header_written = true;
    }

    return true;
}

bool BamWriter::Write(const AlignmentRecord& record,
                      std::string_view query_seq) {
    if (!impl_->file.is_open()) {
        impl_->last_error = "Writer not initialized or already closed";
        return false;
    }

    auto start = std::chrono::steady_clock::now();

    // Handle based on status
    if (record.status == AlignmentStatus::Unmapped) {
        if (!config_.write_unmapped) {
            return true;  // Skip, not an error
        }

        // Unmapped read
        impl_->file << record.read_id << "\t"
                    << SAM_UNMAPPED << "\t"
                    << "*\t0\t0\t*\t*\t0\t0\t";

        if (config_.include_sequences && !query_seq.empty()) {
            impl_->file << query_seq;
        } else {
            impl_->file << "*";
        }
        impl_->file << "\t*";

        // Tags for rejection reason
        if (record.rejection_reason) {
            impl_->file << "\tXR:Z:" << RejectionReasonTag(*record.rejection_reason);
        }

        // WaveCollapse tags
        if (config_.include_wavecollapse_tags) {
            impl_->file << "\tXC:i:" << record.cluster_id
                        << "\tXL:i:" << static_cast<int>(record.collapsed_at_level)
                        << "\tXI:i:" << record.collapsed_at_iteration;
        }

        impl_->file << "\n";
        impl_->stats.unmapped_written++;

    } else if (record.status == AlignmentStatus::Tentative) {
        if (!config_.include_tentative) {
            return true;
        }

        // Write as unmapped with tentative targets in tags
        impl_->file << record.read_id << "\t"
                    << SAM_UNMAPPED << "\t"
                    << "*\t0\t0\t*\t*\t0\t0\t";

        if (config_.include_sequences && !query_seq.empty()) {
            impl_->file << query_seq;
        } else {
            impl_->file << "*";
        }
        impl_->file << "\t*";

        // XR tag with reason
        if (record.rejection_reason) {
            impl_->file << "\tXR:Z:" << RejectionReasonTag(*record.rejection_reason);
        }

        // XT tag: tentative targets
        if (!record.tentative_targets.empty()) {
            impl_->file << "\tXT:Z:";
            for (std::size_t i = 0; i < record.tentative_targets.size(); ++i) {
                const auto& t = record.tentative_targets[i];
                if (i > 0) impl_->file << ";";
                impl_->file << t.target_id << ":"
                            << t.approx_start << "-" << t.approx_end
                            << "(" << t.final_probability << ")";
            }
        }

        // WaveCollapse tags
        if (config_.include_wavecollapse_tags) {
            impl_->file << "\tXC:i:" << record.cluster_id
                        << "\tXL:i:" << static_cast<int>(record.collapsed_at_level)
                        << "\tXI:i:" << record.collapsed_at_iteration;
        }

        impl_->file << "\n";
        impl_->stats.tentative_written++;

    } else {
        // Mapped - has primary alignment
        if (!record.primary) {
            impl_->last_error = "Mapped record without primary alignment";
            return false;
        }

        const auto& hit = *record.primary;
        std::uint8_t mapq = CalculateMapq(record);

        // Get reference index
        auto ref_it = impl_->ref_name_to_idx.find(hit.target_id);
        std::string rname = (ref_it != impl_->ref_name_to_idx.end())
            ? hit.target_id : "*";

        // CIGAR
        std::string cigar_str = hit.cigar.ops.empty() ? "*" : hit.cigar.ops;

        // Primary alignment
        impl_->file << record.read_id << "\t"
                    << 0 << "\t"  // No flags for primary mapped
                    << rname << "\t"
                    << (hit.start + 1) << "\t"  // 1-based
                    << static_cast<int>(mapq) << "\t"
                    << cigar_str << "\t"
                    << "*\t0\t0\t";  // RNEXT, PNEXT, TLEN

        if (config_.include_sequences && !query_seq.empty()) {
            impl_->file << query_seq;
        } else {
            impl_->file << "*";
        }
        impl_->file << "\t*";  // QUAL

        // AS tag (alignment score)
        impl_->file << "\tAS:i:" << hit.score;

        // NM tag (edit distance)
        impl_->file << "\tNM:i:" << hit.nm;

        // RG tag
        impl_->file << "\tRG:Z:" << config_.read_group_id;

        // XA tag for alternatives
        if (config_.include_alternatives && !record.alternatives.empty()) {
            impl_->file << "\tXA:Z:";
            for (std::size_t i = 0; i < record.alternatives.size(); ++i) {
                const auto& alt = record.alternatives[i];
                if (i > 0) impl_->file << ";";
                impl_->file << alt.target_id << ","
                            << "+" << (alt.start + 1) << ","
                            << (alt.cigar.ops.empty() ? "*" : alt.cigar.ops) << ","
                            << alt.nm;
            }
        }

        // WaveCollapse tags
        if (config_.include_wavecollapse_tags) {
            impl_->file << "\tXC:i:" << record.cluster_id
                        << "\tXL:i:" << static_cast<int>(record.collapsed_at_level)
                        << "\tXI:i:" << record.collapsed_at_iteration;
            if (record.is_cluster_representative) {
                impl_->file << "\tXP:i:1";  // Representative flag
            }
        }

        // Paralog tags
        if (config_.include_paralog_tags && record.paralog_assignment) {
            const auto& pa = *record.paralog_assignment;
            impl_->file << "\tPD:i:" << pa.n_discriminating_psvs;
            if (pa.p_canonical) {
                impl_->file << "\tPC:f:" << *pa.p_canonical;
            }
            if (!pa.inter_paralog.empty()) {
                impl_->file << "\tPP:Z:";
                for (std::size_t i = 0; i < pa.inter_paralog.size(); ++i) {
                    const auto& [name, prob] = pa.inter_paralog[i];
                    if (i > 0) impl_->file << ",";
                    impl_->file << name << ":" << prob;
                }
            }
        }

        impl_->file << "\n";
        impl_->stats.mapped_written++;
        impl_->stats.alternatives_written += record.alternatives.size();
    }

    impl_->stats.records_written++;

    auto end = std::chrono::steady_clock::now();
    impl_->stats.write_time_ms +=
        std::chrono::duration<float, std::milli>(end - start).count();

    return true;
}

bool BamWriter::WriteBatch(std::span<const AlignmentRecord> records) {
    for (const auto& record : records) {
        if (!Write(record)) {
            return false;
        }
    }
    return true;
}

bool BamWriter::Close() {
    if (impl_->file.is_open()) {
        impl_->file.close();
        return true;
    }
    return false;
}

BamWriterStats BamWriter::GetStats() const {
    return impl_->stats;
}

std::string BamWriter::LastError() const {
    return impl_->last_error;
}

bool BamWriter::HtslibAvailable() {
#ifdef LLMAP_HAS_HTSLIB
    return true;
#else
    return false;
#endif
}

// ========== Convenience function ==========

bool WriteAlignments(
    const std::filesystem::path& path,
    std::span<const AlignmentRecord> records,
    std::span<const ReferenceSequence> references,
    const BamWriterConfig& config) {

    auto writer = BamWriter::Create(path, references, config);
    if (!writer) {
        return false;
    }

    if (!writer->WriteBatch(records)) {
        return false;
    }

    return writer->Close();
}

}  // namespace llmap::output
