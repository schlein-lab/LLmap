// LLmap — Core validation logic for real reference validation.
//
// Runs alignment pipeline and compares against baseline/ground truth.

#include "validation/real_reference.h"
#include "validation/real_reference_internal.h"

#include <sstream>

#include "classical/classical_pipeline.h"
#include "io/fasta_reader.h"
#include "io/fastq_reader.h"

namespace llmap::validation {

bool RealReferenceConfig::Validate() const {
    return ValidationErrors().empty();
}

std::string RealReferenceConfig::ValidationErrors() const {
    std::ostringstream errors;

    if (reference_fasta.empty()) {
        errors << "reference_fasta is required\n";
    } else if (!std::filesystem::exists(reference_fasta)) {
        errors << "reference_fasta not found: " << reference_fasta << "\n";
    }

    if (reads_fastq.empty()) {
        errors << "reads_fastq is required\n";
    } else if (!std::filesystem::exists(reads_fastq)) {
        errors << "reads_fastq not found: " << reads_fastq << "\n";
    }

    if (!minimap2_bam.empty() && !std::filesystem::exists(minimap2_bam)) {
        errors << "minimap2_bam not found: " << minimap2_bam << "\n";
    }

    if (!ground_truth_bed.empty() && !std::filesystem::exists(ground_truth_bed)) {
        errors << "ground_truth_bed not found: " << ground_truth_bed << "\n";
    }

    return errors.str();
}

std::string RealReferenceResult::Summary() const {
    std::ostringstream ss;
    ss << "=== Real Reference Validation ===\n\n";

    ss << "Reference:\n";
    ss << "  Sequences:    " << n_reference_seqs << "\n";
    ss << "  Total bp:     " << reference_total_bp << "\n\n";

    ss << "Reads:\n";
    ss << "  Total:        " << n_reads << "\n\n";

    ss << "LLmap Results:\n";
    ss << "  Mapped:       " << llmap_mapped << "\n";
    ss << "  Unmapped:     " << llmap_unmapped << "\n";
    ss << "  Rate:         " << (llmap_alignment_rate * 100.0f) << "%\n\n";

    if (has_minimap2_baseline) {
        ss << "Minimap2 Baseline:\n";
        ss << "  Mapped:       " << minimap2_mapped << "\n";
        ss << "  Both mapped:  " << baseline_comparison.both_mapped << "\n";
        ss << "  LLmap only:   " << baseline_comparison.llmap_only << "\n";
        ss << "  MM2 only:     " << baseline_comparison.minimap2_only << "\n";
        ss << "  Recall ratio: " << (baseline_comparison.recall_ratio * 100.0)
           << "%\n\n";
    }

    if (has_ground_truth) {
        ss << "Position Accuracy:\n";
        ss << "  Evaluated:    " << validation.position_evaluated << "\n";
        ss << "  Correct:      " << validation.position_correct << "\n";
        ss << "  Accuracy:     " << (validation.position_accuracy * 100.0)
           << "%\n\n";
    }

    ss << "Timing:\n";
    ss << "  Index load:   " << index_load_time_ms << " ms\n";
    ss << "  Alignment:    " << alignment_time_ms << " ms\n";
    ss << "  Validation:   " << validation_time_ms << " ms\n";
    ss << "  Total:        " << total_time_ms << " ms\n\n";

    ss << "Kill-Switch Verdicts:\n";
    ss << "  Recall (≥99.5% of MM2):  " << (recall_pass ? "PASS" : "FAIL") << "\n";
    ss << "  Lossless (no drops):     " << (lossless_pass ? "PASS" : "FAIL") << "\n";
    ss << "  Overall:                 " << (overall_pass ? "PASS" : "FAIL");
    if (!overall_pass && !verdict_reason.empty()) {
        ss << " (" << verdict_reason << ")";
    }
    ss << "\n";

    return ss.str();
}

RealReferenceResult RunRealReferenceValidation(const RealReferenceConfig& config) {
    internal::Timer total_timer;
    RealReferenceResult result;

    if (!config.Validate()) {
        result.verdict_reason = "Invalid config: " + config.ValidationErrors();
        result.total_time_ms = total_timer.ElapsedMs();
        return result;
    }

    internal::Timer load_timer;

    io::FastaReader ref_reader(config.reference_fasta);
    std::vector<std::string> ref_names;
    std::vector<std::string> ref_seqs;

    while (ref_reader.HasMore()) {
        auto record = ref_reader.Next();
        ref_names.push_back(record.name);
        ref_seqs.push_back(record.sequence);
        result.reference_total_bp += record.sequence.size();
    }

    result.n_reference_seqs = ref_names.size();
    result.index_load_time_ms = load_timer.ElapsedMs();

    if (ref_names.empty()) {
        result.verdict_reason = "No reference sequences loaded";
        result.total_time_ms = total_timer.ElapsedMs();
        return result;
    }

    auto read_reader = io::FastqReader::Open(config.reads_fastq);
    if (!read_reader) {
        result.verdict_reason = "Failed to open reads file";
        result.total_time_ms = total_timer.ElapsedMs();
        return result;
    }

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    std::vector<uint32_t> read_lengths;

    while (read_reader->HasMore()) {
        auto record = read_reader->Next();
        if (record) {
            query_names.push_back(record->id);
            query_seqs.push_back(record->sequence);
            read_lengths.push_back(static_cast<uint32_t>(record->sequence.size()));
        }
    }

    result.n_reads = query_names.size();

    if (query_names.empty()) {
        result.verdict_reason = "No reads loaded";
        result.total_time_ms = total_timer.ElapsedMs();
        return result;
    }

    internal::Timer align_timer;

    classical::ClassicalPipelineConfig pipe_cfg;
    pipe_cfg.minimizer_config.k = config.minimizer_k;
    pipe_cfg.minimizer_config.w = config.minimizer_w;
    pipe_cfg.min_identity = config.min_identity;
    pipe_cfg.min_aligned_bases = 50;
    pipe_cfg.max_chains_to_extend = 5;

    auto align_results = classical::AlignWithClassicalPath(
        ref_names, ref_seqs, query_names, query_seqs, pipe_cfg);

    result.alignment_time_ms = align_timer.ElapsedMs();

    std::vector<AlignmentRecord> records;
    records.reserve(align_results.size());

    for (size_t i = 0; i < align_results.size(); ++i) {
        const auto& ar = align_results[i];
        uint32_t read_len = i < read_lengths.size() ? read_lengths[i] : 0;

        if (ar.HasAlignment()) {
            const auto* primary = ar.PrimaryAlignment();
            if (primary) {
                records.push_back(ConvertToAlignmentRecord(*primary, read_len));
                ++result.llmap_mapped;
            } else {
                records.push_back(make_unmapped(
                    ar.query_name, read_len, RejectionReason::NoSeeds));
                ++result.llmap_unmapped;
            }
        } else {
            records.push_back(make_unmapped(
                ar.query_name, read_len, RejectionReason::NoSeeds));
            ++result.llmap_unmapped;
        }
    }

    result.llmap_alignment_rate = static_cast<float>(result.llmap_mapped) /
                                   static_cast<float>(result.n_reads);

    internal::Timer val_timer;

    if (!config.minimap2_bam.empty() &&
        std::filesystem::exists(config.minimap2_bam)) {

        auto mm2_mapped = ParseMinimap2Bam(config.minimap2_bam);
        result.has_minimap2_baseline = !mm2_mapped.empty();

        if (result.has_minimap2_baseline) {
            for (const auto& [_, is_mapped] : mm2_mapped) {
                if (is_mapped) ++result.minimap2_mapped;
            }
            result.baseline_comparison = CompareToBaseline(records, mm2_mapped);
        }
    }

    if (!config.ground_truth_bed.empty() &&
        std::filesystem::exists(config.ground_truth_bed)) {

        auto real_truths = ParseGroundTruthBed(config.ground_truth_bed);
        result.has_ground_truth = !real_truths.empty();

        if (result.has_ground_truth) {
            std::vector<GroundTruth> truths;
            for (const auto& rt : real_truths) {
                GroundTruth gt;
                gt.read_id = rt.read_id;
                gt.true_start = rt.start;
                gt.true_end = rt.end;
                truths.push_back(std::move(gt));
            }

            KillSwitchValidator validator(config.thresholds);
            validator.SetPositionTolerance(config.position_tolerance);
            validator.LoadGroundTruth(truths);
            validator.AddRecords(records);

            result.validation = validator.Validate();
        }
    }

    result.validation_time_ms = val_timer.ElapsedMs();

    result.lossless_pass = (result.llmap_mapped + result.llmap_unmapped ==
                            result.n_reads);

    if (result.has_minimap2_baseline && result.minimap2_mapped > 0) {
        result.recall_pass = (result.baseline_comparison.recall_ratio >= 0.995);
    } else {
        result.recall_pass = true;
    }

    result.overall_pass = result.lossless_pass && result.recall_pass;

    if (!result.lossless_pass) {
        result.verdict_reason = "Silent read drop detected";
    } else if (!result.recall_pass) {
        std::ostringstream ss;
        ss << "Recall " << (result.baseline_comparison.recall_ratio * 100.0)
           << "% < 99.5% of minimap2";
        result.verdict_reason = ss.str();
    }

    result.total_time_ms = total_timer.ElapsedMs();
    return result;
}

}  // namespace llmap::validation
