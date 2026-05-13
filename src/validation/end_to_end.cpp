// LLmap — End-to-end synthetic validation implementation.

#include "validation/killswitch.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <sstream>

#include "classical/classical_pipeline.h"
#include "synthetic/igh_locus_generator.h"

namespace llmap::validation {

namespace {

class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    float ElapsedMs() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<float, std::milli>(now - start_).count();
    }

private:
    std::chrono::steady_clock::time_point start_;
};

CigarString ConvertCigar(const std::vector<classical::CigarElement>& cigar) {
    std::ostringstream oss;
    for (const auto& elem : cigar) {
        oss << elem.ToString();
    }
    return CigarString{oss.str()};
}

}  // namespace

AlignmentRecord ConvertToAlignmentRecord(
    const classical::ClassicalAlignment& aln,
    uint32_t read_len) {

    AlignmentHit hit;
    hit.target_id = aln.ref_name;
    hit.start = static_cast<uint64_t>(aln.ref_start);
    hit.end = static_cast<uint64_t>(aln.ref_end);
    hit.cigar = ConvertCigar(aln.cigar);
    hit.score = aln.score;
    hit.nm = static_cast<uint32_t>((1.0f - aln.identity) * aln.AlignedBases());

    return make_mapped(aln.query_name, read_len, std::move(hit));
}

std::vector<AlignmentRecord> ConvertResults(
    const std::vector<classical::ReadAlignmentResult>& results,
    const std::vector<uint32_t>& read_lengths) {

    std::vector<AlignmentRecord> records;
    records.reserve(results.size());

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        uint32_t read_len = i < read_lengths.size() ? read_lengths[i] : 0;

        if (result.HasAlignment()) {
            const auto* primary = result.PrimaryAlignment();
            if (primary) {
                records.push_back(ConvertToAlignmentRecord(*primary, read_len));
            }
        } else {
            records.push_back(make_unmapped(
                result.query_name, read_len, RejectionReason::NoSeeds));
        }
    }

    return records;
}

EndToEndConfig EndToEndConfig::Minimal() {
    EndToEndConfig cfg;
    cfg.locus_length = 10000;
    cfg.n_psvs = 5;
    cfg.read_length = 2000;
    cfg.n_reads = 10;
    cfg.min_position_accuracy = 0.80;
    return cfg;
}

EndToEndConfig EndToEndConfig::Standard() {
    EndToEndConfig cfg;
    cfg.locus_length = 30000;
    cfg.n_psvs = 15;
    cfg.read_length = 5000;
    cfg.n_reads = 50;
    cfg.min_position_accuracy = 0.90;
    return cfg;
}

EndToEndConfig EndToEndConfig::Stress() {
    EndToEndConfig cfg;
    cfg.locus_length = 100000;
    cfg.n_psvs = 50;
    cfg.read_length = 10000;
    cfg.n_reads = 200;
    cfg.min_position_accuracy = 0.95;
    return cfg;
}

std::string EndToEndResult::Summary() const {
    std::ostringstream ss;
    ss << "=== End-to-End Synthetic Validation ===\n\n";

    ss << "Synthetic Data:\n";
    ss << "  Reads generated: " << n_reads_generated << "\n";
    ss << "  Canonical:       " << n_canonical << "\n";
    ss << "  Duplicate:       " << n_duplicate << "\n";
    ss << "  Dup fraction:    " << (actual_dup_fraction * 100.0f) << "%\n\n";

    ss << "Alignment Results:\n";
    ss << "  Aligned:         " << n_aligned << "\n";
    ss << "  Unmapped:        " << n_unmapped << "\n";
    ss << "  Alignment rate:  " << (alignment_rate * 100.0f) << "%\n";
    ss << "  Avg identity:    " << (avg_identity * 100.0f) << "%\n\n";

    ss << "Position Accuracy:\n";
    ss << "  Accuracy:        " << (position_accuracy * 100.0) << "%\n";
    ss << "  Mean error:      " << mean_position_error << " bp\n";
    ss << "  Max error:       " << max_position_error << " bp\n\n";

    ss << "Timing:\n";
    ss << "  Generation:      " << generation_time_ms << " ms\n";
    ss << "  Alignment:       " << alignment_time_ms << " ms\n";
    ss << "  Validation:      " << validation_time_ms << " ms\n";
    ss << "  Total:           " << total_time_ms << " ms\n\n";

    ss << "Verdict: " << (passed ? "PASS" : "FAIL");
    if (!verdict_reason.empty()) {
        ss << " (" << verdict_reason << ")";
    }
    ss << "\n";

    return ss.str();
}

EndToEndResult RunEndToEndValidation(const EndToEndConfig& config) {
    Timer total_timer;
    EndToEndResult result;

    // Phase 1: Generate synthetic data
    Timer gen_timer;

    synthetic::IGHLocusConfig locus_cfg;
    locus_cfg.seed = config.seed;
    locus_cfg.locus_length = config.locus_length;
    locus_cfg.n_psvs = config.n_psvs;
    locus_cfg.read_length = config.read_length;
    locus_cfg.mosaic.dup_fraction = config.dup_fraction;

    // Compute approximate depths from dup_fraction and n_reads
    uint32_t total_depth = config.n_reads;
    locus_cfg.mosaic.canonical_depth =
        static_cast<uint32_t>(total_depth * (1.0f - config.dup_fraction));
    locus_cfg.mosaic.dup_depth =
        static_cast<uint32_t>(total_depth * config.dup_fraction);
    if (locus_cfg.mosaic.canonical_depth == 0 && config.dup_fraction < 1.0f) {
        locus_cfg.mosaic.canonical_depth = 1;
    }
    if (locus_cfg.mosaic.dup_depth == 0 && config.dup_fraction > 0.0f) {
        locus_cfg.mosaic.dup_depth = 1;
    }

    synthetic::IGHLocusGenerator generator(locus_cfg);
    auto dataset = generator.generate();

    result.n_reads_generated = dataset.reads.size();
    result.n_canonical = dataset.n_canonical_reads;
    result.n_duplicate = dataset.n_dup_reads;
    result.actual_dup_fraction = dataset.actual_dup_fraction;
    result.generation_time_ms = gen_timer.ElapsedMs();

    if (dataset.reads.empty()) {
        result.verdict_reason = "No reads generated";
        result.passed = false;
        result.total_time_ms = total_timer.ElapsedMs();
        return result;
    }

    // Phase 2: Build reference and align
    Timer align_timer;

    // Use canonical sequence as reference
    std::vector<std::string> ref_names = {"synthetic_locus"};
    std::vector<std::string> ref_seqs = {dataset.locus.canonical_sequence};

    std::vector<std::string> query_names;
    std::vector<std::string> query_seqs;
    std::vector<uint32_t> read_lengths;

    for (const auto& read : dataset.reads) {
        query_names.push_back(read.id);
        query_seqs.push_back(read.sequence);
        read_lengths.push_back(static_cast<uint32_t>(read.sequence.size()));
    }

    classical::ClassicalPipelineConfig pipe_cfg;
    pipe_cfg.minimizer_config.k = config.minimizer_k;
    pipe_cfg.minimizer_config.w = config.minimizer_w;
    pipe_cfg.min_identity = config.min_identity;
    pipe_cfg.min_aligned_bases = 50;
    pipe_cfg.max_chains_to_extend = 5;
    pipe_cfg.report_secondary = false;

    auto align_results = classical::AlignWithClassicalPath(
        ref_names, ref_seqs, query_names, query_seqs, pipe_cfg);

    result.alignment_time_ms = align_timer.ElapsedMs();

    // Compute alignment statistics
    float total_identity = 0.0f;
    for (const auto& ar : align_results) {
        if (ar.HasAlignment()) {
            ++result.n_aligned;
            if (!ar.alignments.empty()) {
                total_identity += ar.alignments[0].identity;
            }
        } else {
            ++result.n_unmapped;
        }
    }

    result.alignment_rate = static_cast<float>(result.n_aligned) /
                           static_cast<float>(result.n_reads_generated);
    result.avg_identity = result.n_aligned > 0
        ? total_identity / static_cast<float>(result.n_aligned)
        : 0.0f;

    // Phase 3: Validate through kill-switch framework
    Timer val_timer;

    auto records = ConvertResults(align_results, read_lengths);

    KillSwitchThresholds thresholds;
    thresholds.require_lossless = config.require_lossless;

    KillSwitchValidator validator(thresholds);

    PositionTolerance pos_tol;
    pos_tol.base_tolerance = static_cast<int64_t>(config.position_tolerance_bp);
    pos_tol.relative_tolerance = 0.05f;
    validator.SetPositionTolerance(pos_tol);

    validator.LoadGroundTruth(dataset);
    validator.AddRecords(records);

    result.validation = validator.Validate();
    result.position_accuracy = result.validation.position_accuracy;

    // Compute detailed position error stats
    auto detailed = validator.ValidateDetailed();
    std::vector<int64_t> errors;
    for (const auto& v : detailed) {
        if (v.has_ground_truth && v.has_alignment_record &&
            v.mapped_status == AlignmentStatus::Mapped) {
            errors.push_back(std::abs(v.position_error));
        }
    }

    if (!errors.empty()) {
        int64_t sum = std::accumulate(errors.begin(), errors.end(), int64_t{0});
        result.mean_position_error = sum / static_cast<int64_t>(errors.size());
        result.max_position_error = *std::max_element(errors.begin(), errors.end());
    }

    result.validation_time_ms = val_timer.ElapsedMs();
    result.total_time_ms = total_timer.ElapsedMs();

    // Determine verdict
    result.passed = true;

    if (config.require_lossless && !result.validation.lossless_pass) {
        result.passed = false;
        result.verdict_reason = "Lossless invariant violated";
    } else if (result.position_accuracy < config.min_position_accuracy) {
        result.passed = false;
        std::ostringstream ss;
        ss << "Position accuracy " << (result.position_accuracy * 100.0)
           << "% < required " << (config.min_position_accuracy * 100.0) << "%";
        result.verdict_reason = ss.str();
    }

    return result;
}

EndToEndResult RunMinimalValidation(uint64_t seed) {
    auto cfg = EndToEndConfig::Minimal();
    cfg.seed = seed;
    return RunEndToEndValidation(cfg);
}

EndToEndResult RunStandardValidation(uint64_t seed) {
    auto cfg = EndToEndConfig::Standard();
    cfg.seed = seed;
    return RunEndToEndValidation(cfg);
}

}  // namespace llmap::validation
