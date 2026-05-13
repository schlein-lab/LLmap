// LLmap — Real reference validation implementation.

#include "validation/real_reference.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "classical/classical_pipeline.h"
#include "io/fasta_reader.h"
#include "io/fastq_reader.h"

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

std::vector<std::string_view> SplitLine(std::string_view line, char delim) {
    std::vector<std::string_view> fields;
    size_t start = 0;
    while (start < line.size()) {
        auto pos = line.find(delim, start);
        if (pos == std::string_view::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

uint64_t ParseUint64(std::string_view s) {
    uint64_t val = 0;
    for (char c : s) {
        if (c >= '0' && c <= '9') {
            val = val * 10 + (c - '0');
        }
    }
    return val;
}

std::string TrimWhitespace(std::string_view s) {
    size_t start = 0;
    while (start < s.size() && (s[start] == ' ' || s[start] == '\t' ||
                                 s[start] == '\n' || s[start] == '\r')) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && (s[end-1] == ' ' || s[end-1] == '\t' ||
                           s[end-1] == '\n' || s[end-1] == '\r')) {
        --end;
    }
    return std::string(s.substr(start, end - start));
}

}  // namespace

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

std::vector<RealGroundTruth> ParseGroundTruthBed(
    const std::filesystem::path& bed_path) {

    std::vector<RealGroundTruth> truths;
    std::ifstream in(bed_path);
    if (!in.is_open()) {
        return truths;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        auto fields = SplitLine(line, '\t');
        if (fields.size() < 4) continue;

        RealGroundTruth truth;
        truth.chrom = std::string(fields[0]);
        truth.start = ParseUint64(fields[1]);
        truth.end = ParseUint64(fields[2]);
        truth.read_id = std::string(fields[3]);

        if (fields.size() > 4) {
            truth.strand = std::string(fields[5]);
        }
        if (fields.size() > 6) {
            truth.gene = std::string(fields[6]);
        }

        truths.push_back(std::move(truth));
    }

    return truths;
}

std::unordered_map<std::string, bool> ParseMinimap2Bam(
    const std::filesystem::path& bam_path) {

    std::unordered_map<std::string, bool> mapped;

    // Shell out to samtools to parse BAM
    // Format: read_name\tflag
    std::string cmd = "samtools view " + bam_path.string() +
                      " 2>/dev/null | cut -f1,2";

    std::array<char, 4096> buffer;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(
        popen(cmd.c_str(), "r"), pclose);

    if (!pipe) {
        return mapped;
    }

    std::string output;
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        output += buffer.data();
    }

    std::istringstream iss(output);
    std::string line;
    while (std::getline(iss, line)) {
        auto fields = SplitLine(line, '\t');
        if (fields.size() < 2) continue;

        std::string read_name(fields[0]);
        uint64_t flag = ParseUint64(fields[1]);

        // Check unmapped flag (0x4)
        bool is_unmapped = (flag & 0x4) != 0;
        mapped[read_name] = !is_unmapped;
    }

    return mapped;
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
    Timer total_timer;
    RealReferenceResult result;

    // Validate config
    if (!config.Validate()) {
        result.verdict_reason = "Invalid config: " + config.ValidationErrors();
        result.total_time_ms = total_timer.ElapsedMs();
        return result;
    }

    // Load reference
    Timer load_timer;

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

    // Load reads
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

    // Run alignment
    Timer align_timer;

    classical::ClassicalPipelineConfig pipe_cfg;
    pipe_cfg.minimizer_config.k = config.minimizer_k;
    pipe_cfg.minimizer_config.w = config.minimizer_w;
    pipe_cfg.min_identity = config.min_identity;
    pipe_cfg.min_aligned_bases = 50;
    pipe_cfg.max_chains_to_extend = 5;

    auto align_results = classical::AlignWithClassicalPath(
        ref_names, ref_seqs, query_names, query_seqs, pipe_cfg);

    result.alignment_time_ms = align_timer.ElapsedMs();

    // Convert to AlignmentRecords
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

    // Load minimap2 baseline if provided
    Timer val_timer;

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

    // Load ground truth if provided
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

    // Determine kill-switch verdicts
    result.lossless_pass = (result.llmap_mapped + result.llmap_unmapped ==
                            result.n_reads);

    if (result.has_minimap2_baseline && result.minimap2_mapped > 0) {
        result.recall_pass = (result.baseline_comparison.recall_ratio >= 0.995);
    } else {
        result.recall_pass = true;  // No baseline to compare
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

std::string GenerateSlurmScript(
    const RealReferenceConfig& validation_config,
    const SlurmJobConfig& slurm_config) {

    std::ostringstream ss;
    ss << "#!/bin/bash\n";
    ss << "#SBATCH --job-name=" << slurm_config.job_name << "\n";
    ss << "#SBATCH --partition=" << slurm_config.partition << "\n";
    ss << "#SBATCH --gres=gpu:" << slurm_config.n_gpus << "\n";
    ss << "#SBATCH --cpus-per-task=" << slurm_config.n_cpus << "\n";
    ss << "#SBATCH --mem=" << slurm_config.memory << "\n";
    ss << "#SBATCH --time=" << slurm_config.time_limit << "\n";

    if (!slurm_config.output_log.empty()) {
        ss << "#SBATCH --output=" << slurm_config.output_log.string() << "\n";
    }
    if (!slurm_config.error_log.empty()) {
        ss << "#SBATCH --error=" << slurm_config.error_log.string() << "\n";
    }

    ss << "\n";
    ss << "# Load modules\n";
    ss << "module load cuda/12.3\n";
    ss << "module load gcc/13.2\n";
    ss << "\n";

    ss << "# Set up environment\n";
    ss << "export OMP_NUM_THREADS=" << slurm_config.n_cpus << "\n";
    ss << "\n";

    if (!slurm_config.work_dir.empty()) {
        ss << "cd " << slurm_config.work_dir.string() << "\n\n";
    }

    ss << "# Run validation\n";
    ss << "./build/llmap validate-real \\\n";
    ss << "    --reference " << validation_config.reference_fasta.string() << " \\\n";
    ss << "    --reads " << validation_config.reads_fastq.string() << " \\\n";

    if (!validation_config.minimap2_bam.empty()) {
        ss << "    --baseline " << validation_config.minimap2_bam.string() << " \\\n";
    }

    if (!validation_config.ground_truth_bed.empty()) {
        ss << "    --truth " << validation_config.ground_truth_bed.string() << " \\\n";
    }

    if (!validation_config.output_dir.empty()) {
        ss << "    --output-dir " << validation_config.output_dir.string() << " \\\n";
    }

    if (validation_config.use_gpu) {
        ss << "    --gpu " << validation_config.gpu_device << " \\\n";
    }

    ss << "    --k " << static_cast<int>(validation_config.minimizer_k) << " \\\n";
    ss << "    --w " << static_cast<int>(validation_config.minimizer_w) << " \\\n";
    ss << "    --min-identity " << validation_config.min_identity << "\n";
    ss << "\n";

    ss << "echo \"Job completed with exit code $?\"\n";

    return ss.str();
}

std::optional<std::string> SubmitSlurmJob(
    const std::filesystem::path& script_path) {

    std::string cmd = "sbatch " + script_path.string() + " 2>&1";

    std::array<char, 256> buffer;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(
        popen(cmd.c_str(), "r"), pclose);

    if (!pipe) {
        return std::nullopt;
    }

    std::string output;
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        output += buffer.data();
    }

    // Parse job ID from "Submitted batch job 12345"
    std::string prefix = "Submitted batch job ";
    auto pos = output.find(prefix);
    if (pos == std::string::npos) {
        return std::nullopt;
    }

    std::string job_id = TrimWhitespace(output.substr(pos + prefix.size()));
    if (job_id.empty()) {
        return std::nullopt;
    }

    return job_id;
}

SlurmJobStatus CheckSlurmJob(const std::string& job_id) {
    SlurmJobStatus status;
    status.job_id = job_id;

    std::string cmd = "squeue -j " + job_id + " -h -o '%T %e' 2>/dev/null";

    std::array<char, 256> buffer;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(
        popen(cmd.c_str(), "r"), pclose);

    if (!pipe) {
        status.state = "UNKNOWN";
        return status;
    }

    std::string output;
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        output += buffer.data();
    }

    output = TrimWhitespace(output);

    if (output.empty()) {
        // Job not in queue - check sacct for completion status
        cmd = "sacct -j " + job_id + " -n -o State,ExitCode -X 2>/dev/null";
        pipe.reset(popen(cmd.c_str(), "r"));

        output.clear();
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            output += buffer.data();
        }

        output = TrimWhitespace(output);
        auto fields = SplitLine(output, ' ');

        if (!fields.empty()) {
            status.state = std::string(fields[0]);
            status.is_complete = (status.state == "COMPLETED" ||
                                   status.state == "FAILED" ||
                                   status.state == "CANCELLED");

            if (fields.size() > 1) {
                auto exit_parts = SplitLine(fields[1], ':');
                if (!exit_parts.empty()) {
                    status.exit_code = static_cast<int>(ParseUint64(exit_parts[0]));
                }
            }
        } else {
            status.state = "COMPLETED";
            status.is_complete = true;
        }
    } else {
        auto fields = SplitLine(output, ' ');
        if (!fields.empty()) {
            status.state = std::string(fields[0]);
        }
        status.is_complete = false;
    }

    return status;
}

std::optional<RealReferenceResult> ParseSlurmOutput(
    const std::filesystem::path& output_path) {

    std::ifstream in(output_path);
    if (!in.is_open()) {
        return std::nullopt;
    }

    RealReferenceResult result;
    std::string line;

    while (std::getline(in, line)) {
        if (line.find("Sequences:") != std::string::npos) {
            auto pos = line.find_last_of(' ');
            if (pos != std::string::npos) {
                result.n_reference_seqs = static_cast<size_t>(
                    ParseUint64(line.substr(pos + 1)));
            }
        }
        else if (line.find("Total bp:") != std::string::npos) {
            auto pos = line.find_last_of(' ');
            if (pos != std::string::npos) {
                result.reference_total_bp = static_cast<size_t>(
                    ParseUint64(line.substr(pos + 1)));
            }
        }
        else if (line.find("LLmap Results:") != std::string::npos) {
            // Skip header line
        }
        else if (line.find("  Mapped:") != std::string::npos) {
            auto pos = line.find_last_of(' ');
            if (pos != std::string::npos) {
                result.llmap_mapped = static_cast<size_t>(
                    ParseUint64(line.substr(pos + 1)));
            }
        }
        else if (line.find("Overall:") != std::string::npos) {
            result.overall_pass = (line.find("PASS") != std::string::npos);
        }
    }

    return result;
}

}  // namespace llmap::validation
