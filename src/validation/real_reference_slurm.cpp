// LLmap — SLURM job management for real reference validation.
//
// Generates scripts, submits jobs, checks status, parses output.

#include "validation/real_reference.h"
#include "validation/real_reference_internal.h"

#include <array>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>

namespace llmap::validation {

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

    std::string prefix = "Submitted batch job ";
    auto pos = output.find(prefix);
    if (pos == std::string::npos) {
        return std::nullopt;
    }

    std::string job_id = internal::TrimWhitespace(output.substr(pos + prefix.size()));
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

    output = internal::TrimWhitespace(output);

    if (output.empty()) {
        cmd = "sacct -j " + job_id + " -n -o State,ExitCode -X 2>/dev/null";
        pipe.reset(popen(cmd.c_str(), "r"));

        output.clear();
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            output += buffer.data();
        }

        output = internal::TrimWhitespace(output);
        auto fields = internal::SplitLine(output, ' ');

        if (!fields.empty()) {
            status.state = std::string(fields[0]);
            status.is_complete = (status.state == "COMPLETED" ||
                                   status.state == "FAILED" ||
                                   status.state == "CANCELLED");

            if (fields.size() > 1) {
                auto exit_parts = internal::SplitLine(fields[1], ':');
                if (!exit_parts.empty()) {
                    status.exit_code = static_cast<int>(internal::ParseUint64(exit_parts[0]));
                }
            }
        } else {
            status.state = "COMPLETED";
            status.is_complete = true;
        }
    } else {
        auto fields = internal::SplitLine(output, ' ');
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
                    internal::ParseUint64(line.substr(pos + 1)));
            }
        }
        else if (line.find("Total bp:") != std::string::npos) {
            auto pos = line.find_last_of(' ');
            if (pos != std::string::npos) {
                result.reference_total_bp = static_cast<size_t>(
                    internal::ParseUint64(line.substr(pos + 1)));
            }
        }
        else if (line.find("LLmap Results:") != std::string::npos) {
            // Skip header line
        }
        else if (line.find("  Mapped:") != std::string::npos) {
            auto pos = line.find_last_of(' ');
            if (pos != std::string::npos) {
                result.llmap_mapped = static_cast<size_t>(
                    internal::ParseUint64(line.substr(pos + 1)));
            }
        }
        else if (line.find("Overall:") != std::string::npos) {
            result.overall_pass = (line.find("PASS") != std::string::npos);
        }
    }

    return result;
}

}  // namespace llmap::validation
