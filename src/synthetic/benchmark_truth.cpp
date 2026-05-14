// LLmap — Benchmark-compatible ground truth generation implementation.

#include "synthetic/benchmark_truth.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace llmap::synthetic {

namespace {

constexpr char BASES[] = "ACGT";

class RandomGenerator {
public:
    explicit RandomGenerator(std::uint64_t seed) : rng_(seed) {}

    char random_base() {
        std::uniform_int_distribution<int> dist(0, 3);
        return BASES[dist(rng_)];
    }

    char mutate_base(char base) {
        char alternatives[3];
        int idx = 0;
        for (char b : BASES) {
            if (b != base) alternatives[idx++] = b;
        }
        std::uniform_int_distribution<int> dist(0, 2);
        return alternatives[dist(rng_)];
    }

    std::string random_sequence(std::uint64_t length) {
        std::string seq;
        seq.reserve(length);
        for (std::uint64_t i = 0; i < length; ++i) {
            seq.push_back(random_base());
        }
        return seq;
    }

    std::string apply_errors(const std::string& seq, float error_rate) {
        if (error_rate <= 0.0f) return seq;
        std::string result = seq;
        std::uniform_real_distribution<float> prob(0.0f, 1.0f);
        for (std::size_t i = 0; i < result.size(); ++i) {
            if (prob(rng_) < error_rate) {
                float et = prob(rng_);
                if (et < 0.7f) {
                    result[i] = mutate_base(result[i]);
                } else if (et < 0.85f && result.size() > 1) {
                    result.erase(i, 1);
                    if (i > 0) --i;
                } else {
                    result.insert(i, 1, random_base());
                    ++i;
                }
            }
        }
        return result;
    }

    std::string quality_string(std::uint64_t length, float error_rate) {
        std::string quals;
        quals.reserve(length);
        int mean_phred = static_cast<int>(-10.0 * std::log10(error_rate + 1e-10));
        mean_phred = std::clamp(mean_phred, 0, 50);
        std::normal_distribution<float> dist(static_cast<float>(mean_phred), 3.0f);
        for (std::uint64_t i = 0; i < length; ++i) {
            int q = std::clamp(static_cast<int>(dist(rng_)), 0, 50);
            quals.push_back(static_cast<char>(q + 33));
        }
        return quals;
    }

    int read_length(std::uint32_t mean, float stddev) {
        std::normal_distribution<float> dist(static_cast<float>(mean), stddev);
        return std::max(100, static_cast<int>(dist(rng_)));
    }

    std::uint64_t uniform_position(std::uint64_t max_val) {
        std::uniform_int_distribution<std::uint64_t> dist(0, max_val);
        return dist(rng_);
    }

    bool coin_flip(float prob = 0.5f) {
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        return dist(rng_) < prob;
    }

    std::size_t pick_index(std::size_t max_exclusive) {
        std::uniform_int_distribution<std::size_t> dist(0, max_exclusive - 1);
        return dist(rng_);
    }

private:
    std::mt19937_64 rng_;
};

std::string reverse_complement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': rc.push_back('T'); break;
            case 'T': rc.push_back('A'); break;
            case 'C': rc.push_back('G'); break;
            case 'G': rc.push_back('C'); break;
            default:  rc.push_back('N'); break;
        }
    }
    return rc;
}

}  // namespace

BenchmarkDataset BenchmarkGenerator::generate_t1(const T1Config& config) {
    BenchmarkDataset dataset;
    dataset.task = BenchmarkTask::T1_WGS;
    RandomGenerator rng(config.seed);

    // Generate reference sequences for each chromosome.
    for (const auto& chrom : config.chromosomes) {
        std::string seq = rng.random_sequence(config.region_length);
        dataset.references.emplace_back(chrom, std::move(seq));
    }

    // Calculate reads per chromosome based on coverage.
    std::uint64_t total_ref_len = config.region_length * config.chromosomes.size();
    std::uint64_t bases_needed = total_ref_len * config.coverage;
    std::uint64_t reads_to_gen = bases_needed / config.read_length;
    reads_to_gen = std::min(reads_to_gen, config.n_reads);

    dataset.reads.reserve(reads_to_gen);

    for (std::uint64_t i = 0; i < reads_to_gen; ++i) {
        BenchmarkRead read;

        // Pick a random chromosome.
        std::size_t chrom_idx = rng.pick_index(dataset.references.size());
        const auto& [chrom_name, chrom_seq] = dataset.references[chrom_idx];

        // Generate read length.
        int len = rng.read_length(config.read_length, config.read_length_stddev);
        len = std::min(len, static_cast<int>(chrom_seq.size()));

        // Pick a random position.
        std::uint64_t max_start = chrom_seq.size() > static_cast<std::size_t>(len)
                                      ? chrom_seq.size() - len : 0;
        std::uint64_t start = rng.uniform_position(max_start);

        // Extract sequence and apply errors.
        std::string raw_seq = chrom_seq.substr(start, len);
        read.is_forward = rng.coin_flip();
        if (!read.is_forward) {
            raw_seq = reverse_complement(raw_seq);
        }

        read.sequence = rng.apply_errors(raw_seq, config.error_rate);
        read.quality = rng.quality_string(read.sequence.size(), config.error_rate);

        // Set ground truth.
        read.chrom = chrom_name;
        read.pos = static_cast<std::int64_t>(start);

        // Generate ID with embedded truth (for debugging).
        std::ostringstream id;
        id << "synth_t1_" << i << "_" << chrom_name << "_" << start;
        read.id = id.str();

        dataset.reads.push_back(std::move(read));
        dataset.total_bases += len;
    }

    dataset.total_reads = dataset.reads.size();
    return dataset;
}

BenchmarkDataset BenchmarkGenerator::generate_t2(const T2Config& config) {
    BenchmarkDataset dataset;
    dataset.task = BenchmarkTask::T2_ParalogStress;
    RandomGenerator rng(config.seed);

    // Generate reference sequences for each paralog family member.
    for (const auto& family : config.families) {
        // Generate a base sequence for the family.
        std::string base_seq = rng.random_sequence(family.region_length);

        // Create divergent copies for each member.
        for (const auto& member : family.members) {
            std::string member_seq = base_seq;

            // Apply divergence (1 - seq_identity fraction of positions mutated).
            float divergence = 1.0f - family.seq_identity;
            for (std::size_t pos = 0; pos < member_seq.size(); ++pos) {
                if (rng.coin_flip(divergence)) {
                    member_seq[pos] = rng.mutate_base(member_seq[pos]);
                }
            }

            dataset.references.emplace_back(member, std::move(member_seq));
        }
    }

    // Early return if no families/references to generate reads from.
    if (dataset.references.empty()) {
        return dataset;
    }

    // Calculate reads per family member.
    std::uint64_t reads_per_member = config.n_reads / dataset.references.size();

    dataset.reads.reserve(config.n_reads);

    std::uint64_t read_idx = 0;
    for (const auto& [member_name, member_seq] : dataset.references) {
        for (std::uint64_t i = 0; i < reads_per_member && read_idx < config.n_reads; ++i) {
            BenchmarkRead read;

            // Generate read length.
            int len = rng.read_length(config.read_length, config.read_length_stddev);
            len = std::min(len, static_cast<int>(member_seq.size()));

            // Pick a random position.
            std::uint64_t max_start = member_seq.size() > static_cast<std::size_t>(len)
                                          ? member_seq.size() - len : 0;
            std::uint64_t start = rng.uniform_position(max_start);

            // Extract sequence and apply errors.
            std::string raw_seq = member_seq.substr(start, len);
            read.is_forward = rng.coin_flip();
            if (!read.is_forward) {
                raw_seq = reverse_complement(raw_seq);
            }

            read.sequence = rng.apply_errors(raw_seq, config.error_rate);
            read.quality = rng.quality_string(read.sequence.size(), config.error_rate);

            // Set ground truth.
            read.chrom = member_name;
            read.pos = static_cast<std::int64_t>(start);
            read.paralog = member_name;

            // Generate ID.
            std::ostringstream id;
            id << "synth_t2_" << read_idx << "_" << member_name << "_" << start;
            read.id = id.str();

            dataset.reads.push_back(std::move(read));
            dataset.total_bases += len;
            ++read_idx;
        }
    }

    dataset.total_reads = dataset.reads.size();
    return dataset;
}

void BenchmarkGenerator::write_fastq(const BenchmarkDataset& dataset,
                                     const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open FASTQ file: " + path);
    }

    for (const auto& read : dataset.reads) {
        out << '@' << read.id << '\n'
            << read.sequence << '\n'
            << '+' << '\n'
            << read.quality << '\n';
    }
}

void BenchmarkGenerator::write_reference(const BenchmarkDataset& dataset,
                                         const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open FASTA file: " + path);
    }

    for (const auto& [name, seq] : dataset.references) {
        out << '>' << name << '\n';
        // Wrap at 80 columns.
        for (std::size_t i = 0; i < seq.size(); i += 80) {
            out << seq.substr(i, 80) << '\n';
        }
    }
}

void BenchmarkGenerator::write_truth_tsv(const BenchmarkDataset& dataset,
                                         const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open truth TSV file: " + path);
    }

    // Header comment (skipped by compute.py).
    out << "# read_id\tchrom\tpos\n";

    for (const auto& read : dataset.reads) {
        out << read.id << '\t' << read.chrom << '\t' << read.pos << '\n';
    }
}

void BenchmarkGenerator::write_paralog_truth_tsv(const BenchmarkDataset& dataset,
                                                 const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open paralog truth TSV file: " + path);
    }

    // Header comment.
    out << "# read_id\tparalog\tchrom\tpos\n";

    for (const auto& read : dataset.reads) {
        out << read.id << '\t'
            << read.paralog << '\t'
            << read.chrom << '\t'
            << read.pos << '\n';
    }
}

namespace paralog_presets {

ParalogFamily igh_constant() {
    return ParalogFamily{
        .name = "IGH_constant",
        .members = {"IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP"},
        .region_length = 50000,
        .n_psvs = 20,
        .seq_identity = 0.95f,
    };
}

ParalogFamily nphp1() {
    return ParalogFamily{
        .name = "NPHP1",
        .members = {"NPHP1", "NPHP1B"},
        .region_length = 80000,
        .n_psvs = 30,
        .seq_identity = 0.98f,
    };
}

ParalogFamily mhc_class1() {
    return ParalogFamily{
        .name = "MHC_class1",
        .members = {"HLA-A", "HLA-A-pseudo"},
        .region_length = 30000,
        .n_psvs = 15,
        .seq_identity = 0.92f,
    };
}

}  // namespace paralog_presets

}  // namespace llmap::synthetic
