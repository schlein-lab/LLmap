// LLmap — Synthetic IGH-locus generator implementation.

#include "synthetic/igh_locus_generator.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace llmap::synthetic {

namespace {

constexpr char BASES[] = "ACGT";

}  // namespace

IGHLocusGenerator::IGHLocusGenerator(IGHLocusConfig config)
    : config_(std::move(config)), rng_(config_.seed) {}

char IGHLocusGenerator::random_base() {
    std::uniform_int_distribution<int> dist(0, 3);
    return BASES[dist(rng_)];
}

char IGHLocusGenerator::mutate_base(char base) {
    char alternatives[3];
    int idx = 0;
    for (char b : BASES) {
        if (b != base) {
            alternatives[idx++] = b;
        }
    }
    std::uniform_int_distribution<int> dist(0, 2);
    return alternatives[dist(rng_)];
}

std::string IGHLocusGenerator::generate_random_sequence(std::uint64_t length) {
    std::string seq;
    seq.reserve(length);
    for (std::uint64_t i = 0; i < length; ++i) {
        seq.push_back(random_base());
    }
    return seq;
}

std::string IGHLocusGenerator::apply_errors(const std::string& seq, float error_rate) {
    if (error_rate <= 0.0f) {
        return seq;
    }

    std::string result = seq;
    std::uniform_real_distribution<float> prob_dist(0.0f, 1.0f);

    for (std::size_t i = 0; i < result.size(); ++i) {
        if (prob_dist(rng_) < error_rate) {
            std::uniform_real_distribution<float> error_type(0.0f, 1.0f);
            float et = error_type(rng_);

            if (et < 0.7f) {
                // Substitution (most common in HiFi)
                result[i] = mutate_base(result[i]);
            } else if (et < 0.85f) {
                // Deletion (skip this base by shifting)
                result.erase(i, 1);
                if (i > 0) --i;
            } else {
                // Insertion
                result.insert(i, 1, random_base());
                ++i;
            }
        }
    }

    return result;
}

std::string IGHLocusGenerator::generate_quality_string(std::uint64_t length, float mean_error_rate) {
    std::string quals;
    quals.reserve(length);

    int mean_phred = static_cast<int>(-10.0 * std::log10(mean_error_rate + 1e-10));
    mean_phred = std::clamp(mean_phred, 0, 50);

    std::normal_distribution<float> qual_dist(
        static_cast<float>(mean_phred), 3.0f);

    for (std::uint64_t i = 0; i < length; ++i) {
        int q = std::clamp(static_cast<int>(qual_dist(rng_)), 0, 50);
        quals.push_back(static_cast<char>(q + 33));  // Phred+33
    }

    return quals;
}

GeneratedLocus IGHLocusGenerator::generate_locus() {
    GeneratedLocus locus;

    // Generate canonical sequence
    locus.canonical_sequence = generate_random_sequence(config_.locus_length);

    // Start with duplicate as copy of canonical
    locus.duplicate_sequence = locus.canonical_sequence;

    // Determine sequence-identical regions
    std::uint64_t identical_length = static_cast<std::uint64_t>(
        config_.locus_length * config_.seq_identical_fraction);

    if (identical_length > 0) {
        // Place identical region in the middle (simulates sequence-identical exons)
        std::uint64_t start = (config_.locus_length - identical_length) / 2;
        locus.identical_regions.push_back({start, start + identical_length});
    }

    // Plant PSVs outside identical regions
    std::vector<std::uint64_t> available_positions;
    for (std::uint64_t pos = 0; pos < config_.locus_length; ++pos) {
        bool in_identical = false;
        for (const auto& region : locus.identical_regions) {
            if (pos >= region.start && pos < region.end) {
                in_identical = true;
                break;
            }
        }
        if (!in_identical) {
            available_positions.push_back(pos);
        }
    }

    // Shuffle and pick PSV positions
    std::shuffle(available_positions.begin(), available_positions.end(), rng_);

    std::uint32_t n_psvs = std::min(config_.n_psvs,
                                     static_cast<std::uint32_t>(available_positions.size()));

    std::vector<std::uint64_t> psv_positions(
        available_positions.begin(),
        available_positions.begin() + n_psvs);
    std::sort(psv_positions.begin(), psv_positions.end());

    // Assign gene contexts based on position
    std::uint64_t segment_size = config_.locus_length / config_.gene_names.size();

    for (std::uint64_t pos : psv_positions) {
        char canonical_allele = locus.canonical_sequence[pos];
        char dup_allele = mutate_base(canonical_allele);

        // Apply the mutation to the duplicate
        locus.duplicate_sequence[pos] = dup_allele;

        std::size_t gene_idx = std::min(
            pos / segment_size,
            config_.gene_names.size() - 1);

        locus.psvs.push_back({
            .position = pos,
            .canonical_allele = canonical_allele,
            .duplicate_allele = dup_allele,
            .gene_context = config_.gene_names[gene_idx],
            .is_discriminating = true,
        });
    }

    // Add non-discriminating "PSVs" in identical regions for tracking
    for (const auto& region : locus.identical_regions) {
        std::uint64_t mid = (region.start + region.end) / 2;
        locus.psvs.push_back({
            .position = mid,
            .canonical_allele = locus.canonical_sequence[mid],
            .duplicate_allele = locus.canonical_sequence[mid],  // same!
            .gene_context = "IDENTICAL_REGION",
            .is_discriminating = false,
        });
    }

    return locus;
}

std::vector<SyntheticRead> IGHLocusGenerator::generate_reads(const GeneratedLocus& locus) {
    std::vector<SyntheticRead> reads;

    // Calculate effective depths
    float dup_fraction = config_.mosaic.dup_fraction;
    std::uint32_t total_depth = config_.mosaic.canonical_depth;

    if (config_.mosaic.dup_depth > 0) {
        total_depth += config_.mosaic.dup_depth;
        dup_fraction = static_cast<float>(config_.mosaic.dup_depth) /
                       static_cast<float>(total_depth);
    }

    // Estimate number of reads needed to achieve coverage
    std::uint64_t total_bases_needed = config_.locus_length * total_depth;
    std::uint64_t approx_reads = total_bases_needed / config_.read_length;

    std::normal_distribution<float> len_dist(
        static_cast<float>(config_.read_length),
        config_.read_length_stddev);

    std::uniform_real_distribution<float> origin_dist(0.0f, 1.0f);
    std::uniform_int_distribution<std::uint64_t> canonical_pos_dist(
        0, config_.locus_length > config_.read_length ?
           config_.locus_length - config_.read_length : 0);
    std::uniform_int_distribution<std::uint64_t> dup_pos_dist(
        0, config_.locus_length > config_.read_length ?
           config_.locus_length - config_.read_length : 0);

    for (std::uint64_t i = 0; i < approx_reads; ++i) {
        SyntheticRead read;

        // Determine origin based on dup_fraction
        bool from_dup = origin_dist(rng_) < dup_fraction;
        read.origin = from_dup ? SyntheticRead::Origin::Duplicate
                               : SyntheticRead::Origin::Canonical;

        // Determine read length
        int read_len = std::clamp(
            static_cast<int>(len_dist(rng_)),
            100,
            static_cast<int>(config_.locus_length));

        // Determine start position
        std::uint64_t max_start = config_.locus_length > static_cast<std::uint64_t>(read_len)
                                      ? config_.locus_length - read_len
                                      : 0;
        std::uniform_int_distribution<std::uint64_t> start_dist(0, max_start);
        std::uint64_t start = start_dist(rng_);

        read.true_start = start;
        read.true_end = start + read_len;

        // Assign gene context
        std::uint64_t segment_size = config_.locus_length / config_.gene_names.size();
        std::size_t gene_idx = std::min(
            start / segment_size,
            config_.gene_names.size() - 1);
        read.source_gene = config_.gene_names[gene_idx];

        // Extract sequence from appropriate copy
        const std::string& source = from_dup ? locus.duplicate_sequence
                                             : locus.canonical_sequence;
        std::string raw_seq = source.substr(start, read_len);

        // Apply sequencing errors
        read.sequence = apply_errors(raw_seq, config_.error_rate);

        // Generate quality string matching the (post-error) sequence length
        read.quality = generate_quality_string(read.sequence.size(), config_.error_rate);

        // Track covered PSVs
        for (const auto& psv : locus.psvs) {
            if (psv.position >= start && psv.position < start + read_len) {
                read.covered_psv_positions.push_back(psv.position);
            }
        }

        // Generate read ID encoding ground truth for validation
        std::ostringstream id;
        id << "synth_" << i
           << "_" << (from_dup ? "DUP" : "CAN")
           << "_" << start << "-" << (start + read_len)
           << "_" << read.source_gene;
        read.id = id.str();

        reads.push_back(std::move(read));
    }

    return reads;
}

GeneratedDataset IGHLocusGenerator::generate() {
    GeneratedDataset dataset;

    dataset.locus = generate_locus();
    dataset.reads = generate_reads(dataset.locus);

    // Compute summary statistics
    for (const auto& read : dataset.reads) {
        if (read.origin == SyntheticRead::Origin::Canonical) {
            ++dataset.n_canonical_reads;
        } else {
            ++dataset.n_dup_reads;
        }
        dataset.total_psv_observations += read.covered_psv_positions.size();
    }

    if (!dataset.reads.empty()) {
        dataset.actual_dup_fraction = static_cast<float>(dataset.n_dup_reads) /
                                      static_cast<float>(dataset.reads.size());
    }

    return dataset;
}

void IGHLocusGenerator::write_fastq(const std::vector<SyntheticRead>& reads,
                                    const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open FASTQ file for writing: " + path);
    }

    for (const auto& read : reads) {
        out << '@' << read.id << '\n'
            << read.sequence << '\n'
            << '+' << '\n'
            << read.quality << '\n';
    }
}

void IGHLocusGenerator::write_ground_truth(const GeneratedDataset& dataset,
                                           const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open ground truth file for writing: " + path);
    }

    out << "read_id\torigin\ttrue_start\ttrue_end\tsource_gene\tn_psvs_covered\n";

    for (const auto& read : dataset.reads) {
        out << read.id << '\t'
            << (read.origin == SyntheticRead::Origin::Duplicate ? "DUP" : "CAN") << '\t'
            << read.true_start << '\t'
            << read.true_end << '\t'
            << read.source_gene << '\t'
            << read.covered_psv_positions.size() << '\n';
    }
}

// Preset implementations

namespace presets {

IGHLocusConfig canonical_only(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.0f;
    cfg.mosaic.canonical_depth = 30;
    return cfg;
}

IGHLocusConfig balanced_mosaic(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.5f;
    cfg.mosaic.canonical_depth = 15;
    cfg.mosaic.dup_depth = 15;
    return cfg;
}

IGHLocusConfig dup_fraction_5(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.05f;
    return cfg;
}

IGHLocusConfig dup_fraction_10(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.10f;
    return cfg;
}

IGHLocusConfig dup_fraction_30(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.30f;
    return cfg;
}

IGHLocusConfig dup_fraction_50(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.50f;
    return cfg;
}

IGHLocusConfig dup_fraction_100(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 1.0f;
    return cfg;
}

IGHLocusConfig seq_identical_stress(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.seq_identical_fraction = 0.5f;  // half the locus is identical
    cfg.mosaic.dup_fraction = 0.5f;
    cfg.n_psvs = 10;  // fewer PSVs because less space
    return cfg;
}

IGHLocusConfig tiny_test(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.locus_length = 5000;
    cfg.n_psvs = 5;
    cfg.read_length = 1000;
    cfg.read_length_stddev = 100.0f;
    cfg.mosaic.canonical_depth = 10;
    cfg.mosaic.dup_fraction = 0.3f;
    return cfg;
}

}  // namespace presets

}  // namespace llmap::synthetic
