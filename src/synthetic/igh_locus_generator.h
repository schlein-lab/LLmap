// LLmap — Synthetic IGH-locus generator.
//
// Generates synthetic IGH locus sequences with planted PSVs, controlled
// mosaic duplication fractions, and sequence-identical exon edge cases.
// Deterministic by seed for reproducible benchmarking.

#pragma once

#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace llmap::synthetic {

// A PSV (Paralog-Specific Variant) planted in the synthetic locus.
struct PlantedPSV {
    std::uint64_t position;         // position in canonical sequence
    char canonical_allele;
    char duplicate_allele;
    std::string gene_context;       // e.g. "IGHG1", "IGHG4"
    bool is_discriminating{true};   // false for sequence-identical regions
};

// A synthetic read with ground truth.
struct SyntheticRead {
    std::string id;
    std::string sequence;
    std::string quality;            // Phred+33 encoded

    // Ground truth for validation
    enum class Origin { Canonical, Duplicate };
    Origin origin;
    std::uint64_t true_start;
    std::uint64_t true_end;
    std::string source_gene;

    std::vector<std::uint64_t> covered_psv_positions;
};

// Configuration for mosaic duplication simulation.
struct MosaicConfig {
    // Fraction of reads that come from the duplicated copy (0.0 to 1.0).
    // 0.0 = pure canonical, 1.0 = pure duplicate, 0.5 = balanced mosaic
    float dup_fraction{0.0f};

    // Coverage depth for each copy; final coverage = canonical_depth + dup_depth
    std::uint32_t canonical_depth{30};
    std::uint32_t dup_depth{0};  // computed from dup_fraction if not set
};

// Configuration for the synthetic locus generator.
struct IGHLocusConfig {
    // Seed for deterministic generation.
    std::uint64_t seed{42};

    // Locus length in bp (real IGH ~1.2Mbp, we use smaller for tests).
    std::uint64_t locus_length{50000};

    // Number of PSVs to plant.
    std::uint32_t n_psvs{20};

    // Fraction of locus that is sequence-identical between copies (edge case).
    float seq_identical_fraction{0.1f};

    // Mosaic duplication configuration.
    MosaicConfig mosaic;

    // Read parameters.
    std::uint32_t read_length{10000};    // HiFi-length default
    float read_length_stddev{1000.0f};
    float error_rate{0.001f};            // HiFi Q30

    // Gene segment definitions (simplified).
    std::vector<std::string> gene_names{"IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP"};
};

// Output from generation.
struct GeneratedLocus {
    std::string canonical_sequence;
    std::string duplicate_sequence;
    std::vector<PlantedPSV> psvs;

    // Regions that are sequence-identical (edge case for coverage-based CNV).
    struct IdenticalRegion {
        std::uint64_t start;
        std::uint64_t end;
    };
    std::vector<IdenticalRegion> identical_regions;
};

struct GeneratedDataset {
    GeneratedLocus locus;
    std::vector<SyntheticRead> reads;

    // Summary statistics for validation.
    std::uint64_t n_canonical_reads{0};
    std::uint64_t n_dup_reads{0};
    float actual_dup_fraction{0.0f};
    std::uint64_t total_psv_observations{0};
};

// The main generator class.
class IGHLocusGenerator {
public:
    explicit IGHLocusGenerator(IGHLocusConfig config);

    // Generate the locus sequences with planted PSVs.
    [[nodiscard]] GeneratedLocus generate_locus();

    // Generate reads from a locus.
    [[nodiscard]] std::vector<SyntheticRead> generate_reads(const GeneratedLocus& locus);

    // Convenience: generate everything.
    [[nodiscard]] GeneratedDataset generate();

    // Write reads to FASTQ format.
    static void write_fastq(const std::vector<SyntheticRead>& reads,
                            const std::string& path);

    // Write ground truth to TSV format.
    static void write_ground_truth(const GeneratedDataset& dataset,
                                   const std::string& path);

private:
    IGHLocusConfig config_;
    std::mt19937_64 rng_;

    std::string generate_random_sequence(std::uint64_t length);
    char random_base();
    char mutate_base(char base);
    std::string apply_errors(const std::string& seq, float error_rate);
    std::string generate_quality_string(std::uint64_t length, float mean_error_rate);
};

// Preset configurations for common test scenarios.
namespace presets {

// Pure canonical (no duplication).
[[nodiscard]] IGHLocusConfig canonical_only(std::uint64_t seed = 42);

// 50/50 balanced mosaic.
[[nodiscard]] IGHLocusConfig balanced_mosaic(std::uint64_t seed = 42);

// Various dup fractions for sweep.
[[nodiscard]] IGHLocusConfig dup_fraction_5(std::uint64_t seed = 42);
[[nodiscard]] IGHLocusConfig dup_fraction_10(std::uint64_t seed = 42);
[[nodiscard]] IGHLocusConfig dup_fraction_30(std::uint64_t seed = 42);
[[nodiscard]] IGHLocusConfig dup_fraction_50(std::uint64_t seed = 42);
[[nodiscard]] IGHLocusConfig dup_fraction_100(std::uint64_t seed = 42);

// Edge case: high sequence-identical fraction.
[[nodiscard]] IGHLocusConfig seq_identical_stress(std::uint64_t seed = 42);

// Small fast test config.
[[nodiscard]] IGHLocusConfig tiny_test(std::uint64_t seed = 42);

}  // namespace presets

}  // namespace llmap::synthetic
