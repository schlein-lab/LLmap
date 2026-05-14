// LLmap — Benchmark-compatible ground truth generation for Phase 11.
//
// Generates synthetic datasets in formats compatible with benchmarks/metrics/compute.py:
// - T1 (WGS): read_id<TAB>chrom<TAB>pos
// - T2 (paralog stress): read_id<TAB>paralog<TAB>chrom<TAB>pos

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace llmap::synthetic {

// Forward declarations.
struct SyntheticRead;
struct GeneratedDataset;

// Task types for benchmark generation.
enum class BenchmarkTask {
    T1_WGS,           // Synthetic-truth WGS: standard positional truth
    T2_ParalogStress  // Paralog stress: per-read paralog assignment truth
};

// Configuration for T1 (WGS) benchmark generation.
struct T1Config {
    std::uint64_t seed{42};
    std::uint64_t n_reads{1000000};
    std::uint32_t coverage{30};
    std::uint32_t read_length{10000};
    float read_length_stddev{1000.0f};
    float error_rate{0.001f};

    // Reference regions to generate reads from.
    // For T1: chr20 + chr14 (IGH locus)
    std::vector<std::string> chromosomes{"chr20", "chr14"};
    std::uint64_t region_length{5000000};  // per-chromosome length
};

// Paralog family definition for T2.
struct ParalogFamily {
    std::string name;                        // e.g. "IGHG"
    std::vector<std::string> members;        // e.g. {"IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP"}
    std::uint32_t region_length{50000};      // total region length
    std::uint32_t n_psvs{20};                // PSVs discriminating members
    float seq_identity{0.95f};               // member sequence identity
};

// Configuration for T2 (paralog stress) benchmark generation.
struct T2Config {
    std::uint64_t seed{42};
    std::uint64_t n_reads{500000};
    std::uint32_t coverage{30};
    std::uint32_t read_length{10000};
    float read_length_stddev{1000.0f};
    float error_rate{0.001f};

    // Paralog families to include.
    std::vector<ParalogFamily> families;
};

// Benchmark-compatible read with positional ground truth.
struct BenchmarkRead {
    std::string id;
    std::string sequence;
    std::string quality;

    // Ground truth for compute.py compatibility.
    std::string chrom;              // chromosome/contig name
    std::int64_t pos;               // 0-based position
    std::string paralog;            // T2 only: assigned paralog member
    bool is_forward{true};          // strand
};

// Dataset generated for benchmarks.
struct BenchmarkDataset {
    BenchmarkTask task;
    std::vector<BenchmarkRead> reads;

    // Reference sequences (name -> sequence).
    std::vector<std::pair<std::string, std::string>> references;

    // Summary statistics.
    std::uint64_t total_reads{0};
    std::uint64_t total_bases{0};
};

// Main benchmark generator class.
class BenchmarkGenerator {
public:
    // Generate T1 (WGS) dataset.
    [[nodiscard]] static BenchmarkDataset generate_t1(const T1Config& config);

    // Generate T2 (paralog stress) dataset.
    [[nodiscard]] static BenchmarkDataset generate_t2(const T2Config& config);

    // Write reads to FASTQ.
    static void write_fastq(const BenchmarkDataset& dataset, const std::string& path);

    // Write reference to FASTA.
    static void write_reference(const BenchmarkDataset& dataset, const std::string& path);

    // Write ground truth TSV compatible with compute.py.
    // T1 format: read_id<TAB>chrom<TAB>pos
    static void write_truth_tsv(const BenchmarkDataset& dataset, const std::string& path);

    // Write paralog truth TSV for T2.
    // Format: read_id<TAB>paralog<TAB>chrom<TAB>pos
    static void write_paralog_truth_tsv(const BenchmarkDataset& dataset,
                                        const std::string& path);
};

// Preset T2 paralog families.
namespace paralog_presets {

// IGH constant region paralogs.
[[nodiscard]] ParalogFamily igh_constant();

// NPHP1 / NPHP1B duplication.
[[nodiscard]] ParalogFamily nphp1();

// MHC class I pseudogenes (simplified).
[[nodiscard]] ParalogFamily mhc_class1();

}  // namespace paralog_presets

}  // namespace llmap::synthetic
