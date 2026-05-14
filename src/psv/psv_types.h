// LLmap — PSV (Paralog-Specific Variant) types.
//
// PSVs are genomic positions where different paralogs have distinct alleles,
// enabling read-level disambiguation. These types define the PSV catalog
// and observation structures.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::psv {

// A PSV site in the reference
struct PsvSite {
    std::uint64_t psv_id{0};
    std::string chrom;
    std::uint64_t position{0};
    char ref_allele{'N'};

    // paralog_id -> expected allele at this position
    std::unordered_map<std::string, char> paralog_alleles;

    // Informativeness score: how well this site discriminates paralogs
    // Range [0,1]: 1.0 = perfect discriminator (each paralog has unique allele)
    float informativeness{0.0f};
};

// Observed allele at a PSV site from a read
struct PsvObservation {
    std::uint64_t psv_id{0};
    char observed_allele{'N'};
    std::uint8_t base_quality{0};
    std::uint32_t read_position{0};
};

// Likelihood of a read originating from a specific paralog
struct ParalogLikelihood {
    std::string paralog_id;
    float log_likelihood{0.0f};
    float posterior{0.0f};
    std::uint32_t supporting_psvs{0};
    std::uint32_t conflicting_psvs{0};
};

// Result of PSV-based paralog assignment for a single read
struct PsvAssignmentResult {
    std::string read_id;
    std::vector<ParalogLikelihood> likelihoods;
    std::string best_paralog;
    float best_posterior{0.0f};
    float entropy{0.0f};
    std::uint32_t total_psvs_observed{0};
    std::uint32_t informative_psvs{0};
    bool is_confident{false};
};

// Configuration for PSV-based assignment
struct PsvAssignmentConfig {
    // Minimum base quality to consider a PSV observation
    std::uint8_t min_base_quality{20};

    // Minimum informativeness score to use a PSV site
    float min_informativeness{0.5f};

    // Posterior threshold for confident assignment
    float confidence_threshold{0.9f};

    // Prior probability of sequencing error
    float error_rate{0.01f};

    // Maximum entropy for confident call (bits)
    float max_entropy{0.5f};

    // Minimum PSVs required for assignment
    std::uint32_t min_psvs{1};
};

// Statistics for PSV assignment
struct PsvStats {
    std::uint64_t reads_processed{0};
    std::uint64_t reads_with_psvs{0};
    std::uint64_t reads_assigned{0};
    std::uint64_t reads_ambiguous{0};
    std::uint64_t total_psv_observations{0};
    std::uint64_t informative_observations{0};
    float mean_psvs_per_read{0.0f};
    float assignment_rate{0.0f};
};

// Convert allele char to numeric (A=0, C=1, G=2, T=3, N=4)
[[nodiscard]] inline std::uint8_t AlleleToNum(char allele) noexcept {
    switch (allele) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// Convert numeric to allele char
[[nodiscard]] inline char NumToAllele(std::uint8_t num) noexcept {
    constexpr char alleles[] = {'A', 'C', 'G', 'T', 'N'};
    return num < 5 ? alleles[num] : 'N';
}

// Check if allele is valid (A/C/G/T)
[[nodiscard]] inline bool IsValidAllele(char allele) noexcept {
    return AlleleToNum(allele) < 4;
}

}  // namespace llmap::psv
