// LLmap — lossless alignment record.
//
// The lossless guarantee is enforced at construction: every record must have
// status ∈ {Mapped, Tentative, Unmapped} and Unmapped implies a RejectionReason.
// `make_*` factory functions are the only sanctioned way to construct.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>
#include <stdexcept>

namespace llmap {

// --- enums ---------------------------------------------------------------

enum class AlignmentStatus : std::uint8_t {
    Mapped,
    Tentative,
    Unmapped,
};

enum class RejectionReason : std::uint8_t {
    NoSeeds,
    LowSeedDensity,
    ChainScoreBelowThreshold,
    AmbiguousNoAnchor,
    FailedWfa2Extension,
    DidNotConverge,
    HostContamination,
    LowComplexity,
};

// --- per-hit records -----------------------------------------------------

struct CigarString {
    std::string ops;
};

struct PsvObservation {
    std::uint64_t psv_id;
    char allele;
    float quality;
};

struct AlignmentHit {
    std::string target_id;
    std::uint64_t start{0};
    std::uint64_t end{0};
    CigarString cigar;
    std::int32_t score{0};
    std::uint32_t nm{0};
    bool is_reverse{false};
    std::vector<PsvObservation> psv_calls;
};

struct TentativeTarget {
    std::string target_id;
    std::uint64_t approx_start{0};
    std::uint64_t approx_end{0};
    std::uint32_t n_seeds{0};
    std::int32_t partial_chain_score{0};
    float sequence_identity_estimate{0.0f};
    float final_probability{0.0f};
};

struct AnchorEvidence {
    bool five_prime_anchored{false};
    bool three_prime_anchored{false};
    float anchor_confidence{0.0f};
};

struct ParalogCall {
    // inter_paralog: P(paralog_k) keyed by paralog id; vector ordered for cache friendliness
    std::vector<std::pair<std::string, float>> inter_paralog;
    std::optional<float> p_canonical;
    std::optional<float> p_dup;
    std::uint32_t n_discriminating_psvs{0};
};

// --- the lossless record -------------------------------------------------

struct AlignmentRecord {
    std::string read_id;
    std::uint32_t read_len{0};
    AlignmentStatus status{AlignmentStatus::Unmapped};

    std::optional<AlignmentHit> primary;
    std::vector<AlignmentHit> alternatives;
    std::vector<TentativeTarget> tentative_targets;
    std::vector<float> confidence_scores;
    AnchorEvidence anchor_resolution;
    std::optional<ParalogCall> paralog_assignment;
    std::optional<std::string> cell_barcode;
    std::optional<std::uint8_t> haplotype;
    std::optional<RejectionReason> rejection_reason;

    // WaveCollapse provenance
    std::uint32_t collapsed_at_iteration{0};
    std::uint8_t collapsed_at_level{0};
    std::uint32_t cluster_id{0};
    bool is_cluster_representative{false};

    // Invariant: status implies which optional fields must be set.
    // Returns true iff the record is internally consistent.
    [[nodiscard]] bool is_lossless_consistent() const noexcept;
};

// --- factory functions (only sanctioned construction path) ---------------

[[nodiscard]] AlignmentRecord make_mapped(
    std::string read_id,
    std::uint32_t read_len,
    AlignmentHit primary,
    std::vector<AlignmentHit> alternatives = {}
);

[[nodiscard]] AlignmentRecord make_tentative(
    std::string read_id,
    std::uint32_t read_len,
    std::vector<TentativeTarget> targets,
    RejectionReason reason
);

[[nodiscard]] AlignmentRecord make_unmapped(
    std::string read_id,
    std::uint32_t read_len,
    RejectionReason reason
);

}  // namespace llmap
