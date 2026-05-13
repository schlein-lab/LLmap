// LLmap — Collapse Check for WaveCollapse convergence detection.
//
// Determines when reads have converged (max P(b|r) >= tau_collapse) and
// marks them as collapsed. Tracks dropout statistics for monitoring.
//
// Dropout types:
//   - Converged: max probability >= tau, marked as collapsed
//   - Ambiguous: multiple buckets with similar probability, needs refinement
//   - Low-confidence: max probability too low after many iterations
//   - Orphaned: no candidates remaining (empty probability distribution)

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "core/wave_state.h"

namespace llmap {

// Configuration for collapse checking
struct CollapseCheckConfig {
    float tau_collapse = 0.99f;      // Threshold for collapse (max P >= tau)
    float tau_ambiguous = 0.05f;     // Gap between top-2 below this is ambiguous
    float tau_low_confidence = 0.5f; // Max P below this is low-confidence
    bool mark_ambiguous_as_dropout = false;  // Whether ambiguous reads are dropouts
    bool auto_collapse = true;       // Actually collapse reads, or just detect
};

// Dropout classification for a single read
enum class DropoutType : std::uint8_t {
    None = 0,           // Not a dropout, still active
    Converged = 1,      // Collapsed successfully
    Ambiguous = 2,      // Multiple near-equal candidates
    LowConfidence = 3,  // Max probability too low
    Orphaned = 4,       // No candidates
};

// Per-read collapse analysis result
struct ReadCollapseInfo {
    std::uint32_t read_idx{0};
    DropoutType dropout_type{DropoutType::None};
    float max_probability{0.0f};
    float second_probability{0.0f};  // For ambiguity detection
    std::uint32_t top_bucket{0};
    std::uint32_t second_bucket{0};
    float entropy{0.0f};  // Shannon entropy of distribution
};

// Aggregate statistics from collapse check
struct CollapseCheckStats {
    std::uint32_t total_reads{0};
    std::uint32_t already_collapsed{0};
    std::uint32_t newly_collapsed{0};
    std::uint32_t active_after{0};
    std::uint32_t ambiguous{0};
    std::uint32_t low_confidence{0};
    std::uint32_t orphaned{0};

    // Probability statistics
    float mean_max_prob{0.0f};
    float min_max_prob{1.0f};
    float max_max_prob{0.0f};
    float mean_entropy{0.0f};

    // Convergence quality
    float collapse_rate() const {
        return total_reads > 0
            ? static_cast<float>(already_collapsed + newly_collapsed) / total_reads
            : 0.0f;
    }

    float dropout_rate() const {
        std::uint32_t dropouts = ambiguous + low_confidence + orphaned;
        return total_reads > 0
            ? static_cast<float>(dropouts) / total_reads
            : 0.0f;
    }
};

// Collapse checker class
class CollapseChecker {
public:
    explicit CollapseChecker(const CollapseCheckConfig& config = {});
    ~CollapseChecker() = default;

    CollapseChecker(const CollapseChecker&) = default;
    CollapseChecker& operator=(const CollapseChecker&) = default;
    CollapseChecker(CollapseChecker&&) = default;
    CollapseChecker& operator=(CollapseChecker&&) = default;

    // Perform collapse check on all reads
    // If auto_collapse is enabled, collapses converged reads in-place
    CollapseCheckStats Check(WaveState& state);

    // Perform collapse check and return per-read details
    // Useful for detailed analysis and debugging
    CollapseCheckStats CheckDetailed(
        WaveState& state,
        std::vector<ReadCollapseInfo>& read_info);

    // Check a single read without modifying state
    ReadCollapseInfo AnalyzeRead(
        const WaveState& state,
        std::uint32_t read_idx) const;

    // Get/set configuration
    const CollapseCheckConfig& Config() const { return config_; }
    void SetConfig(const CollapseCheckConfig& config) { config_ = config; }

    // Set collapse threshold
    void SetTauCollapse(float tau) { config_.tau_collapse = tau; }

private:
    CollapseCheckConfig config_;

    // Classify dropout type for a read
    DropoutType ClassifyDropout(const ReadCollapseInfo& info) const;
};

// Convenience function: check and collapse with default config
CollapseCheckStats CheckAndCollapse(
    WaveState& state,
    float tau_collapse = 0.99f);

// Convenience function: get dropout classification for each read
std::vector<DropoutType> GetDropoutTypes(
    const WaveState& state,
    float tau_collapse = 0.99f,
    float tau_ambiguous = 0.05f);

// Compute Shannon entropy of probability distribution
float ComputeReadEntropy(std::span<const float> probs);

// Get indices of reads by dropout type
std::vector<std::uint32_t> GetReadsByDropoutType(
    const std::vector<DropoutType>& types,
    DropoutType target);

}  // namespace llmap
