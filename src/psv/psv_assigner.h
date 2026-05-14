// LLmap — PSV-based paralog assigner.
//
// Uses Paralog-Specific Variants to compute posterior probabilities for
// each paralog given observed alleles in a read. Integrates with
// AlignmentRecord to populate ParalogCall fields.

#pragma once

#include "core/alignment_record.h"
#include "psv/psv_catalog.h"
#include "psv/psv_types.h"

#include <cstdint>
#include <string>
#include <vector>

namespace llmap::psv {

// PSV-based paralog assignment engine
class PsvAssigner {
public:
    explicit PsvAssigner(const PsvCatalog& catalog,
                         PsvAssignmentConfig config = {});

    // Assign paralog probabilities for a single read
    // observations: PSV alleles observed at each PSV site
    [[nodiscard]] PsvAssignmentResult Assign(
        std::string_view read_id,
        const std::vector<PsvObservation>& observations) const;

    // Extract PSV observations from an alignment
    [[nodiscard]] std::vector<PsvObservation> ExtractObservations(
        const AlignmentHit& hit,
        std::string_view read_sequence) const;

    // Update AlignmentRecord with PSV-based paralog assignment
    void UpdateRecord(AlignmentRecord& record,
                      std::string_view read_sequence) const;

    // Batch update multiple records
    void UpdateRecords(std::vector<AlignmentRecord>& records,
                       const std::vector<std::string>& sequences) const;

    // Get assignment statistics
    [[nodiscard]] const PsvStats& GetStats() const noexcept { return stats_; }

    // Reset statistics
    void ResetStats() noexcept;

    // Get configuration
    [[nodiscard]] const PsvAssignmentConfig& GetConfig() const noexcept {
        return config_;
    }

private:
    const PsvCatalog& catalog_;
    PsvAssignmentConfig config_;
    mutable PsvStats stats_;

    // Compute log-likelihood of observation given paralog
    [[nodiscard]] float ComputeLogLikelihood(
        const PsvObservation& obs,
        const PsvSite& site,
        std::string_view paralog_id) const;

    // Compute posterior from log-likelihoods using log-sum-exp
    void ComputePosteriors(std::vector<ParalogLikelihood>& likelihoods) const;

    // Compute entropy of posterior distribution
    [[nodiscard]] float ComputeEntropy(
        const std::vector<ParalogLikelihood>& likelihoods) const;
};

// Convert PsvAssignmentResult to ParalogCall for AlignmentRecord
[[nodiscard]] ParalogCall ResultToParalogCall(const PsvAssignmentResult& result);

// Merge PSV-based assignment with existing probabilistic assignment
// Uses Bayesian combination: P(paralog|psv,prob) ∝ P(psv|paralog) * P(paralog|prob)
[[nodiscard]] ParalogCall MergeAssignments(
    const ParalogCall& existing,
    const PsvAssignmentResult& psv_result,
    float psv_weight = 0.5f);

}  // namespace llmap::psv
