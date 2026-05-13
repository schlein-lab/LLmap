// LLmap — Parquet writer: AlignmentRecord to ProbabilityEntry conversion.

#include "output/parquet_writer.h"

namespace llmap::output {

std::vector<ProbabilityEntry> RecordToEntries(
    const AlignmentRecord& record,
    float min_probability) {

    std::vector<ProbabilityEntry> entries;

    // Get confidence (first score if available)
    float base_confidence = 0.0f;
    if (!record.confidence_scores.empty()) {
        base_confidence = record.confidence_scores[0];
    }

    if (record.status == AlignmentStatus::Mapped && record.primary) {
        // Primary alignment gets probability 1.0 (collapsed)
        ProbabilityEntry primary_entry{
            .read_id = record.read_id,
            .bucket_id = record.primary->target_id,
            .probability = 1.0f,
            .confidence = base_confidence,
            .level = record.collapsed_at_level,
            .iteration = record.collapsed_at_iteration,
            .is_collapsed = true,
        };
        entries.push_back(std::move(primary_entry));

        // Alternatives get lower probabilities
        for (std::size_t i = 0; i < record.alternatives.size(); ++i) {
            const auto& alt = record.alternatives[i];
            float alt_prob = 0.9f / static_cast<float>(i + 2);  // Decreasing
            if (alt_prob >= min_probability) {
                ProbabilityEntry alt_entry{
                    .read_id = record.read_id,
                    .bucket_id = alt.target_id,
                    .probability = alt_prob,
                    .confidence = base_confidence * 0.5f,  // Lower confidence
                    .level = record.collapsed_at_level,
                    .iteration = record.collapsed_at_iteration,
                    .is_collapsed = false,
                };
                entries.push_back(std::move(alt_entry));
            }
        }
    } else if (record.status == AlignmentStatus::Tentative) {
        // Tentative targets have explicit probabilities
        for (const auto& target : record.tentative_targets) {
            if (target.final_probability >= min_probability) {
                ProbabilityEntry entry{
                    .read_id = record.read_id,
                    .bucket_id = target.target_id,
                    .probability = target.final_probability,
                    .confidence = target.sequence_identity_estimate,
                    .level = record.collapsed_at_level,
                    .iteration = record.collapsed_at_iteration,
                    .is_collapsed = false,
                };
                entries.push_back(std::move(entry));
            }
        }
    } else {
        // Unmapped: single entry with probability 0
        ProbabilityEntry unmapped_entry{
            .read_id = record.read_id,
            .bucket_id = "*",  // No bucket
            .probability = 0.0f,
            .confidence = 0.0f,
            .level = 0,
            .iteration = 0,
            .is_collapsed = false,
        };
        entries.push_back(std::move(unmapped_entry));
    }

    return entries;
}

}  // namespace llmap::output
