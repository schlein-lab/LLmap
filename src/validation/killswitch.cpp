// LLmap — Kill-switch validation framework implementation.

#include "validation/killswitch.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <sstream>
#include <unordered_set>

namespace llmap::validation {

namespace {

std::vector<std::string_view> SplitBy(std::string_view s, char delim) {
    std::vector<std::string_view> parts;
    size_t start = 0;
    while (start < s.size()) {
        auto pos = s.find(delim, start);
        if (pos == std::string_view::npos) {
            parts.push_back(s.substr(start));
            break;
        }
        parts.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return parts;
}

uint64_t ParseUint64(std::string_view s) {
    uint64_t val = 0;
    std::from_chars(s.data(), s.data() + s.size(), val);
    return val;
}

}  // namespace

std::optional<GroundTruth> GroundTruth::ParseFromReadId(std::string_view read_id) {
    auto parts = SplitBy(read_id, '_');
    if (parts.size() < 4) {
        return std::nullopt;
    }
    if (parts[0] != "synth") {
        return std::nullopt;
    }

    GroundTruth truth;
    truth.read_id = std::string(read_id);
    truth.origin = std::string(parts[2]);

    auto range_parts = SplitBy(parts[3], '-');
    if (range_parts.size() >= 2) {
        truth.true_start = ParseUint64(range_parts[0]);
        truth.true_end = ParseUint64(range_parts[1]);
    }

    if (parts.size() > 4) {
        truth.source_gene = std::string(parts[4]);
    }

    return truth;
}

KillSwitchValidator::KillSwitchValidator() = default;

KillSwitchValidator::KillSwitchValidator(KillSwitchThresholds thresholds)
    : thresholds_(std::move(thresholds)) {}

void KillSwitchValidator::LoadGroundTruth(const synthetic::GeneratedDataset& dataset) {
    ground_truth_.clear();
    for (const auto& read : dataset.reads) {
        GroundTruth truth;
        truth.read_id = read.id;
        truth.origin = (read.origin == synthetic::SyntheticRead::Origin::Duplicate)
                           ? "DUP" : "CAN";
        truth.true_start = read.true_start;
        truth.true_end = read.true_end;
        truth.source_gene = read.source_gene;
        truth.n_psvs_covered = read.covered_psv_positions.size();
        ground_truth_[truth.read_id] = std::move(truth);
    }
}

void KillSwitchValidator::LoadGroundTruth(const std::vector<GroundTruth>& truths) {
    ground_truth_.clear();
    for (const auto& t : truths) {
        ground_truth_[t.read_id] = t;
    }
}

void KillSwitchValidator::AddRecords(const std::vector<AlignmentRecord>& records) {
    for (const auto& r : records) {
        records_[r.read_id] = r;
    }
}

void KillSwitchValidator::AddRecord(const AlignmentRecord& record) {
    records_[record.read_id] = record;
}

bool KillSwitchValidator::CheckPositionCorrect(
    const AlignmentRecord& record,
    const GroundTruth& truth) const {

    if (!record.primary.has_value()) {
        return false;
    }

    const auto& hit = *record.primary;
    int64_t mapped_start = static_cast<int64_t>(hit.start);
    int64_t true_start = static_cast<int64_t>(truth.true_start);
    int64_t error = std::abs(mapped_start - true_start);

    int64_t read_len = static_cast<int64_t>(truth.true_end - truth.true_start);
    int64_t relative_tolerance = static_cast<int64_t>(
        read_len * position_tolerance_.relative_tolerance);
    int64_t tolerance = std::max(position_tolerance_.base_tolerance, relative_tolerance);

    return error <= tolerance;
}

ReadValidation KillSwitchValidator::ValidateRead(
    const std::string& read_id,
    const GroundTruth* truth,
    const AlignmentRecord* record) const {

    ReadValidation v;
    v.read_id = read_id;
    v.has_ground_truth = (truth != nullptr);
    v.has_alignment_record = (record != nullptr);

    if (!record) {
        return v;
    }

    v.is_lossless_consistent = record->is_lossless_consistent();
    v.mapped_status = record->status;

    if (!truth) {
        return v;
    }

    if (record->status == AlignmentStatus::Mapped && record->primary.has_value()) {
        const auto& hit = *record->primary;
        int64_t mapped_start = static_cast<int64_t>(hit.start);
        int64_t true_start = static_cast<int64_t>(truth->true_start);
        v.position_error = mapped_start - true_start;
        v.position_correct = CheckPositionCorrect(*record, *truth);
    }

    v.origin_known = !truth->origin.empty();
    if (v.origin_known && record->paralog_assignment.has_value()) {
        const auto& pa = *record->paralog_assignment;
        if (pa.p_canonical.has_value() && pa.p_dup.has_value()) {
            bool is_canonical = *pa.p_canonical > *pa.p_dup;
            std::string inferred_origin = is_canonical ? "CAN" : "DUP";
            v.origin_correct = (inferred_origin == truth->origin);
        }
    }

    return v;
}

std::vector<ReadValidation> KillSwitchValidator::ValidateDetailed() const {
    std::vector<ReadValidation> results;

    std::unordered_set<std::string> all_ids;
    for (const auto& [id, _] : ground_truth_) {
        all_ids.insert(id);
    }
    for (const auto& [id, _] : records_) {
        all_ids.insert(id);
    }

    for (const auto& id : all_ids) {
        const GroundTruth* truth = nullptr;
        const AlignmentRecord* record = nullptr;

        auto truth_it = ground_truth_.find(id);
        if (truth_it != ground_truth_.end()) {
            truth = &truth_it->second;
        }

        auto record_it = records_.find(id);
        if (record_it != records_.end()) {
            record = &record_it->second;
        }

        results.push_back(ValidateRead(id, truth, record));
    }

    return results;
}

ValidationStats KillSwitchValidator::Validate() const {
    ValidationStats stats;
    stats.input_reads = ground_truth_.size();
    stats.output_records = records_.size();

    auto detailed = ValidateDetailed();

    for (const auto& v : detailed) {
        if (v.has_ground_truth && !v.has_alignment_record) {
            ++stats.missing_records;
        }

        if (v.has_alignment_record) {
            switch (v.mapped_status) {
                case AlignmentStatus::Mapped:
                    ++stats.mapped;
                    break;
                case AlignmentStatus::Tentative:
                    ++stats.tentative;
                    break;
                case AlignmentStatus::Unmapped:
                    ++stats.unmapped;
                    break;
            }

            if (v.is_lossless_consistent) {
                ++stats.lossless_consistent;
            } else {
                ++stats.lossless_inconsistent;
            }
        }

        if (v.has_ground_truth && v.has_alignment_record &&
            v.mapped_status == AlignmentStatus::Mapped) {
            ++stats.position_evaluated;
            if (v.position_correct) {
                ++stats.position_correct;
            }
        }

        if (v.origin_known && v.has_alignment_record) {
            ++stats.origin_evaluated;
            if (v.origin_correct) {
                ++stats.origin_correct;
            }
        }
    }

    stats.Compute();
    return stats;
}

void ValidationStats::Compute() {
    lossless_pass = (missing_records == 0) && (lossless_inconsistent == 0);

    if (position_evaluated > 0) {
        position_accuracy = static_cast<double>(position_correct) /
                            static_cast<double>(position_evaluated);
    }

    if (origin_evaluated > 0) {
        origin_accuracy = static_cast<double>(origin_correct) /
                          static_cast<double>(origin_evaluated);
    }

    kill_switch_pass = lossless_pass;
    if (!lossless_pass) {
        kill_reason = "Lossless invariant violated";
    }
}

std::string ValidationStats::Summary() const {
    std::ostringstream ss;
    ss << "=== Kill-Switch Validation Summary ===\n\n";

    ss << "Lossless Guarantee:\n";
    ss << "  Input reads:     " << input_reads << "\n";
    ss << "  Output records:  " << output_records << "\n";
    ss << "  Missing:         " << missing_records << "\n";
    ss << "  Pass:            " << (lossless_pass ? "YES" : "NO") << "\n\n";

    ss << "Status Breakdown:\n";
    ss << "  Mapped:          " << mapped << "\n";
    ss << "  Tentative:       " << tentative << "\n";
    ss << "  Unmapped:        " << unmapped << "\n\n";

    ss << "Position Accuracy:\n";
    ss << "  Evaluated:       " << position_evaluated << "\n";
    ss << "  Correct:         " << position_correct << "\n";
    ss << "  Accuracy:        " << (position_accuracy * 100.0) << "%\n\n";

    ss << "Origin Accuracy (synthetic only):\n";
    ss << "  Evaluated:       " << origin_evaluated << "\n";
    ss << "  Correct:         " << origin_correct << "\n";
    ss << "  Accuracy:        " << (origin_accuracy * 100.0) << "%\n\n";

    ss << "Consistency:\n";
    ss << "  Consistent:      " << lossless_consistent << "\n";
    ss << "  Inconsistent:    " << lossless_inconsistent << "\n\n";

    ss << "Kill-Switch Verdict: " << (kill_switch_pass ? "PASS" : "FAIL");
    if (!kill_switch_pass) {
        ss << " (" << kill_reason << ")";
    }
    ss << "\n";

    return ss.str();
}

bool KillSwitchValidator::PassKillSwitch() const {
    auto stats = Validate();
    return stats.kill_switch_pass;
}

std::string KillSwitchValidator::KillSwitchVerdict() const {
    auto stats = Validate();
    return stats.Summary();
}

ValidationStats ValidateSyntheticRun(
    const synthetic::GeneratedDataset& dataset,
    const std::vector<AlignmentRecord>& records,
    const KillSwitchThresholds& thresholds) {

    KillSwitchValidator validator(thresholds);
    validator.LoadGroundTruth(dataset);
    validator.AddRecords(records);
    return validator.Validate();
}

BaselineComparison CompareToBaseline(
    const std::vector<AlignmentRecord>& llmap_records,
    const std::unordered_map<std::string, bool>& minimap2_mapped) {

    BaselineComparison cmp;

    std::unordered_set<std::string> all_ids;
    for (const auto& r : llmap_records) {
        all_ids.insert(r.read_id);
    }
    for (const auto& [id, _] : minimap2_mapped) {
        all_ids.insert(id);
    }

    std::unordered_map<std::string, bool> llmap_mapped;
    for (const auto& r : llmap_records) {
        llmap_mapped[r.read_id] = (r.status == AlignmentStatus::Mapped);
    }

    for (const auto& id : all_ids) {
        bool llm = llmap_mapped.count(id) && llmap_mapped.at(id);
        bool mm2 = minimap2_mapped.count(id) && minimap2_mapped.at(id);

        if (llm) ++cmp.llmap_mapped;
        if (mm2) ++cmp.minimap2_mapped;
        if (llm && mm2) ++cmp.both_mapped;
        if (llm && !mm2) ++cmp.llmap_only;
        if (!llm && mm2) ++cmp.minimap2_only;
        if (!llm && !mm2) ++cmp.neither_mapped;
    }

    size_t total = all_ids.size();
    if (total > 0) {
        cmp.llmap_recall = static_cast<double>(cmp.llmap_mapped) / total;
        cmp.minimap2_recall = static_cast<double>(cmp.minimap2_mapped) / total;
    }
    if (cmp.minimap2_mapped > 0) {
        cmp.recall_ratio = static_cast<double>(cmp.llmap_mapped) /
                           static_cast<double>(cmp.minimap2_mapped);
    }

    return cmp;
}

}  // namespace llmap::validation
