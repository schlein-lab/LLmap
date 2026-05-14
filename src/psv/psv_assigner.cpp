// LLmap — PSV-based paralog assigner implementation.

#include "psv/psv_assigner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

namespace llmap::psv {

namespace {

constexpr float kLogZero = -1e9f;

float LogSumExp(const std::vector<float>& log_probs) {
    if (log_probs.empty()) {
        return kLogZero;
    }

    const float max_log = *std::max_element(log_probs.begin(), log_probs.end());
    if (max_log <= kLogZero) {
        return kLogZero;
    }

    float sum = 0.0f;
    for (const float lp : log_probs) {
        sum += std::exp(lp - max_log);
    }

    return max_log + std::log(sum);
}

}  // namespace

PsvAssigner::PsvAssigner(const PsvCatalog& catalog, PsvAssignmentConfig config)
    : catalog_(catalog), config_(std::move(config)) {}

PsvAssignmentResult PsvAssigner::Assign(
    std::string_view read_id,
    const std::vector<PsvObservation>& observations) const {

    PsvAssignmentResult result;
    result.read_id = std::string{read_id};
    result.total_psvs_observed = static_cast<std::uint32_t>(observations.size());

    stats_.reads_processed++;

    if (observations.empty()) {
        return result;
    }

    stats_.reads_with_psvs++;

    const auto paralogs = catalog_.GetParalogs();
    if (paralogs.empty()) {
        return result;
    }

    result.likelihoods.reserve(paralogs.size());
    for (const auto& paralog_id : paralogs) {
        ParalogLikelihood pl;
        pl.paralog_id = paralog_id;
        pl.log_likelihood = 0.0f;
        result.likelihoods.push_back(std::move(pl));
    }

    for (const auto& obs : observations) {
        if (obs.base_quality < config_.min_base_quality) {
            continue;
        }

        const auto* site = catalog_.GetSiteById(obs.psv_id);
        if (!site || site->informativeness < config_.min_informativeness) {
            continue;
        }

        result.informative_psvs++;
        stats_.informative_observations++;

        for (auto& pl : result.likelihoods) {
            const float ll = ComputeLogLikelihood(obs, *site, pl.paralog_id);
            pl.log_likelihood += ll;

            auto it = site->paralog_alleles.find(pl.paralog_id);
            if (it != site->paralog_alleles.end()) {
                if (obs.observed_allele == it->second) {
                    pl.supporting_psvs++;
                } else {
                    pl.conflicting_psvs++;
                }
            }
        }
    }

    stats_.total_psv_observations += observations.size();

    if (result.informative_psvs < config_.min_psvs) {
        return result;
    }

    ComputePosteriors(result.likelihoods);
    result.entropy = ComputeEntropy(result.likelihoods);

    auto best_it = std::max_element(
        result.likelihoods.begin(),
        result.likelihoods.end(),
        [](const auto& a, const auto& b) { return a.posterior < b.posterior; });

    if (best_it != result.likelihoods.end()) {
        result.best_paralog = best_it->paralog_id;
        result.best_posterior = best_it->posterior;
        result.is_confident =
            result.best_posterior >= config_.confidence_threshold &&
            result.entropy <= config_.max_entropy;
    }

    if (result.is_confident) {
        stats_.reads_assigned++;
    } else if (result.informative_psvs >= config_.min_psvs) {
        stats_.reads_ambiguous++;
    }

    return result;
}

std::vector<PsvObservation> PsvAssigner::ExtractObservations(
    const AlignmentHit& hit,
    std::string_view read_sequence) const {

    std::vector<PsvObservation> observations;

    const auto sites = catalog_.GetSitesInRegion(hit.target_id, hit.start, hit.end);
    if (sites.empty()) {
        return observations;
    }

    std::uint64_t ref_pos = hit.start;
    std::uint64_t read_pos = 0;
    const auto& cigar = hit.cigar.ops;

    std::size_t i = 0;
    while (i < cigar.size()) {
        std::uint32_t len = 0;
        while (i < cigar.size() && std::isdigit(cigar[i])) {
            len = len * 10 + (cigar[i] - '0');
            i++;
        }

        if (i >= cigar.size()) break;
        const char op = cigar[i++];

        switch (op) {
            case 'M':
            case '=':
            case 'X':
                for (std::uint32_t j = 0; j < len && read_pos < read_sequence.size(); ++j) {
                    for (const auto* site : sites) {
                        if (site->position == ref_pos) {
                            PsvObservation obs;
                            obs.psv_id = site->psv_id;
                            obs.observed_allele = read_sequence[read_pos];
                            obs.read_position = static_cast<std::uint32_t>(read_pos);
                            obs.base_quality = 30;
                            observations.push_back(obs);
                        }
                    }
                    ref_pos++;
                    read_pos++;
                }
                break;
            case 'I':
            case 'S':
                read_pos += len;
                break;
            case 'D':
            case 'N':
                ref_pos += len;
                break;
            case 'H':
            case 'P':
                break;
            default:
                break;
        }
    }

    return observations;
}

void PsvAssigner::UpdateRecord(
    AlignmentRecord& record,
    std::string_view read_sequence) const {

    if (!record.primary.has_value()) {
        return;
    }

    auto observations = ExtractObservations(*record.primary, read_sequence);

    for (const auto& alt : record.alternatives) {
        auto alt_obs = ExtractObservations(alt, read_sequence);
        observations.insert(observations.end(), alt_obs.begin(), alt_obs.end());
    }

    if (observations.empty()) {
        return;
    }

    auto result = Assign(record.read_id, observations);

    if (record.paralog_assignment.has_value()) {
        record.paralog_assignment = MergeAssignments(*record.paralog_assignment, result);
    } else {
        record.paralog_assignment = ResultToParalogCall(result);
    }
}

void PsvAssigner::UpdateRecords(
    std::vector<AlignmentRecord>& records,
    const std::vector<std::string>& sequences) const {

    const auto n = std::min(records.size(), sequences.size());
    for (std::size_t i = 0; i < n; ++i) {
        UpdateRecord(records[i], sequences[i]);
    }
}

void PsvAssigner::ResetStats() noexcept {
    stats_ = PsvStats{};
}

float PsvAssigner::ComputeLogLikelihood(
    const PsvObservation& obs,
    const PsvSite& site,
    std::string_view paralog_id) const {

    const auto it = site.paralog_alleles.find(std::string{paralog_id});
    if (it == site.paralog_alleles.end()) {
        return 0.0f;
    }

    const char expected = it->second;
    const char observed = obs.observed_allele;

    if (observed == expected) {
        return std::log(1.0f - config_.error_rate);
    } else {
        return std::log(config_.error_rate / 3.0f);
    }
}

void PsvAssigner::ComputePosteriors(
    std::vector<ParalogLikelihood>& likelihoods) const {

    if (likelihoods.empty()) {
        return;
    }

    std::vector<float> log_probs;
    log_probs.reserve(likelihoods.size());
    for (const auto& pl : likelihoods) {
        log_probs.push_back(pl.log_likelihood);
    }

    const float log_sum = LogSumExp(log_probs);

    for (auto& pl : likelihoods) {
        pl.posterior = std::exp(pl.log_likelihood - log_sum);
    }
}

float PsvAssigner::ComputeEntropy(
    const std::vector<ParalogLikelihood>& likelihoods) const {

    float entropy = 0.0f;
    for (const auto& pl : likelihoods) {
        if (pl.posterior > 0.0f) {
            entropy -= pl.posterior * std::log2(pl.posterior);
        }
    }
    return entropy;
}

ParalogCall ResultToParalogCall(const PsvAssignmentResult& result) {
    ParalogCall call;
    call.n_discriminating_psvs = result.informative_psvs;

    for (const auto& pl : result.likelihoods) {
        call.inter_paralog.emplace_back(pl.paralog_id, pl.posterior);
    }

    std::sort(call.inter_paralog.begin(), call.inter_paralog.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    return call;
}

ParalogCall MergeAssignments(
    const ParalogCall& existing,
    const PsvAssignmentResult& psv_result,
    float psv_weight) {

    ParalogCall merged;
    merged.n_discriminating_psvs =
        existing.n_discriminating_psvs + psv_result.informative_psvs;

    std::unordered_map<std::string, float> existing_probs;
    for (const auto& [paralog, prob] : existing.inter_paralog) {
        existing_probs[paralog] = prob;
    }

    std::unordered_map<std::string, float> psv_probs;
    for (const auto& pl : psv_result.likelihoods) {
        psv_probs[pl.paralog_id] = pl.posterior;
    }

    std::set<std::string> all_paralogs;
    for (const auto& [p, _] : existing_probs) all_paralogs.insert(p);
    for (const auto& [p, _] : psv_probs) all_paralogs.insert(p);

    float sum = 0.0f;
    for (const auto& paralog : all_paralogs) {
        float existing_p = 0.0f;
        float psv_p = 0.0f;

        auto it1 = existing_probs.find(paralog);
        if (it1 != existing_probs.end()) {
            existing_p = it1->second;
        }

        auto it2 = psv_probs.find(paralog);
        if (it2 != psv_probs.end()) {
            psv_p = it2->second;
        }

        float combined = (1.0f - psv_weight) * existing_p + psv_weight * psv_p;

        if (psv_result.informative_psvs > 0 && existing.n_discriminating_psvs > 0) {
            combined = existing_p * psv_p;
        }

        merged.inter_paralog.emplace_back(paralog, combined);
        sum += combined;
    }

    if (sum > 0.0f) {
        for (auto& [_, prob] : merged.inter_paralog) {
            prob /= sum;
        }
    }

    std::sort(merged.inter_paralog.begin(), merged.inter_paralog.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    return merged;
}

}  // namespace llmap::psv
