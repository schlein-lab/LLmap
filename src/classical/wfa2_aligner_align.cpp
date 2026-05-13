#include "classical/wfa2_aligner.h"
#include "classical/wfa2_aligner_impl.h"

#include <chrono>

namespace llmap::classical {

WFA2Aligner::WFA2Aligner(const WFA2Config& config)
    : config_(config),
      impl_(std::make_unique<WFA2AlignerImpl>(config)) {}

WFA2Aligner::~WFA2Aligner() = default;

WFA2Aligner::WFA2Aligner(WFA2Aligner&&) noexcept = default;
WFA2Aligner& WFA2Aligner::operator=(WFA2Aligner&&) noexcept = default;

bool WFA2Aligner::IsNativeWFA2() const {
#ifdef LLMAP_HAS_WFA2
    return true;
#else
    return false;
#endif
}

std::optional<WFA2Result> WFA2Aligner::Align(
    std::string_view query,
    std::string_view reference) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef LLMAP_HAS_WFA2
    // Native WFA2-lib implementation would go here
    // For now, fall through to fallback
#endif

    // Fallback to Gotoh DP
    FallbackAligner fallback(config_);
    auto result = fallback.Align(query, reference);

    auto end = std::chrono::high_resolution_clock::now();
    float time_us = std::chrono::duration<float, std::micro>(end - start).count();

    if (result) {
        result->alignment_time_us = time_us;
    }

    return result;
}

std::optional<WFA2Result> WFA2Aligner::ExtendRight(
    std::string_view query,
    std::string_view reference,
    int32_t query_start,
    int32_t ref_start) const {

    if (query_start < 0 || static_cast<size_t>(query_start) >= query.size() ||
        ref_start < 0 || static_cast<size_t>(ref_start) >= reference.size()) {
        return std::nullopt;
    }

    auto sub_query = query.substr(query_start);
    auto sub_ref = reference.substr(ref_start);

    auto result = Align(sub_query, sub_ref);
    if (result) {
        result->query_start += query_start;
        result->query_end += query_start;
        result->ref_start += ref_start;
        result->ref_end += ref_start;
    }

    return result;
}

std::optional<WFA2Result> WFA2Aligner::ExtendLeft(
    std::string_view query,
    std::string_view reference,
    int32_t query_end,
    int32_t ref_end) const {

    if (query_end <= 0 || static_cast<size_t>(query_end) > query.size() ||
        ref_end <= 0 || static_cast<size_t>(ref_end) > reference.size()) {
        return std::nullopt;
    }

    // Reverse the sequences and align
    std::string rev_query(query.substr(0, query_end));
    std::string rev_ref(reference.substr(0, ref_end));
    std::reverse(rev_query.begin(), rev_query.end());
    std::reverse(rev_ref.begin(), rev_ref.end());

    auto result = Align(rev_query, rev_ref);
    if (result) {
        // Reverse the CIGAR
        std::reverse(result->cigar.begin(), result->cigar.end());

        // Adjust coordinates
        result->query_start = query_end - result->query_end;
        result->query_end = query_end - result->query_start;
        result->ref_start = ref_end - result->ref_end;
        result->ref_end = ref_end - result->ref_start;
    }

    return result;
}

std::vector<std::optional<WFA2Result>> WFA2Aligner::AlignBatch(
    std::span<const std::string_view> queries,
    std::span<const std::string_view> references) const {

    if (queries.size() != references.size()) {
        return {};
    }

    std::vector<std::optional<WFA2Result>> results;
    results.reserve(queries.size());

    for (size_t i = 0; i < queries.size(); ++i) {
        results.push_back(Align(queries[i], references[i]));
    }

    return results;
}

}  // namespace llmap::classical
