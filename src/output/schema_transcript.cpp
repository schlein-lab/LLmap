// LLmap — Transcript-Mode output schema implementation.
//
// Each tag builder is a tiny pure helper. The Parquet sidecar row
// structs are header-only POD; no separate emit logic here. The
// parquet writer that consumes them lands in src/output/
// parquet_writer_transcript.cpp once arrow is wired in Block 9.

#include "output/schema_transcript.h"

#include <cstdio>
#include <sstream>

namespace llmap::output::transcript_schema {

namespace {

template <typename T>
std::string Join(std::span<const T> values, char sep = ',') {
    std::ostringstream os;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) os << sep;
        os << values[i];
    }
    return os.str();
}

}  // namespace

char XsTag(char strand) noexcept {
    return (strand == '+' || strand == '-') ? strand : '?';
}

std::string JiTag(
    std::span<const std::pair<std::uint64_t, std::uint64_t>> junctions) {
    if (junctions.empty()) return {};
    std::ostringstream os;
    for (std::size_t i = 0; i < junctions.size(); ++i) {
        if (i > 0) os << ',';
        os << junctions[i].first << ',' << junctions[i].second;
    }
    return os.str();
}

std::string JmTag(std::span<const float> per_junction_conf) {
    if (per_junction_conf.empty()) return {};
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(2);
    for (std::size_t i = 0; i < per_junction_conf.size(); ++i) {
        if (i > 0) os << ',';
        os << per_junction_conf[i];
    }
    return os.str();
}

const char* XkTag(AlignmentStatus s) noexcept {
    return AlignmentStatusName(s);
}

std::string XcTag(std::uint32_t cluster_id) {
    if (cluster_id == 0) return "-";
    return std::to_string(cluster_id);
}

std::string XaTag(std::span<const std::string> sources) {
    return Join(sources, ',');
}

std::string XqTag(std::int32_t legacy_mapq) {
    if (legacy_mapq < 0) return "-1";
    return std::to_string(legacy_mapq);
}

std::string XmTag(std::span<const ModCallView> calls) {
    if (calls.empty()) return {};
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(2);
    for (std::size_t i = 0; i < calls.size(); ++i) {
        if (i > 0) os << ',';
        os << calls[i].kind << ':' << calls[i].pos << ':' << calls[i].confidence;
    }
    return os.str();
}

std::string XfTag(std::string_view splicing_state_name) {
    return std::string(splicing_state_name);
}

}  // namespace llmap::output::transcript_schema
