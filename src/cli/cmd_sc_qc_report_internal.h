// LLmap — cmd_sc_qc_report internal header.
// Shared types and declarations for split cmd_sc_qc_report modules.

#pragma once

#include <optional>
#include <regex>
#include <string>

#include "singlecell/per_cell_paralog.h"
#include "singlecell/sc_paralog_qc.h"

namespace llmap::cli::sc_qc_internal {

struct ScQcReportArgs {
    std::string parquet_input;
    std::string qc_json;
    std::string qc_tsv_dir;
    std::string filtered_matrix;
    std::string cb_tag{"CB:Z:"};
    std::string cb_pattern;
    std::string cb_file;
    std::string sample_id;
    float min_prob{0.01f};
    std::string aggregation{"mean"};
    singlecell::QcThresholds thresholds;
    bool verbose{false};
    bool help{false};
};

void PrintScQcReportUsage();
bool ParseScQcReportArgs(int argc, char** argv, ScQcReportArgs& args);

std::optional<std::string> ExtractCbFromTag(
    std::string_view read_id,
    std::string_view tag);

std::optional<std::string> ExtractCbFromRegex(
    std::string_view read_id,
    const std::regex& pattern);

singlecell::CellParalogConfig::AggregationMethod ParseAggregation(
    const std::string& method);

}  // namespace llmap::cli::sc_qc_internal
