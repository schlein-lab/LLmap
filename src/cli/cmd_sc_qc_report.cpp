// LLmap — `llmap sc-qc-report` CLI command.
//
// Generates QC reports for single-cell paralog matrices with configurable
// thresholds for filtering and quality assessment.

#include "cli/commands.h"
#include "cli/cmd_sc_qc_report_internal.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <regex>
#include <string>
#include <unordered_map>

#include "core/alignment_record.h"
#include "output/parquet_reader.h"
#include "singlecell/per_cell_paralog.h"
#include "singlecell/sc_paralog_qc.h"

namespace llmap::cli {

using namespace sc_qc_internal;

int run_sc_qc_report(int argc, char** argv) {
    ScQcReportArgs args;
    if (!ParseScQcReportArgs(argc, argv, args)) {
        PrintScQcReportUsage();
        return 1;
    }

    if (args.help) {
        PrintScQcReportUsage();
        return 0;
    }

    if (args.parquet_input.empty()) {
        std::fprintf(stderr, "Error: --parquet is required\n");
        PrintScQcReportUsage();
        return 1;
    }

    if (args.qc_json.empty() && args.qc_tsv_dir.empty() &&
        args.filtered_matrix.empty()) {
        std::fprintf(stderr,
            "Error: at least one output is required "
            "(--qc-json, --qc-tsv, or --filtered-matrix)\n");
        PrintScQcReportUsage();
        return 1;
    }

    if (!std::filesystem::exists(args.parquet_input)) {
        std::fprintf(stderr, "Error: input file not found: %s\n",
                     args.parquet_input.c_str());
        return 1;
    }

    auto start_time = std::chrono::steady_clock::now();

    if (args.verbose) {
        std::fprintf(stderr, "Reading probability entries: %s\n",
                     args.parquet_input.c_str());
    }

    output::ParquetReaderConfig reader_cfg;
    reader_cfg.min_probability = args.min_prob;
    auto entries = output::ReadParquet(args.parquet_input, reader_cfg);

    if (entries.empty()) {
        std::fprintf(stderr, "Error: no entries read from input file\n");
        return 1;
    }

    if (args.verbose) {
        std::fprintf(stderr, "Read %zu probability entries\n", entries.size());
    }

    std::optional<std::regex> cb_regex;
    if (!args.cb_pattern.empty()) {
        try {
            cb_regex = std::regex(args.cb_pattern);
        } catch (const std::regex_error& e) {
            std::fprintf(stderr, "Error: invalid regex pattern: %s\n", e.what());
            return 1;
        }
    }

    std::unordered_map<std::string, std::string> cb_map;
    if (!args.cb_file.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Loading barcode map: %s\n", args.cb_file.c_str());
        }
        std::FILE* f = std::fopen(args.cb_file.c_str(), "r");
        if (!f) {
            std::fprintf(stderr, "Error: cannot open barcode file: %s\n",
                         args.cb_file.c_str());
            return 1;
        }
        char line[4096];
        while (std::fgets(line, sizeof(line), f)) {
            std::string_view sv(line);
            if (sv.empty() || sv[0] == '#') continue;
            auto tab = sv.find('\t');
            if (tab != std::string_view::npos) {
                auto read_id = std::string(sv.substr(0, tab));
                auto end = sv.find_first_of("\r\n", tab + 1);
                if (end == std::string_view::npos) end = sv.size();
                auto barcode = std::string(sv.substr(tab + 1, end - tab - 1));
                cb_map[read_id] = barcode;
            }
        }
        std::fclose(f);
        if (args.verbose) {
            std::fprintf(stderr, "Loaded %zu barcode mappings\n", cb_map.size());
        }
    }

    if (args.verbose) {
        std::fprintf(stderr, "Building cell x paralog matrix...\n");
    }

    singlecell::CellParalogMatrix matrix;
    std::size_t n_with_cb = 0;
    std::size_t n_skipped = 0;

    for (const auto& entry : entries) {
        std::optional<std::string> cell_barcode;

        if (!cb_map.empty()) {
            auto it = cb_map.find(entry.read_id);
            if (it != cb_map.end()) {
                cell_barcode = it->second;
            }
        } else if (cb_regex.has_value()) {
            cell_barcode = ExtractCbFromRegex(entry.read_id, *cb_regex);
        } else {
            cell_barcode = ExtractCbFromTag(entry.read_id, args.cb_tag);
        }

        if (!cell_barcode.has_value() || cell_barcode->empty()) {
            ++n_skipped;
            continue;
        }

        ++n_with_cb;

        AlignmentRecord record;
        record.read_id = entry.read_id;
        record.cell_barcode = std::move(cell_barcode);
        record.confidence_scores = {entry.confidence};

        ParalogCall paralog_call;
        paralog_call.inter_paralog.emplace_back(entry.bucket_id, entry.probability);
        record.paralog_assignment = std::move(paralog_call);

        matrix.AddRecord(record);
    }

    if (args.verbose) {
        std::fprintf(stderr, "Entries with cell barcode: %zu\n", n_with_cb);
        std::fprintf(stderr, "Entries skipped (no barcode): %zu\n", n_skipped);
    }

    singlecell::CellParalogConfig config;
    config.min_probability = args.min_prob;
    config.aggregation = ParseAggregation(args.aggregation);
    config.normalize_rows = false;

    matrix.Finalize(config);

    auto stats = matrix.GetStats();

    if (args.verbose) {
        std::fprintf(stderr, "Matrix statistics:\n");
        std::fprintf(stderr, "  Unique cells:    %zu\n",
                     static_cast<size_t>(stats.unique_cells));
        std::fprintf(stderr, "  Unique paralogs: %zu\n",
                     static_cast<size_t>(stats.unique_paralogs));
        std::fprintf(stderr, "  Matrix entries:  %zu\n",
                     static_cast<size_t>(stats.matrix_entries));
    }

    if (args.verbose) {
        std::fprintf(stderr, "Generating QC report...\n");
    }

    auto report = singlecell::GenerateQcReport(
        matrix, args.thresholds, args.sample_id);

    bool success = true;

    if (!args.qc_json.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Writing JSON report: %s\n", args.qc_json.c_str());
        }
        if (!singlecell::ExportQcReportJson(report, args.qc_json)) {
            std::fprintf(stderr, "Error: failed to write JSON report\n");
            success = false;
        }
    }

    if (!args.qc_tsv_dir.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Writing TSV reports: %s/\n", args.qc_tsv_dir.c_str());
        }
        std::filesystem::create_directories(args.qc_tsv_dir);
        if (!singlecell::ExportQcReportTsv(report, args.qc_tsv_dir)) {
            std::fprintf(stderr, "Error: failed to write TSV reports\n");
            success = false;
        }
    }

    if (!args.filtered_matrix.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Computing QC-filtered matrix...\n");
        }

        auto passing_cells = singlecell::GetCellsPassingQc(
            report.cells, args.thresholds);

        if (args.verbose) {
            std::fprintf(stderr, "Cells passing QC: %zu / %zu\n",
                         passing_cells.size(), report.cells.size());
        }

        auto filtered = singlecell::FilterMatrixByQc(matrix, passing_cells);

        if (args.verbose) {
            std::fprintf(stderr, "Writing filtered matrix: %s\n",
                         args.filtered_matrix.c_str());
        }

        std::filesystem::path out_path(args.filtered_matrix);
        std::string ext = out_path.extension().string();

        bool write_ok = false;
        if (ext == ".tsv") {
            write_ok = singlecell::ExportToTSV(filtered, out_path);
        } else {
            write_ok = singlecell::ExportToCSV(filtered, out_path);
        }

        if (!write_ok) {
            std::fprintf(stderr, "Error: failed to write filtered matrix\n");
            success = false;
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    float total_time_ms = std::chrono::duration<float, std::milli>(
        end_time - start_time).count();

    std::printf("Single-cell QC report complete:\n");
    std::printf("  Input entries:       %zu\n", entries.size());
    std::printf("  With barcode:        %zu\n", n_with_cb);
    std::printf("  Unique cells:        %zu\n",
                static_cast<size_t>(stats.unique_cells));
    std::printf("  Cells passing QC:    %zu\n",
                static_cast<size_t>(report.global.cells_passing_qc));
    std::printf("  Global assignment:   %.1f%%\n",
                100.0f * report.global.global_assignment_rate);
    std::printf("  Mean entropy:        %.3f\n", report.global.mean_entropy);
    std::printf("  Total time:          %.2f s\n", total_time_ms / 1000.0f);

    if (!args.qc_json.empty()) {
        std::printf("  JSON output:         %s\n", args.qc_json.c_str());
    }
    if (!args.qc_tsv_dir.empty()) {
        std::printf("  TSV output:          %s/\n", args.qc_tsv_dir.c_str());
    }
    if (!args.filtered_matrix.empty()) {
        std::printf("  Filtered matrix:     %s\n", args.filtered_matrix.c_str());
    }

    return success ? 0 : 1;
}

}  // namespace llmap::cli
