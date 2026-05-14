// LLmap — cmd_sc_qc_report argument parsing and barcode extraction.

#include "cli/cmd_sc_qc_report_internal.h"

#include <cstdio>
#include <string>

namespace llmap::cli::sc_qc_internal {

void PrintScQcReportUsage() {
    std::puts(
        "Usage: llmap sc-qc-report [options]\n"
        "\n"
        "Generate QC reports for single-cell paralog assignment matrices.\n"
        "\n"
        "Required:\n"
        "  --parquet FILE        Input Parquet/CSV probability file\n"
        "\n"
        "Output (at least one required):\n"
        "  --qc-json FILE        Output QC report as JSON\n"
        "  --qc-tsv DIR          Output QC report as TSV files in directory\n"
        "  --filtered-matrix FILE  Output QC-filtered matrix (CSV/TSV)\n"
        "\n"
        "Cell barcode extraction (one of):\n"
        "  --cb-tag TAG          SAM-style tag in read name [CB:Z:]\n"
        "  --cb-pattern REGEX    Regex with capture group for barcode\n"
        "  --cb-file FILE        TSV with read_id<tab>cell_barcode mapping\n"
        "\n"
        "QC thresholds:\n"
        "  --min-assignment-rate FLOAT  Minimum assignment rate per cell [0.1]\n"
        "  --min-confidence FLOAT       Minimum mean confidence per cell [0.5]\n"
        "  --min-reads-per-cell INT     Minimum reads per cell [10]\n"
        "  --min-paralogs-per-cell INT  Minimum paralogs per cell [1]\n"
        "  --max-entropy FLOAT          Maximum assignment entropy [1.5]\n"
        "  --min-detection-rate FLOAT   Minimum paralog detection rate [0.01]\n"
        "\n"
        "Optional:\n"
        "  --sample-id STRING    Sample identifier for report\n"
        "  --min-prob FLOAT      Minimum probability threshold [0.01]\n"
        "  --aggregation MODE    mean, max, sum, weighted [mean]\n"
        "  --verbose             Print progress information\n"
        "  -h, --help            Show this help\n"
        "\n"
        "Example:\n"
        "  llmap sc-qc-report --parquet align.parquet --qc-json report.json\n"
        "  llmap sc-qc-report --parquet align.csv --qc-tsv qc_output/ \\\n"
        "    --min-reads-per-cell 50 --min-confidence 0.7\n"
        "  llmap sc-qc-report --parquet align.parquet \\\n"
        "    --qc-json report.json --filtered-matrix filtered.csv\n"
    );
}

bool ParseScQcReportArgs(int argc, char** argv, ScQcReportArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if (arg == "--parquet" && i + 1 < argc) {
            args.parquet_input = argv[++i];
        } else if (arg == "--qc-json" && i + 1 < argc) {
            args.qc_json = argv[++i];
        } else if (arg == "--qc-tsv" && i + 1 < argc) {
            args.qc_tsv_dir = argv[++i];
        } else if (arg == "--filtered-matrix" && i + 1 < argc) {
            args.filtered_matrix = argv[++i];
        } else if (arg == "--cb-tag" && i + 1 < argc) {
            args.cb_tag = argv[++i];
        } else if (arg == "--cb-pattern" && i + 1 < argc) {
            args.cb_pattern = argv[++i];
        } else if (arg == "--cb-file" && i + 1 < argc) {
            args.cb_file = argv[++i];
        } else if (arg == "--sample-id" && i + 1 < argc) {
            args.sample_id = argv[++i];
        } else if (arg == "--min-prob" && i + 1 < argc) {
            args.min_prob = std::stof(argv[++i]);
        } else if (arg == "--aggregation" && i + 1 < argc) {
            args.aggregation = argv[++i];
        } else if (arg == "--min-assignment-rate" && i + 1 < argc) {
            args.thresholds.min_assignment_rate = std::stof(argv[++i]);
        } else if (arg == "--min-confidence" && i + 1 < argc) {
            args.thresholds.min_confidence = std::stof(argv[++i]);
        } else if (arg == "--min-reads-per-cell" && i + 1 < argc) {
            args.thresholds.min_reads_per_cell =
                static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--min-paralogs-per-cell" && i + 1 < argc) {
            args.thresholds.min_paralogs_per_cell =
                static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--max-entropy" && i + 1 < argc) {
            args.thresholds.max_entropy = std::stof(argv[++i]);
        } else if (arg == "--min-detection-rate" && i + 1 < argc) {
            args.thresholds.min_detection_rate = std::stof(argv[++i]);
        } else if (arg == "--verbose" || arg == "-v") {
            args.verbose = true;
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

std::optional<std::string> ExtractCbFromTag(
    std::string_view read_id,
    std::string_view tag) {
    auto pos = read_id.find(tag);
    if (pos == std::string_view::npos) return std::nullopt;

    pos += tag.size();
    auto end = read_id.find_first_of(" \t_", pos);
    if (end == std::string_view::npos) end = read_id.size();

    return std::string(read_id.substr(pos, end - pos));
}

std::optional<std::string> ExtractCbFromRegex(
    std::string_view read_id,
    const std::regex& pattern) {
    std::string id_str(read_id);
    std::smatch match;
    if (std::regex_search(id_str, match, pattern) && match.size() > 1) {
        return match[1].str();
    }
    return std::nullopt;
}

singlecell::CellParalogConfig::AggregationMethod ParseAggregation(
    const std::string& method) {
    if (method == "mean") {
        return singlecell::CellParalogConfig::AggregationMethod::Mean;
    } else if (method == "max") {
        return singlecell::CellParalogConfig::AggregationMethod::Max;
    } else if (method == "sum") {
        return singlecell::CellParalogConfig::AggregationMethod::Sum;
    } else if (method == "weighted") {
        return singlecell::CellParalogConfig::AggregationMethod::Weighted;
    }
    return singlecell::CellParalogConfig::AggregationMethod::Mean;
}

}  // namespace llmap::cli::sc_qc_internal
