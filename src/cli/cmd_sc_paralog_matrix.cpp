// LLmap — `llmap sc-paralog-matrix` CLI command.
//
// Reads parquet/CSV probability entries and builds a cell × paralog matrix.
// Outputs sparse CSV or dense CSV for downstream single-cell analysis.

#include "cli/commands.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <regex>
#include <string>

#include "core/alignment_record.h"
#include "output/parquet_reader.h"
#include "singlecell/per_cell_paralog.h"

namespace llmap::cli {

namespace {

void print_sc_paralog_matrix_usage() {
    std::puts(
        "Usage: llmap sc-paralog-matrix [options]\n"
        "\n"
        "Build cell × paralog probability matrix from alignment output.\n"
        "\n"
        "Required:\n"
        "  --parquet FILE      Input Parquet/CSV probability file\n"
        "  --output FILE       Output CSV file (or .h5ad for AnnData)\n"
        "\n"
        "Cell barcode extraction (one of):\n"
        "  --cb-tag TAG        SAM-style tag in read name, e.g. 'CB:Z:' [CB:Z:]\n"
        "  --cb-pattern REGEX  Regex with capture group for barcode\n"
        "  --cb-file FILE      TSV with read_id<tab>cell_barcode mapping\n"
        "\n"
        "Optional:\n"
        "  --min-prob FLOAT    Minimum probability threshold [0.01]\n"
        "  --min-reads INT     Minimum reads per cell [1]\n"
        "  --aggregation MODE  mean, max, sum, weighted [mean]\n"
        "  --normalize         Normalize rows to sum to 1.0\n"
        "  --dense             Output dense matrix instead of sparse\n"
        "  --verbose           Print progress information\n"
        "  -h, --help          Show this help\n"
        "\n"
        "Output formats:\n"
        "  .csv    Sparse CSV: cell_barcode,paralog_id,probability,read_count\n"
        "  .tsv    Sparse TSV\n"
        "  .h5ad   AnnData format (requires HDF5)\n"
        "  --dense Dense matrix: rows=cells, columns=paralogs\n"
        "\n"
        "Example:\n"
        "  llmap sc-paralog-matrix --parquet align.parquet --output matrix.csv\n"
        "  llmap sc-paralog-matrix --parquet align.csv --output matrix.h5ad \\\n"
        "    --cb-pattern '([ACGT]{16})_' --min-reads 10\n"
    );
}

struct ScParalogMatrixArgs {
    std::string parquet_input;
    std::string output;
    std::string cb_tag{"CB:Z:"};
    std::string cb_pattern;
    std::string cb_file;
    float min_prob{0.01f};
    std::uint32_t min_reads{1};
    std::string aggregation{"mean"};
    bool normalize{false};
    bool dense{false};
    bool verbose{false};
    bool help{false};
};

bool parse_args(int argc, char** argv, ScParalogMatrixArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if (arg == "--parquet" && i + 1 < argc) {
            args.parquet_input = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            args.output = argv[++i];
        } else if (arg == "--cb-tag" && i + 1 < argc) {
            args.cb_tag = argv[++i];
        } else if (arg == "--cb-pattern" && i + 1 < argc) {
            args.cb_pattern = argv[++i];
        } else if (arg == "--cb-file" && i + 1 < argc) {
            args.cb_file = argv[++i];
        } else if (arg == "--min-prob" && i + 1 < argc) {
            args.min_prob = std::stof(argv[++i]);
        } else if (arg == "--min-reads" && i + 1 < argc) {
            args.min_reads = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--aggregation" && i + 1 < argc) {
            args.aggregation = argv[++i];
        } else if (arg == "--normalize") {
            args.normalize = true;
        } else if (arg == "--dense") {
            args.dense = true;
        } else if (arg == "--verbose" || arg == "-v") {
            args.verbose = true;
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

std::optional<std::string> extract_cb_from_tag(
    std::string_view read_id,
    std::string_view tag) {
    auto pos = read_id.find(tag);
    if (pos == std::string_view::npos) return std::nullopt;

    pos += tag.size();
    auto end = read_id.find_first_of(" \t_", pos);
    if (end == std::string_view::npos) end = read_id.size();

    return std::string(read_id.substr(pos, end - pos));
}

std::optional<std::string> extract_cb_from_regex(
    std::string_view read_id,
    const std::regex& pattern) {
    std::string id_str(read_id);
    std::smatch match;
    if (std::regex_search(id_str, match, pattern) && match.size() > 1) {
        return match[1].str();
    }
    return std::nullopt;
}

singlecell::CellParalogConfig::AggregationMethod parse_aggregation(
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

}  // namespace

int run_sc_paralog_matrix(int argc, char** argv) {
    ScParalogMatrixArgs args;
    if (!parse_args(argc, argv, args)) {
        print_sc_paralog_matrix_usage();
        return 1;
    }

    if (args.help) {
        print_sc_paralog_matrix_usage();
        return 0;
    }

    if (args.parquet_input.empty()) {
        std::fprintf(stderr, "Error: --parquet is required\n");
        print_sc_paralog_matrix_usage();
        return 1;
    }

    if (args.output.empty()) {
        std::fprintf(stderr, "Error: --output is required\n");
        print_sc_paralog_matrix_usage();
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
        std::fprintf(stderr, "Building cell × paralog matrix...\n");
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
            cell_barcode = extract_cb_from_regex(entry.read_id, *cb_regex);
        } else {
            cell_barcode = extract_cb_from_tag(entry.read_id, args.cb_tag);
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
    config.min_reads_per_cell = args.min_reads;
    config.aggregation = parse_aggregation(args.aggregation);
    config.normalize_rows = args.normalize;

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
        std::fprintf(stderr, "  Sparsity:        %.2f%%\n", 100.0f * stats.sparsity);
    }

    if (stats.matrix_entries == 0) {
        std::fprintf(stderr, "Warning: no matrix entries after filtering\n");
    }

    if (args.verbose) {
        std::fprintf(stderr, "Writing output: %s\n", args.output.c_str());
    }

    bool success = false;
    std::filesystem::path out_path(args.output);
    std::string ext = out_path.extension().string();

    if (args.dense) {
        success = singlecell::ExportToDenseCSV(matrix, out_path);
    } else if (ext == ".tsv") {
        success = singlecell::ExportToTSV(matrix, out_path);
    } else if (ext == ".h5ad") {
        success = singlecell::ExportToH5AD(matrix, out_path);
        if (!success) {
            std::fprintf(stderr, "Warning: h5ad export failed (HDF5 not available?), "
                         "falling back to CSV\n");
            out_path.replace_extension(".csv");
            success = singlecell::ExportToCSV(matrix, out_path);
        }
    } else {
        success = singlecell::ExportToCSV(matrix, out_path);
    }

    if (!success) {
        std::fprintf(stderr, "Error: failed to write output file\n");
        return 1;
    }

    auto end_time = std::chrono::steady_clock::now();
    float total_time_ms = std::chrono::duration<float, std::milli>(
        end_time - start_time).count();

    std::printf("Cell × Paralog matrix complete:\n");
    std::printf("  Input entries:   %zu\n", entries.size());
    std::printf("  With barcode:    %zu\n", n_with_cb);
    std::printf("  Unique cells:    %zu\n", static_cast<size_t>(stats.unique_cells));
    std::printf("  Unique paralogs: %zu\n", static_cast<size_t>(stats.unique_paralogs));
    std::printf("  Matrix entries:  %zu\n", static_cast<size_t>(stats.matrix_entries));
    std::printf("  Sparsity:        %.2f%%\n", 100.0f * stats.sparsity);
    std::printf("  Total time:      %.2f s\n", total_time_ms / 1000.0f);
    std::printf("  Output:          %s\n", out_path.c_str());

    return 0;
}

}  // namespace llmap::cli
