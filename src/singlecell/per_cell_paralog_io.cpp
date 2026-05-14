// LLmap — Per-cell paralog matrix: I/O (CSV/TSV/H5AD export/import).

#include "per_cell_paralog.h"

#include <fstream>
#include <sstream>
#include <charconv>

namespace llmap::singlecell {

namespace {

void WriteCSVLine(std::ostream& os, const std::vector<std::string>& fields, char sep) {
    bool first = true;
    for (const auto& field : fields) {
        if (!first) os << sep;
        bool needs_quote = field.find(sep) != std::string::npos ||
                           field.find('"') != std::string::npos ||
                           field.find('\n') != std::string::npos;
        if (needs_quote) {
            os << '"';
            for (char c : field) {
                if (c == '"') os << '"';
                os << c;
            }
            os << '"';
        } else {
            os << field;
        }
        first = false;
    }
    os << '\n';
}

std::vector<std::string> ParseCSVLine(const std::string& line, char sep) {
    std::vector<std::string> fields;
    std::string current;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (in_quotes) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    current += '"';
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                current += c;
            }
        } else {
            if (c == '"') {
                in_quotes = true;
            } else if (c == sep) {
                fields.push_back(std::move(current));
                current.clear();
            } else {
                current += c;
            }
        }
    }
    fields.push_back(std::move(current));
    return fields;
}

}  // namespace

bool ExportToCSV(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        return false;
    }

    WriteCSVLine(file, {"cell_barcode", "paralog_id", "probability",
                        "read_count", "mean_confidence"}, ',');

    for (const auto& entry : matrix.GetEntries()) {
        WriteCSVLine(file, {
            entry.cell_barcode,
            entry.paralog_id,
            std::to_string(entry.probability),
            std::to_string(entry.read_count),
            std::to_string(entry.mean_confidence)
        }, ',');
    }

    return file.good();
}

bool ExportToTSV(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        return false;
    }

    WriteCSVLine(file, {"cell_barcode", "paralog_id", "probability",
                        "read_count", "mean_confidence"}, '\t');

    for (const auto& entry : matrix.GetEntries()) {
        WriteCSVLine(file, {
            entry.cell_barcode,
            entry.paralog_id,
            std::to_string(entry.probability),
            std::to_string(entry.read_count),
            std::to_string(entry.mean_confidence)
        }, '\t');
    }

    return file.good();
}

bool ExportToDenseCSV(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        return false;
    }

    auto [dense, cells, paralogs] = matrix.GetDenseMatrix();

    std::vector<std::string> header = {"cell_barcode"};
    header.insert(header.end(), paralogs.begin(), paralogs.end());
    WriteCSVLine(file, header, ',');

    for (size_t i = 0; i < cells.size(); ++i) {
        std::vector<std::string> row = {cells[i]};
        for (size_t j = 0; j < paralogs.size(); ++j) {
            row.push_back(std::to_string(dense[i][j]));
        }
        WriteCSVLine(file, row, ',');
    }

    return file.good();
}

bool ExportToH5AD(
    const CellParalogMatrix& matrix,
    const std::filesystem::path& path) {
    // HDF5 support requires libhdf5 + HighFive or similar
    // For V1.0, we indicate HDF5 is not available and return false
    // The caller should fall back to CSV
    (void)matrix;
    (void)path;
    return false;
}

std::optional<CellParalogMatrix> ImportFromCSV(
    const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        return std::nullopt;
    }

    std::string line;
    if (!std::getline(file, line)) {
        return std::nullopt;
    }

    auto header = ParseCSVLine(line, ',');
    if (header.size() < 3) {
        return std::nullopt;
    }

    CellParalogMatrix matrix;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto fields = ParseCSVLine(line, ',');
        if (fields.size() < 3) continue;

        AlignmentRecord record;
        record.read_id = "imported";
        record.read_len = 100;
        record.status = AlignmentStatus::Mapped;
        record.cell_barcode = fields[0];

        float prob = 0.0f;
        std::from_chars(fields[2].data(), fields[2].data() + fields[2].size(), prob);

        ParalogCall call;
        call.inter_paralog.emplace_back(fields[1], prob);
        record.paralog_assignment = call;

        matrix.AddRecord(record);
    }

    CellParalogConfig config;
    config.min_probability = 0.0f;
    config.min_reads_per_cell = 1;
    config.min_reads_per_paralog = 1;
    matrix.Finalize(config);

    return matrix;
}

}  // namespace llmap::singlecell
