// LLmap — Single-cell barcode and UMI preservation: validation and whitelist.

#include "cb_preservation.h"
#include "cb_preservation_internal.h"

#include <fstream>

namespace llmap::singlecell {

bool ValidateBarcode(
    std::string_view barcode,
    const BarcodeExtractionConfig& config) {
    if (barcode.empty()) return false;

    if (config.expected_barcode_length > 0 &&
        barcode.size() != config.expected_barcode_length &&
        !config.allow_partial_barcode) {
        return false;
    }

    for (char c : barcode) {
        if (config.allow_n_in_barcode) {
            if (!internal::IsValidDnaBaseOrN(c)) return false;
        } else {
            if (!internal::IsValidDnaBase(c)) return false;
        }
    }

    return true;
}

bool ValidateUmi(
    std::string_view umi,
    const BarcodeExtractionConfig& config) {
    if (umi.empty()) return false;

    if (config.expected_umi_length > 0 &&
        umi.size() != config.expected_umi_length &&
        !config.allow_partial_barcode) {
        return false;
    }

    for (char c : umi) {
        if (config.allow_n_in_barcode) {
            if (!internal::IsValidDnaBaseOrN(c)) return false;
        } else {
            if (!internal::IsValidDnaBase(c)) return false;
        }
    }

    return true;
}

std::optional<CellBarcodeWhitelist> CellBarcodeWhitelist::LoadFromFile(
    const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        return std::nullopt;
    }

    CellBarcodeWhitelist whitelist;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        size_t sep = line.find_first_of("\t ");
        if (sep != std::string::npos) {
            line = line.substr(0, sep);
        }
        whitelist.Add(std::move(line));
    }

    return whitelist;
}

void CellBarcodeWhitelist::Add(std::string barcode) {
    barcodes_.insert(std::move(barcode));
}

bool CellBarcodeWhitelist::Contains(std::string_view barcode) const {
    return barcodes_.contains(std::string(barcode));
}

std::optional<std::string> CellBarcodeWhitelist::FindNearest(
    std::string_view barcode,
    std::uint32_t max_edit_distance) const {
    if (Contains(barcode)) {
        return std::string(barcode);
    }

    if (max_edit_distance == 0) {
        return std::nullopt;
    }

    std::string best_match;
    std::uint32_t best_distance = max_edit_distance + 1;

    for (const auto& candidate : barcodes_) {
        if (candidate.size() != barcode.size()) continue;

        std::uint32_t distance = 0;
        for (size_t i = 0; i < barcode.size() && distance <= max_edit_distance; ++i) {
            if (candidate[i] != barcode[i]) {
                ++distance;
            }
        }

        if (distance < best_distance) {
            best_distance = distance;
            best_match = candidate;
        }
    }

    if (best_distance <= max_edit_distance) {
        return best_match;
    }

    return std::nullopt;
}

std::size_t CellBarcodeWhitelist::Size() const noexcept {
    return barcodes_.size();
}

bool CellBarcodeWhitelist::Empty() const noexcept {
    return barcodes_.empty();
}

}  // namespace llmap::singlecell
