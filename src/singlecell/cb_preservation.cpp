// LLmap — Single-cell barcode and UMI preservation implementation.

#include "cb_preservation.h"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <sstream>
#include <unordered_set>

namespace llmap::singlecell {

std::string_view TagTypeToString(TagType type) noexcept {
    switch (type) {
        case TagType::CB: return "CB";
        case TagType::CR: return "CR";
        case TagType::CY: return "CY";
        case TagType::UB: return "UB";
        case TagType::UR: return "UR";
        case TagType::UY: return "UY";
        case TagType::RG: return "RG";
        case TagType::BC: return "BC";
        case TagType::QT: return "QT";
    }
    return "??";
}

std::optional<TagType> StringToTagType(std::string_view tag) noexcept {
    if (tag == "CB") return TagType::CB;
    if (tag == "CR") return TagType::CR;
    if (tag == "CY") return TagType::CY;
    if (tag == "UB") return TagType::UB;
    if (tag == "UR") return TagType::UR;
    if (tag == "UY") return TagType::UY;
    if (tag == "RG") return TagType::RG;
    if (tag == "BC") return TagType::BC;
    if (tag == "QT") return TagType::QT;
    return std::nullopt;
}

TagValue TagValue::FromString(std::string value) {
    TagValue tv;
    tv.string_value = std::move(value);
    tv.is_string = true;
    return tv;
}

TagValue TagValue::FromInt(std::int64_t value) {
    TagValue tv;
    tv.int_value = value;
    tv.is_string = false;
    return tv;
}

bool TagValue::operator==(const TagValue& other) const noexcept {
    if (is_string != other.is_string) return false;
    if (is_string) return string_value == other.string_value;
    return int_value == other.int_value;
}

bool TagValue::operator!=(const TagValue& other) const noexcept {
    return !(*this == other);
}

bool SingleCellTags::HasCellBarcode() const noexcept {
    return cell_barcode.has_value() || cell_barcode_raw.has_value();
}

bool SingleCellTags::HasUmi() const noexcept {
    return umi.has_value() || umi_raw.has_value();
}

bool SingleCellTags::HasAnyTags() const noexcept {
    return HasCellBarcode() || HasUmi() ||
           read_group.has_value() ||
           sample_barcode.has_value() ||
           !custom_tags.empty();
}

std::optional<std::string_view> SingleCellTags::GetEffectiveCellBarcode() const noexcept {
    if (cell_barcode.has_value()) {
        return std::string_view{*cell_barcode};
    }
    if (cell_barcode_raw.has_value()) {
        return std::string_view{*cell_barcode_raw};
    }
    return std::nullopt;
}

std::optional<std::string_view> SingleCellTags::GetEffectiveUmi() const noexcept {
    if (umi.has_value()) {
        return std::string_view{*umi};
    }
    if (umi_raw.has_value()) {
        return std::string_view{*umi_raw};
    }
    return std::nullopt;
}

void SingleCellTags::MergeFrom(const SingleCellTags& other) {
    if (other.cell_barcode) cell_barcode = other.cell_barcode;
    if (other.cell_barcode_raw) cell_barcode_raw = other.cell_barcode_raw;
    if (other.cell_barcode_quality) cell_barcode_quality = other.cell_barcode_quality;
    if (other.umi) umi = other.umi;
    if (other.umi_raw) umi_raw = other.umi_raw;
    if (other.umi_quality) umi_quality = other.umi_quality;
    if (other.read_group) read_group = other.read_group;
    if (other.sample_barcode) sample_barcode = other.sample_barcode;
    if (other.sample_barcode_quality) sample_barcode_quality = other.sample_barcode_quality;
    for (const auto& [key, val] : other.custom_tags) {
        custom_tags[key] = val;
    }
}

void SingleCellTags::Clear() noexcept {
    cell_barcode.reset();
    cell_barcode_raw.reset();
    cell_barcode_quality.reset();
    umi.reset();
    umi_raw.reset();
    umi_quality.reset();
    read_group.reset();
    sample_barcode.reset();
    sample_barcode_quality.reset();
    custom_tags.clear();
}

bool SingleCellTags::operator==(const SingleCellTags& other) const noexcept {
    return cell_barcode == other.cell_barcode &&
           cell_barcode_raw == other.cell_barcode_raw &&
           cell_barcode_quality == other.cell_barcode_quality &&
           umi == other.umi &&
           umi_raw == other.umi_raw &&
           umi_quality == other.umi_quality &&
           read_group == other.read_group &&
           sample_barcode == other.sample_barcode &&
           sample_barcode_quality == other.sample_barcode_quality &&
           custom_tags == other.custom_tags;
}

bool SingleCellTags::operator!=(const SingleCellTags& other) const noexcept {
    return !(*this == other);
}

namespace {

// Parse a single SAM tag field: "XX:T:VALUE"
std::tuple<std::string, char, std::string> ParseSamTag(std::string_view field) {
    if (field.size() < 5 || field[2] != ':' || field[4] != ':') {
        return {"", '\0', ""};
    }
    return {
        std::string(field.substr(0, 2)),
        field[3],
        std::string(field.substr(5))
    };
}

// Check if character is valid DNA base
bool IsValidDnaBase(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
           c == 'a' || c == 'c' || c == 'g' || c == 't';
}

// Check if character is valid DNA base or N
bool IsValidDnaBaseOrN(char c) {
    return IsValidDnaBase(c) || c == 'N' || c == 'n';
}

}  // namespace

SingleCellTags ExtractTagsFromAux(std::string_view aux_string) {
    SingleCellTags tags;

    size_t pos = 0;
    while (pos < aux_string.size()) {
        size_t end = aux_string.find('\t', pos);
        if (end == std::string_view::npos) {
            end = aux_string.size();
        }

        auto field = aux_string.substr(pos, end - pos);
        auto [tag, type, value] = ParseSamTag(field);

        if (!tag.empty()) {
            if (tag == "CB") {
                tags.cell_barcode = value;
            } else if (tag == "CR") {
                tags.cell_barcode_raw = value;
            } else if (tag == "CY") {
                tags.cell_barcode_quality = value;
            } else if (tag == "UB") {
                tags.umi = value;
            } else if (tag == "UR") {
                tags.umi_raw = value;
            } else if (tag == "UY") {
                tags.umi_quality = value;
            } else if (tag == "RG") {
                tags.read_group = value;
            } else if (tag == "BC") {
                tags.sample_barcode = value;
            } else if (tag == "QT") {
                tags.sample_barcode_quality = value;
            } else {
                if (type == 'Z' || type == 'A' || type == 'H') {
                    tags.custom_tags[tag] = TagValue::FromString(value);
                } else if (type == 'i' || type == 'I' || type == 's' || type == 'S' ||
                           type == 'c' || type == 'C') {
                    std::int64_t int_val = 0;
                    std::from_chars(value.data(), value.data() + value.size(), int_val);
                    tags.custom_tags[tag] = TagValue::FromInt(int_val);
                }
            }
        }

        pos = end + 1;
    }

    return tags;
}

SingleCellTags ExtractTagsFromPairs(
    const std::vector<std::pair<std::string, std::string>>& tag_pairs) {
    SingleCellTags tags;

    for (const auto& [tag, value] : tag_pairs) {
        if (tag == "CB") {
            tags.cell_barcode = value;
        } else if (tag == "CR") {
            tags.cell_barcode_raw = value;
        } else if (tag == "CY") {
            tags.cell_barcode_quality = value;
        } else if (tag == "UB") {
            tags.umi = value;
        } else if (tag == "UR") {
            tags.umi_raw = value;
        } else if (tag == "UY") {
            tags.umi_quality = value;
        } else if (tag == "RG") {
            tags.read_group = value;
        } else if (tag == "BC") {
            tags.sample_barcode = value;
        } else if (tag == "QT") {
            tags.sample_barcode_quality = value;
        } else {
            tags.custom_tags[tag] = TagValue::FromString(value);
        }
    }

    return tags;
}

SingleCellTags ExtractTagsFromReadName(
    std::string_view read_name,
    const BarcodeExtractionConfig& config) {
    SingleCellTags tags;

    // Look for underscore-separated fields: @readname_BARCODE_UMI
    // Also handle 10x style: @SAMPLE:BARCODE:UMI
    size_t last_underscore = read_name.rfind('_');
    if (last_underscore == std::string_view::npos) {
        // Try colon separator
        last_underscore = read_name.rfind(':');
    }

    if (last_underscore == std::string_view::npos) {
        return tags;
    }

    // Potential UMI is after last separator
    auto potential_umi = read_name.substr(last_underscore + 1);
    if (potential_umi.empty()) {
        return tags;
    }

    // Check if it looks like a UMI (all DNA bases)
    bool valid_umi = std::all_of(potential_umi.begin(), potential_umi.end(),
        [&config](char c) {
            return config.allow_n_in_barcode ? IsValidDnaBaseOrN(c) : IsValidDnaBase(c);
        });

    if (!valid_umi) {
        return tags;
    }

    // Check length constraints
    if (config.expected_umi_length > 0 &&
        potential_umi.size() != config.expected_umi_length &&
        !config.allow_partial_barcode) {
        return tags;
    }

    // Found UMI, now look for barcode
    auto prefix = read_name.substr(0, last_underscore);
    size_t second_last = prefix.rfind('_');
    if (second_last == std::string_view::npos) {
        second_last = prefix.rfind(':');
    }

    if (second_last != std::string_view::npos) {
        auto potential_cb = prefix.substr(second_last + 1);

        bool valid_cb = std::all_of(potential_cb.begin(), potential_cb.end(),
            [&config](char c) {
                return config.allow_n_in_barcode ? IsValidDnaBaseOrN(c) : IsValidDnaBase(c);
            });

        if (valid_cb) {
            if (config.expected_barcode_length == 0 ||
                potential_cb.size() == config.expected_barcode_length ||
                config.allow_partial_barcode) {
                tags.cell_barcode_raw = std::string(potential_cb);
            }
        }
    }

    tags.umi_raw = std::string(potential_umi);
    return tags;
}

std::string FormatTagsAsAux(const SingleCellTags& tags) {
    std::ostringstream oss;
    bool first = true;

    auto emit = [&](std::string_view tag, const std::string& value) {
        if (!first) oss << '\t';
        oss << tag << ":Z:" << value;
        first = false;
    };

    if (tags.cell_barcode) emit("CB", *tags.cell_barcode);
    if (tags.cell_barcode_raw) emit("CR", *tags.cell_barcode_raw);
    if (tags.cell_barcode_quality) emit("CY", *tags.cell_barcode_quality);
    if (tags.umi) emit("UB", *tags.umi);
    if (tags.umi_raw) emit("UR", *tags.umi_raw);
    if (tags.umi_quality) emit("UY", *tags.umi_quality);
    if (tags.read_group) emit("RG", *tags.read_group);
    if (tags.sample_barcode) emit("BC", *tags.sample_barcode);
    if (tags.sample_barcode_quality) emit("QT", *tags.sample_barcode_quality);

    for (const auto& [tag, value] : tags.custom_tags) {
        if (!first) oss << '\t';
        if (value.is_string) {
            oss << tag << ":Z:" << value.string_value;
        } else {
            oss << tag << ":i:" << value.int_value;
        }
        first = false;
    }

    return oss.str();
}

std::vector<std::pair<std::string, std::string>> FormatTagsAsPairs(
    const SingleCellTags& tags) {
    std::vector<std::pair<std::string, std::string>> pairs;

    if (tags.cell_barcode) pairs.emplace_back("CB", *tags.cell_barcode);
    if (tags.cell_barcode_raw) pairs.emplace_back("CR", *tags.cell_barcode_raw);
    if (tags.cell_barcode_quality) pairs.emplace_back("CY", *tags.cell_barcode_quality);
    if (tags.umi) pairs.emplace_back("UB", *tags.umi);
    if (tags.umi_raw) pairs.emplace_back("UR", *tags.umi_raw);
    if (tags.umi_quality) pairs.emplace_back("UY", *tags.umi_quality);
    if (tags.read_group) pairs.emplace_back("RG", *tags.read_group);
    if (tags.sample_barcode) pairs.emplace_back("BC", *tags.sample_barcode);
    if (tags.sample_barcode_quality) pairs.emplace_back("QT", *tags.sample_barcode_quality);

    for (const auto& [tag, value] : tags.custom_tags) {
        if (value.is_string) {
            pairs.emplace_back(tag, value.string_value);
        } else {
            pairs.emplace_back(tag, std::to_string(value.int_value));
        }
    }

    return pairs;
}

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
            if (!IsValidDnaBaseOrN(c)) return false;
        } else {
            if (!IsValidDnaBase(c)) return false;
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
            if (!IsValidDnaBaseOrN(c)) return false;
        } else {
            if (!IsValidDnaBase(c)) return false;
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
        // Handle tab/space separated (first column only)
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

    // Hamming distance for same-length barcodes
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
