// LLmap — Single-cell barcode and UMI preservation: core types implementation.

#include "cb_preservation.h"
#include "cb_preservation_internal.h"

#include <algorithm>

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

namespace internal {

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

bool IsValidDnaBase(char c) noexcept {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
           c == 'a' || c == 'c' || c == 'g' || c == 't';
}

bool IsValidDnaBaseOrN(char c) noexcept {
    return IsValidDnaBase(c) || c == 'N' || c == 'n';
}

}  // namespace internal

}  // namespace llmap::singlecell
