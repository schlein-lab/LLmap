// LLmap — Single-cell barcode and UMI preservation: extraction functions.

#include "cb_preservation.h"
#include "cb_preservation_internal.h"

#include <algorithm>
#include <charconv>
#include <sstream>

namespace llmap::singlecell {

SingleCellTags ExtractTagsFromAux(std::string_view aux_string) {
    SingleCellTags tags;

    size_t pos = 0;
    while (pos < aux_string.size()) {
        size_t end = aux_string.find('\t', pos);
        if (end == std::string_view::npos) {
            end = aux_string.size();
        }

        auto field = aux_string.substr(pos, end - pos);
        auto [tag, type, value] = internal::ParseSamTag(field);

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
        last_underscore = read_name.rfind(':');
    }

    if (last_underscore == std::string_view::npos) {
        return tags;
    }

    auto potential_umi = read_name.substr(last_underscore + 1);
    if (potential_umi.empty()) {
        return tags;
    }

    bool valid_umi = std::all_of(potential_umi.begin(), potential_umi.end(),
        [&config](char c) {
            return config.allow_n_in_barcode ?
                   internal::IsValidDnaBaseOrN(c) : internal::IsValidDnaBase(c);
        });

    if (!valid_umi) {
        return tags;
    }

    if (config.expected_umi_length > 0 &&
        potential_umi.size() != config.expected_umi_length &&
        !config.allow_partial_barcode) {
        return tags;
    }

    auto prefix = read_name.substr(0, last_underscore);
    size_t second_last = prefix.rfind('_');
    if (second_last == std::string_view::npos) {
        second_last = prefix.rfind(':');
    }

    if (second_last != std::string_view::npos) {
        auto potential_cb = prefix.substr(second_last + 1);

        bool valid_cb = std::all_of(potential_cb.begin(), potential_cb.end(),
            [&config](char c) {
                return config.allow_n_in_barcode ?
                       internal::IsValidDnaBaseOrN(c) : internal::IsValidDnaBase(c);
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

}  // namespace llmap::singlecell
