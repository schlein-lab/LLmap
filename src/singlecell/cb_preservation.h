// LLmap — Single-cell barcode and UMI preservation.
//
// Preserves cell barcodes (CB), corrected barcodes (CB), UMI tags (UB/UY/UR),
// and cell Ranger-style tags through the alignment pipeline. This ensures
// no information loss when processing single-cell data through LLmap.
//
// Supports:
// - 10x Genomics (CB/UB tags)
// - Parse Biosciences (CB/UB/UY tags)
// - PacBio Kinnex/MAS-Seq (CB/UB tags)
// - Generic SAM/BAM single-cell tags

#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace llmap::singlecell {

// Single-cell tag types following SAM specification
enum class TagType : std::uint8_t {
    CB,  // Cell Barcode (corrected)
    CR,  // Cell Barcode (raw/uncorrected)
    CY,  // Cell Barcode quality (Phred+33)
    UB,  // UMI (corrected)
    UR,  // UMI (raw/uncorrected)
    UY,  // UMI quality (Phred+33)
    RG,  // Read group
    BC,  // Sample barcode sequence
    QT,  // Sample barcode quality
};

// Convert TagType to SAM tag string
[[nodiscard]] std::string_view TagTypeToString(TagType type) noexcept;

// Parse SAM tag string to TagType (returns nullopt if not a recognized SC tag)
[[nodiscard]] std::optional<TagType> StringToTagType(std::string_view tag) noexcept;

// A single-cell tag value (string or int)
struct TagValue {
    std::string string_value;
    std::int64_t int_value{0};
    bool is_string{true};

    // Factory methods
    static TagValue FromString(std::string value);
    static TagValue FromInt(std::int64_t value);

    // Comparison
    bool operator==(const TagValue& other) const noexcept;
    bool operator!=(const TagValue& other) const noexcept;
};

// Container for all single-cell tags associated with a read
struct SingleCellTags {
    // Core cell barcode (CB tag - corrected barcode)
    std::optional<std::string> cell_barcode;

    // Raw/uncorrected cell barcode (CR tag)
    std::optional<std::string> cell_barcode_raw;

    // Cell barcode quality (CY tag)
    std::optional<std::string> cell_barcode_quality;

    // Unique Molecular Identifier - corrected (UB tag)
    std::optional<std::string> umi;

    // UMI raw/uncorrected (UR tag)
    std::optional<std::string> umi_raw;

    // UMI quality (UY tag)
    std::optional<std::string> umi_quality;

    // Read group (RG tag)
    std::optional<std::string> read_group;

    // Sample barcode (BC tag)
    std::optional<std::string> sample_barcode;

    // Sample barcode quality (QT tag)
    std::optional<std::string> sample_barcode_quality;

    // Additional custom tags (key -> value)
    std::unordered_map<std::string, TagValue> custom_tags;

    // Check if any single-cell tags are present
    [[nodiscard]] bool HasCellBarcode() const noexcept;
    [[nodiscard]] bool HasUmi() const noexcept;
    [[nodiscard]] bool HasAnyTags() const noexcept;

    // Get the effective cell barcode (CB if present, else CR)
    [[nodiscard]] std::optional<std::string_view> GetEffectiveCellBarcode() const noexcept;

    // Get the effective UMI (UB if present, else UR)
    [[nodiscard]] std::optional<std::string_view> GetEffectiveUmi() const noexcept;

    // Merge tags from another instance (other takes precedence for conflicts)
    void MergeFrom(const SingleCellTags& other);

    // Clear all tags
    void Clear() noexcept;

    // Comparison
    bool operator==(const SingleCellTags& other) const noexcept;
    bool operator!=(const SingleCellTags& other) const noexcept;
};

// Configuration for barcode extraction
struct BarcodeExtractionConfig {
    // Expected barcode length (0 = auto-detect)
    std::uint32_t expected_barcode_length{16};

    // Expected UMI length (0 = auto-detect)
    std::uint32_t expected_umi_length{12};

    // Minimum quality for barcode bases (Phred scale)
    std::uint8_t min_barcode_quality{10};

    // Minimum quality for UMI bases
    std::uint8_t min_umi_quality{10};

    // Allow partial barcodes (shorter than expected)
    bool allow_partial_barcode{false};

    // Allow N bases in barcode
    bool allow_n_in_barcode{false};

    // Platform preset
    enum class Platform : std::uint8_t {
        Generic,     // No assumptions
        TenX,        // 10x Genomics
        Parse,       // Parse Biosciences
        Kinnex,      // PacBio Kinnex/MAS-Seq
    };
    Platform platform{Platform::Generic};
};

// Statistics for barcode extraction
struct BarcodeExtractionStats {
    std::uint64_t total_reads{0};
    std::uint64_t reads_with_cb{0};
    std::uint64_t reads_with_umi{0};
    std::uint64_t reads_with_both{0};
    std::uint64_t low_quality_barcodes{0};
    std::uint64_t short_barcodes{0};
    std::uint64_t n_containing_barcodes{0};

    // Unique cell count
    std::uint64_t unique_cells{0};

    // Mean reads per cell
    float mean_reads_per_cell{0.0f};
};

// Extract single-cell tags from a SAM/BAM aux tag string
// Format: "CB:Z:ACGT\tUB:Z:TGCA\t..."
[[nodiscard]] SingleCellTags ExtractTagsFromAux(std::string_view aux_string);

// Extract single-cell tags from individual tag-value pairs
[[nodiscard]] SingleCellTags ExtractTagsFromPairs(
    const std::vector<std::pair<std::string, std::string>>& tag_pairs);

// Extract barcode from read name (10x-style: @readname_BARCODE_UMI)
[[nodiscard]] SingleCellTags ExtractTagsFromReadName(
    std::string_view read_name,
    const BarcodeExtractionConfig& config = {});

// Format single-cell tags as SAM aux string
[[nodiscard]] std::string FormatTagsAsAux(const SingleCellTags& tags);

// Format single-cell tags as individual SAM fields
[[nodiscard]] std::vector<std::pair<std::string, std::string>> FormatTagsAsPairs(
    const SingleCellTags& tags);

// Validate barcode against expected format
[[nodiscard]] bool ValidateBarcode(
    std::string_view barcode,
    const BarcodeExtractionConfig& config);

// Validate UMI against expected format
[[nodiscard]] bool ValidateUmi(
    std::string_view umi,
    const BarcodeExtractionConfig& config);

// Cell barcode whitelist for barcode correction
class CellBarcodeWhitelist {
public:
    // Load whitelist from file (one barcode per line)
    static std::optional<CellBarcodeWhitelist> LoadFromFile(
        const std::filesystem::path& path);

    // Create empty whitelist
    CellBarcodeWhitelist() = default;

    // Add a barcode to the whitelist
    void Add(std::string barcode);

    // Check if barcode is in whitelist
    [[nodiscard]] bool Contains(std::string_view barcode) const;

    // Find nearest barcode within edit distance (returns nullopt if none found)
    [[nodiscard]] std::optional<std::string> FindNearest(
        std::string_view barcode,
        std::uint32_t max_edit_distance = 1) const;

    // Number of barcodes in whitelist
    [[nodiscard]] std::size_t Size() const noexcept;

    // Check if empty
    [[nodiscard]] bool Empty() const noexcept;

private:
    std::unordered_set<std::string> barcodes_;
};

}  // namespace llmap::singlecell
