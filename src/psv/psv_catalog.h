// LLmap — PSV (Paralog-Specific Variant) catalog.
//
// Manages a collection of PSV sites indexed for fast lookup by genomic
// position. Supports loading from BED/VCF-like formats and serialization.

#pragma once

#include "psv/psv_types.h"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::psv {

// PSV catalog: collection of PSV sites with position-based lookup
class PsvCatalog {
public:
    PsvCatalog() = default;

    // Add a PSV site to the catalog
    void AddSite(PsvSite site);

    // Look up PSV sites overlapping a genomic region
    [[nodiscard]] std::vector<const PsvSite*> GetSitesInRegion(
        std::string_view chrom,
        std::uint64_t start,
        std::uint64_t end) const;

    // Look up a single PSV site by ID
    [[nodiscard]] const PsvSite* GetSiteById(std::uint64_t psv_id) const;

    // Look up a single PSV site by position
    [[nodiscard]] const PsvSite* GetSiteAtPosition(
        std::string_view chrom,
        std::uint64_t position) const;

    // Get all paralogs known in this catalog
    [[nodiscard]] std::vector<std::string> GetParalogs() const;

    // Get all chromosomes with PSV sites
    [[nodiscard]] std::vector<std::string> GetChromosomes() const;

    // Get total number of PSV sites
    [[nodiscard]] std::size_t Size() const noexcept { return sites_.size(); }

    // Check if catalog is empty
    [[nodiscard]] bool Empty() const noexcept { return sites_.empty(); }

    // Clear all sites
    void Clear() noexcept;

    // Build spatial index for fast region queries
    void BuildIndex();

    // Check if index is built
    [[nodiscard]] bool IsIndexed() const noexcept { return indexed_; }

    // Get sites for iteration
    [[nodiscard]] const std::vector<PsvSite>& Sites() const noexcept {
        return sites_;
    }

private:
    std::vector<PsvSite> sites_;

    // ID -> index in sites_
    std::unordered_map<std::uint64_t, std::size_t> id_index_;

    // chrom -> (position -> index in sites_)
    std::unordered_map<std::string, std::unordered_map<std::uint64_t, std::size_t>> pos_index_;

    // chrom -> sorted positions for range queries
    std::unordered_map<std::string, std::vector<std::uint64_t>> sorted_positions_;

    bool indexed_{false};
};

// Load PSV catalog from a BED-like file
// Format: chrom\tpos\tref\tparalog1:allele1,paralog2:allele2,...\tinformativeness
[[nodiscard]] std::optional<PsvCatalog> LoadPsvCatalogFromBed(
    const std::filesystem::path& path);

// Load PSV catalog from VCF with paralog annotations
[[nodiscard]] std::optional<PsvCatalog> LoadPsvCatalogFromVcf(
    const std::filesystem::path& path);

// Save PSV catalog to BED-like file
bool SavePsvCatalogToBed(
    const PsvCatalog& catalog,
    const std::filesystem::path& path);

// Compute informativeness for a PSV site based on allele distribution
[[nodiscard]] float ComputeInformativeness(const PsvSite& site);

// Generate a PSV catalog from reference sequences and known paralog alignments
class PsvCatalogBuilder {
public:
    // Add a paralog sequence
    void AddParalogSequence(std::string paralog_id, std::string_view sequence);

    // Set the reference/canonical sequence
    void SetReferenceSequence(std::string_view sequence);

    // Build catalog by comparing sequences and finding discriminating positions
    [[nodiscard]] PsvCatalog Build(
        std::string_view chrom,
        std::uint64_t start_position = 0,
        float min_informativeness = 0.5f);

private:
    std::string reference_;
    std::unordered_map<std::string, std::string> paralogs_;
};

}  // namespace llmap::psv
