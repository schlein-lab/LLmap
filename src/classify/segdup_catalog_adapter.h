// LLmap — Adapter from llmap::catalog::SegDupCatalog (real loader) to
// llmap::classify::igh::ISegDupCatalog (read-only classifier interface).
//
// Bridges the two parallel agents:
//   * src/catalog/   owns the JSON loader + storage
//   * src/classify/  owns the IGH-mode upstream-anchor classifier
//
// The adapter holds a non-owning pointer/reference to the catalog and
// translates lookups into the `CatalogEntryView` shape the classifier reads.
// Adapter lifetime must not exceed the wrapped catalog's lifetime.
//
// Naming: `_adapter` (not `_bridge`) — Gang-of-Four adapter pattern.

#pragma once

#include "classify/segdup_catalog_iface.h"

#include <optional>
#include <string_view>

namespace llmap {
namespace catalog { class SegDupCatalog; }       // fwd-decl

namespace classify::igh {

/// @brief Adapts a `llmap::catalog::SegDupCatalog` to `ISegDupCatalog`.
///
/// Constructor takes a const reference; assembly defaults to "GRCh38" because
/// every T1 entry must declare GRCh38 coordinates per catalog schema v0.2.
/// Pass another assembly to query CHM13 / pangenome assemblies.
class SegDupCatalogAdapter final : public ISegDupCatalog {
public:
    explicit SegDupCatalogAdapter(const catalog::SegDupCatalog& cat,
                                   std::string_view assembly = "GRCh38")
        : cat_(cat), assembly_(assembly) {}

    [[nodiscard]] std::optional<CatalogEntryView> LookupByPosition(
        std::string_view target_id, std::uint64_t pos) const override;

    [[nodiscard]] std::optional<CatalogEntryView>
    LookupByHaplotypeClass(std::string_view haplotype_class) const override;

private:
    const catalog::SegDupCatalog& cat_;
    std::string assembly_;
};

}  // namespace classify::igh
}  // namespace llmap
