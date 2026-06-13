// LLmap — Exogenous contaminant reference panel — Block 3.
//
// Feeds the Layer-1 `exo:*` / `spikein:*` origin classes: a read whose best
// alignment is to a contaminant reference (EBV episome in LCL samples, PhiX /
// Lambda / ERCC spike-ins, Mycoplasma in cell-line RNA-seq, kitome taxa) is
// exogenous, not host. The panel (a) defines which references to add to the
// alignment index, and (b) maps a hit's ref_id → its provenance PV tag.
//
// Same loader pattern as pseudogene_catalog / numt_catalog: a small built-in
// STARTER of the common contaminants (verified=false accessions) + a loadable
// TSV for a production panel (e.g. the full kitome list from Salter et al.).
// Dependency-light (stdlib only).

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::provenance {

enum class ExoCategory : std::uint8_t {
    Unknown = 0,
    Virus,
    Bacteria,
    SpikeIn,
    Fungus,
};

[[nodiscard]] const char* ExoCategoryName(ExoCategory c) noexcept;

struct ExoReference {
    std::string taxon;     // human label, e.g. "EBV" / "PhiX" / "Mycoplasma"
    std::string ref_id;    // reference accession a read may map to, e.g. "NC_007605"
    std::string pv_tag;    // provenance tag, e.g. "exo:ebv" / "spikein:phix"
    ExoCategory category{ExoCategory::Unknown};
    bool        verified{false};  // false ⇒ built-in starter accession
};

class ExoPanel {
public:
    // Built-in common contaminants (EBV, PhiX, Lambda, Mycoplasma, ERCC).
    // verified=false; for production prefer LoadTsv() with a curated panel.
    void LoadBuiltinStarter();

    // TSV: taxon, ref_id, pv_tag, category ("virus"|"bacteria"|"spikein"|"fungus").
    [[nodiscard]] bool LoadTsv(const std::string& path);

    // The panel entry whose reference is `ref_id` (a read's best hit), or null.
    [[nodiscard]] const ExoReference* Lookup(std::string_view ref_id) const;

    [[nodiscard]] std::size_t Size() const noexcept { return refs_.size(); }
    [[nodiscard]] const std::vector<ExoReference>& Entries() const noexcept {
        return refs_;
    }

private:
    std::vector<ExoReference> refs_;
};

}  // namespace llmap::provenance
