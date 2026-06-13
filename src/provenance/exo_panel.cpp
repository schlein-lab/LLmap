// LLmap — Exogenous contaminant reference panel implementation.

#include "provenance/exo_panel.h"

#include <fstream>
#include <sstream>

namespace llmap::provenance {

const char* ExoCategoryName(ExoCategory c) noexcept {
    switch (c) {
        case ExoCategory::Virus:    return "virus";
        case ExoCategory::Bacteria: return "bacteria";
        case ExoCategory::SpikeIn:  return "spikein";
        case ExoCategory::Fungus:   return "fungus";
        case ExoCategory::Unknown:  return "unknown";
    }
    return "unknown";
}

namespace {

ExoCategory ParseCategory(std::string_view s) {
    if (s == "virus")    return ExoCategory::Virus;
    if (s == "bacteria") return ExoCategory::Bacteria;
    if (s == "spikein")  return ExoCategory::SpikeIn;
    if (s == "fungus")   return ExoCategory::Fungus;
    return ExoCategory::Unknown;
}

ExoReference Make(std::string taxon, std::string ref_id, std::string pv_tag,
                  ExoCategory cat) {
    ExoReference r;
    r.taxon = std::move(taxon);
    r.ref_id = std::move(ref_id);
    r.pv_tag = std::move(pv_tag);
    r.category = cat;
    r.verified = false;  // built-in starter accession
    return r;
}

}  // namespace

void ExoPanel::LoadBuiltinStarter() {
    refs_ = {
        // EBV — the LCL artefact: most 1000G/HapMap samples are EBV-immortalised.
        Make("EBV", "NC_007605", "exo:ebv", ExoCategory::Virus),
        // PhiX — Illumina sequencing control, often not removed.
        Make("PhiX", "NC_001422", "spikein:phix", ExoCategory::SpikeIn),
        // Lambda — ONT DCS / control.
        Make("Lambda", "J02459", "spikein:lambda", ExoCategory::SpikeIn),
        // ERCC — RNA spike-in controls.
        Make("ERCC", "ERCC", "spikein:ercc", ExoCategory::SpikeIn),
        // Mycoplasma — systematic cell-line RNA-seq contaminant.
        Make("Mycoplasma", "NC_006360", "exo:mycoplasma", ExoCategory::Bacteria),
    };
}

bool ExoPanel::LoadTsv(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::vector<ExoReference> loaded;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        ExoReference r;
        std::string cat;
        if (!(ss >> r.taxon >> r.ref_id >> r.pv_tag >> cat)) continue;
        r.category = ParseCategory(cat);
        r.verified = true;  // curated panel is authoritative
        loaded.push_back(std::move(r));
    }
    if (loaded.empty()) return false;
    refs_ = std::move(loaded);
    return true;
}

const ExoReference* ExoPanel::Lookup(std::string_view ref_id) const {
    for (const auto& r : refs_) {
        if (r.ref_id == ref_id) return &r;
    }
    return nullptr;
}

}  // namespace llmap::provenance
