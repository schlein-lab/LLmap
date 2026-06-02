// LLmap — Bulk (T2) JSONL loader for SegDup catalog.
//
// T2 records carry only the minimum LLmap needs to (a) know "this region
// is a segdup" and (b) route to coordinate lookup. Schema (per line):
//   {
//     "locus_id": "<id>",
//     "structural_architecture": "<enum>",
//     "coords": { "<assembly>": { "chrom":..., "start":..., "end":... } }
//   }
// Extra fields are tolerated. Empty lines and `#`-prefixed comments are
// skipped.

#include "catalog/segdup_catalog.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <sstream>

namespace llmap::catalog {

using json = nlohmann::json;

namespace {

template <typename T>
T get_or(const json& j, std::string_view key, T fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    try { return it->get<T>(); }
    catch (const std::exception&) { return fallback; }
}

GenomicCoords parse_coords(std::string_view assembly, const json& j) {
    GenomicCoords c;
    c.assembly = std::string(assembly);
    c.chrom = get_or<std::string>(j, "chrom", "");
    c.start = get_or<std::int64_t>(j, "start", 0);
    c.end   = get_or<std::int64_t>(j, "end", 0);
    std::string strand = get_or<std::string>(j, "strand", ".");
    c.strand = strand.empty() ? '.' : strand[0];
    return c;
}

}  // namespace

std::optional<SegDupCatalogEntry>
parse_bulk_jsonl_line(std::string_view json_text,
                      const std::filesystem::path& source_path,
                      std::string& err) {
    json j;
    try {
        j = json::parse(json_text);
    } catch (const std::exception& ex) {
        err = std::string("JSON parse error: ") + ex.what();
        return std::nullopt;
    }
    if (!j.is_object()) {
        err = "JSONL line is not an object";
        return std::nullopt;
    }
    SegDupCatalogEntry e;
    e.tier = SegDupCatalogEntry::Tier::T2_Bulk;
    e.source_path = source_path;
    e.locus_id = get_or<std::string>(j, "locus_id", "");
    e.structural_architecture =
        get_or<std::string>(j, "structural_architecture", "");
    if (e.locus_id.empty()) {
        err = "missing locus_id";
        return std::nullopt;
    }
    auto coords_it = j.find("coords");
    if (coords_it == j.end() || !coords_it->is_object()) {
        err = "missing coords";
        return std::nullopt;
    }
    for (auto it = coords_it->begin(); it != coords_it->end(); ++it) {
        if (!it->is_object()) continue;
        auto c = parse_coords(it.key(), *it);
        if (c.chrom.empty() || c.end <= c.start) continue;
        e.coords_by_assembly.emplace(it.key(), std::move(c));
    }
    if (e.coords_by_assembly.empty()) {
        err = "no valid coordinate entries";
        return std::nullopt;
    }
    return e;
}

LoadStatus SegDupCatalog::load_bulk_jsonl(const std::filesystem::path& path) {
    LoadStatus st;
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        st.error = "cannot open " + path.string();
        return st;
    }
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        // Skip blank / comment lines (lenient JSONL).
        std::string_view sv = line;
        std::size_t i = 0;
        while (i < sv.size() && (sv[i] == ' ' || sv[i] == '\t')) ++i;
        if (i == sv.size()) continue;
        if (sv[i] == '#') continue;

        std::string err;
        auto e = parse_bulk_jsonl_line(line, path, err);
        if (!e) {
            st.records_skipped += 1;
            if (st.error.empty()) {
                std::ostringstream msg;
                msg << "line " << line_no << ": " << err;
                st.error = msg.str();
            }
            continue;
        }
        add_entry(std::move(*e));
        st.records_loaded += 1;
    }
    st.ok = true;
    return st;
}

}  // namespace llmap::catalog
