// LLmap — Curated (T1) JSON parser for SegDup catalog.
//
// Uses nlohmann/json. Tolerant of missing optional fields: only locus_id
// + structural_architecture + coords are required. Schema version 0.1
// and 0.2 both parse through the same path because v0.2 only adds new
// fields without renaming existing ones.

#include "catalog/segdup_catalog.h"

#include <nlohmann/json.hpp>

#include <sstream>

namespace llmap::catalog {

using json = nlohmann::json;

namespace {

// --- small helpers --------------------------------------------------------

template <typename T>
T get_or(const json& j, std::string_view key, T fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    try { return it->get<T>(); }
    catch (const std::exception&) { return fallback; }
}

std::vector<std::string> get_string_array(const json& j, std::string_view key) {
    std::vector<std::string> out;
    auto it = j.find(key);
    if (it == j.end() || !it->is_array()) return out;
    for (const auto& v : *it) {
        if (v.is_string()) out.push_back(v.get<std::string>());
    }
    return out;
}

GenomicCoords parse_coords(std::string_view assembly, const json& j) {
    GenomicCoords c;
    c.assembly = std::string(assembly);
    c.chrom = get_or<std::string>(j, "chrom", "");
    c.start = get_or<std::int64_t>(j, "start", 0);
    c.end   = get_or<std::int64_t>(j, "end", 0);
    std::string strand = get_or<std::string>(j, "strand", ".");
    c.strand = strand.empty() ? '.' : strand[0];
    c.note = get_or<std::string>(j, "note", "");
    return c;
}

DiagnosticSnp parse_snp(const json& j) {
    DiagnosticSnp s;
    s.id          = get_or<std::string>(j, "id", "");
    s.ref         = get_or<std::string>(j, "ref", "");
    s.alt         = get_or<std::string>(j, "alt", "");
    s.consequence = get_or<std::string>(j, "consequence", "");
    s.snp_class   = get_or<std::string>(j, "class", "");
    s.freq_canonical = get_or<double>(j, "frequency_in_canonical_haps", -1.0);
    s.freq_dup       = get_or<double>(j, "frequency_in_dup_haps", -1.0);

    // position is a sub-object; tolerant of either nested or flat layout.
    auto pos_it = j.find("position");
    if (pos_it != j.end() && pos_it->is_object()) {
        s.transcript = get_or<std::string>(*pos_it, "transcript", "");
        s.cds_offset = get_or<std::int64_t>(*pos_it, "cds_offset", -1);
        // Accept either "grch38_chr14"-style flat key or generic
        // "grch38_pos" / "chrom" + "position".
        for (auto it = pos_it->begin(); it != pos_it->end(); ++it) {
            const std::string& key = it.key();
            if (key.rfind("grch38_chr", 0) == 0 && it->is_number_integer()) {
                s.grch38_pos = it->get<std::int64_t>();
                s.chrom = key.substr(7);  // "chr14"
            }
        }
        if (s.grch38_pos < 0) {
            s.grch38_pos = get_or<std::int64_t>(*pos_it, "grch38_pos", -1);
            s.chrom      = get_or<std::string>(*pos_it, "chrom", s.chrom);
        }
    }
    return s;
}

MappingPrimary parse_primary(const json& j) {
    MappingPrimary p;
    p.kmer_size           = get_or<std::int32_t>(j, "kmer_size", 0);
    p.max_mismatch        = get_or<std::int32_t>(j, "max_mismatch", 0);
    p.include_flanking_bp = get_or<std::int32_t>(j, "include_flanking_bp", 0);
    p.include_flanking_anchor =
        get_or<std::string>(j, "include_flanking_anchor", "");
    p.require_unique_chain = get_or<bool>(j, "require_unique_chain", false);
    p.note = get_or<std::string>(j, "note", "");
    return p;
}

MappingFallbackStage parse_stage(const json& j) {
    MappingFallbackStage st;
    st.stage     = get_or<std::int32_t>(j, "stage", 0);
    st.name      = get_or<std::string>(j, "name", "");
    st.trigger   = get_or<std::string>(j, "trigger", "");
    st.rationale = get_or<std::string>(j, "rationale", "");
    if (j.contains("kmer_size") && j["kmer_size"].is_number_integer()) {
        st.kmer_size = j["kmer_size"].get<std::int32_t>();
    }
    if (j.contains("max_mismatch") && j["max_mismatch"].is_number_integer()) {
        st.max_mismatch = j["max_mismatch"].get<std::int32_t>();
    }
    if (j.contains("k") && j["k"].is_number_integer()) {
        st.top_k = j["k"].get<std::int32_t>();
    }
    if (j.contains("use_extension") && j["use_extension"].is_boolean()) {
        st.use_extension = j["use_extension"].get<bool>();
    }
    if (j.contains("emit_warning") && j["emit_warning"].is_boolean()) {
        st.emit_warning = j["emit_warning"].get<bool>();
    }
    if (j.contains("opt_in_flag") && j["opt_in_flag"].is_string()) {
        st.opt_in_flag = j["opt_in_flag"].get<std::string>();
    }
    return st;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public parse_curated_json
// ---------------------------------------------------------------------------

std::optional<SegDupCatalogEntry>
parse_curated_json(std::string_view json_text,
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
        err = "top-level JSON is not an object";
        return std::nullopt;
    }

    SegDupCatalogEntry e;
    e.tier = SegDupCatalogEntry::Tier::T1_Curated;
    e.source_path = source_path;
    e.locus_id       = get_or<std::string>(j, "locus_id", "");
    e.human_name     = get_or<std::string>(j, "human_name", "");
    e.version        = get_or<std::string>(j, "version", "");
    e.schema_version = get_or<std::string>(j, "schema_version", "");
    e.structural_architecture =
        get_or<std::string>(j, "structural_architecture", "");
    e.haplotype_class = get_or<std::string>(j, "haplotype_class", "");
    e.mechanism         = get_string_array(j, "mechanism");
    e.clinical_function = get_string_array(j, "clinical_function");
    e.nahr_status = get_or<std::string>(j, "nahr_status", "");

    if (e.locus_id.empty()) {
        err = "missing required field locus_id";
        return std::nullopt;
    }
    if (e.structural_architecture.empty()) {
        err = "missing required field structural_architecture";
        return std::nullopt;
    }

    // coords
    auto coords_it = j.find("coords");
    if (coords_it == j.end() || !coords_it->is_object()) {
        err = "missing or non-object coords";
        return std::nullopt;
    }
    for (auto it = coords_it->begin(); it != coords_it->end(); ++it) {
        if (!it->is_object()) continue;  // null assembly entries permitted
        auto coords = parse_coords(it.key(), *it);
        if (coords.chrom.empty()) continue;
        e.coords_by_assembly.emplace(it.key(), std::move(coords));
    }
    if (e.coords_by_assembly.empty()) {
        err = "coords object has no non-null assemblies";
        return std::nullopt;
    }

    // diagnostic_features.discriminating_snps (optional)
    auto df_it = j.find("diagnostic_features");
    if (df_it != j.end() && df_it->is_object()) {
        auto snps_it = df_it->find("discriminating_snps");
        if (snps_it != df_it->end() && snps_it->is_array()) {
            for (const auto& s : *snps_it) {
                if (s.is_object()) e.discriminating_snps.push_back(parse_snp(s));
            }
        }

        // diagnostic_features.promoter_signature (optional, schema v0.2)
        auto ps_it = df_it->find("promoter_signature");
        if (ps_it != df_it->end() && ps_it->is_object()) {
            PromoterSignature ps;
            ps.window_relative_to_anchor =
                get_or<std::string>(*ps_it, "window_relative_to_anchor", "");
            ps.anchor = get_or<std::string>(*ps_it, "anchor", "");
            auto cm = ps_it->find("canonical_motif");
            if (cm != ps_it->end() && cm->is_string()) {
                ps.canonical_motif = cm->get<std::string>();
            }
            auto dm = ps_it->find("duplicate_motif");
            if (dm != ps_it->end() && dm->is_string()) {
                ps.duplicate_motif = dm->get<std::string>();
            }
            ps.note = get_or<std::string>(*ps_it, "note", "");
            e.promoter_signature = std::move(ps);
        }
    }

    // mapping_strategy (optional)
    auto ms_it = j.find("mapping_strategy");
    if (ms_it != j.end() && ms_it->is_object()) {
        auto pr_it = ms_it->find("primary");
        if (pr_it != ms_it->end() && pr_it->is_object()) {
            e.mapping_primary = parse_primary(*pr_it);
        }
        auto fc_it = ms_it->find("fallback_chain");
        if (fc_it != ms_it->end() && fc_it->is_array()) {
            for (const auto& stage : *fc_it) {
                if (stage.is_object()) {
                    e.fallback_chain.push_back(parse_stage(stage));
                }
            }
        }
    }
    return e;
}

}  // namespace llmap::catalog
