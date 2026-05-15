#include "annot/annotation_store.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

namespace llmap::annot {

namespace {

void EncodeParamPair(std::ostringstream& oss, const char* key, bool& first,
                     auto value) {
    if (!first) oss << ',';
    oss << key << '=' << value;
    first = false;
}

std::string EncodeParams(const ParamOverride& p) {
    std::ostringstream oss;
    bool first = true;
    if (p.k)                       EncodeParamPair(oss, "k",                   first, static_cast<int>(*p.k));
    if (p.w)                       EncodeParamPair(oss, "w",                   first, static_cast<int>(*p.w));
    if (p.max_occ)                 EncodeParamPair(oss, "max_occ",             first, *p.max_occ);
    if (p.lambda_scale)            EncodeParamPair(oss, "lambda_scale",        first, *p.lambda_scale);
    if (p.identity_threshold)      EncodeParamPair(oss, "identity_threshold",  first, *p.identity_threshold);
    if (p.anchor_weight_scale)     EncodeParamPair(oss, "anchor_weight_scale", first, *p.anchor_weight_scale);
    if (p.report_multi_position)   EncodeParamPair(oss, "report_multi_position", first, *p.report_multi_position ? 1 : 0);
    if (p.require_psv_disambig)    EncodeParamPair(oss, "require_psv_disambig", first, *p.require_psv_disambig ? 1 : 0);
    if (p.allow_high_mismatch)     EncodeParamPair(oss, "allow_high_mismatch",  first, *p.allow_high_mismatch ? 1 : 0);
    if (p.require_llm_at_runtime)  EncodeParamPair(oss, "require_llm_at_runtime", first, *p.require_llm_at_runtime ? 1 : 0);
    return oss.str();
}

void ApplyParam(ParamOverride& p, const std::string& key,
                const std::string& value) {
    auto as_int = [&]() { return std::stoi(value); };
    auto as_uint = [&]() {
        return static_cast<uint32_t>(std::stoul(value));
    };
    auto as_float = [&]() { return std::stof(value); };
    auto as_bool = [&]() { return value == "1" || value == "true"; };
    if      (key == "k")                      p.k = static_cast<uint8_t>(as_int());
    else if (key == "w")                      p.w = static_cast<uint8_t>(as_int());
    else if (key == "max_occ")                p.max_occ = as_uint();
    else if (key == "lambda_scale")           p.lambda_scale = as_float();
    else if (key == "identity_threshold")     p.identity_threshold = as_float();
    else if (key == "anchor_weight_scale")    p.anchor_weight_scale = as_float();
    else if (key == "report_multi_position")  p.report_multi_position = as_bool();
    else if (key == "require_psv_disambig")   p.require_psv_disambig = as_bool();
    else if (key == "allow_high_mismatch")    p.allow_high_mismatch = as_bool();
    else if (key == "require_llm_at_runtime") p.require_llm_at_runtime = as_bool();
    // Unknown key: silently ignore for forward compatibility.
}

ParamOverride DecodeParams(const std::string& s) {
    ParamOverride p;
    if (s.empty() || s == "-") return p;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        auto eq = item.find('=');
        if (eq == std::string::npos) continue;
        ApplyParam(p, item.substr(0, eq), item.substr(eq + 1));
    }
    return p;
}

AnnotationLayer DecodeLayer(int v) {
    switch (v) {
        case 1: return AnnotationLayer::Taxonomy;
        case 2: return AnnotationLayer::SpecificLocus;
        case 3: return AnnotationLayer::AgentDecision;
        case 4: return AnnotationLayer::Stochastic;
        default: return AnnotationLayer::Default;
    }
}

}  // namespace

std::unique_ptr<AnnotationStore> AnnotationStore::Create(
    std::vector<AnnotationInterval> intervals,
    std::vector<std::string> contig_names) {
    auto store = std::make_unique<AnnotationStore>();
    store->contig_names_ = std::move(contig_names);
    store->index_.Build(std::move(intervals));
    return store;
}

std::unique_ptr<AnnotationStore> AnnotationStore::Load(
    const std::filesystem::path& path,
    const std::vector<std::string>& contig_names) {

    std::ifstream f(path);
    if (!f) return nullptr;

    std::unordered_map<std::string, uint32_t> name_to_id;
    for (size_t i = 0; i < contig_names.size(); ++i) {
        name_to_id[contig_names[i]] = static_cast<uint32_t>(i);
    }

    auto store = std::make_unique<AnnotationStore>();
    store->contig_names_ = contig_names;

    std::vector<AnnotationInterval> intervals;
    std::string line;
    size_t line_no = 0;
    size_t skipped = 0;
    while (std::getline(f, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;

        // tab-split: ref_name, start, end, region_name, source, layer, params
        std::vector<std::string> cols;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, '\t')) cols.push_back(std::move(tok));
        if (cols.size() < 6) {
            ++skipped;
            continue;
        }

        auto it = name_to_id.find(cols[0]);
        if (it == name_to_id.end()) {
            ++skipped;
            continue;
        }

        AnnotationInterval iv;
        iv.ref_id = it->second;
        try {
            iv.start = static_cast<uint32_t>(std::stoul(cols[1]));
            iv.end = static_cast<uint32_t>(std::stoul(cols[2]));
        } catch (...) {
            ++skipped;
            continue;
        }
        iv.region_name = cols[3];
        iv.source = cols[4];
        try {
            iv.layer = DecodeLayer(std::stoi(cols[5]));
        } catch (...) {
            iv.layer = AnnotationLayer::Default;
        }
        if (cols.size() > 6) iv.params = DecodeParams(cols[6]);
        intervals.push_back(std::move(iv));
    }

    if (skipped > 0) {
        std::cerr << "[annot] skipped " << skipped
                  << " malformed lines in " << path.string() << "\n";
    }

    store->index_.Build(std::move(intervals));
    return store;
}

bool AnnotationStore::Save(const std::filesystem::path& path) const {
    std::ofstream f(path);
    if (!f) return false;
    f << "# LLmap region annotation -- BED-like, layer 0=default 1=taxonomy "
         "2=specific_locus 3=agent 4=stochastic\n";
    f << "# ref_name\tstart\tend\tregion_name\tsource\tlayer\tparams\n";
    for (const auto& iv : index_.Intervals()) {
        const std::string& cn = iv.ref_id < contig_names_.size()
            ? contig_names_[iv.ref_id]
            : ("ref" + std::to_string(iv.ref_id));
        std::string params = EncodeParams(iv.params);
        f << cn << '\t' << iv.start << '\t' << iv.end << '\t'
          << iv.region_name << '\t' << iv.source << '\t'
          << static_cast<int>(iv.layer) << '\t'
          << (params.empty() ? "-" : params) << '\n';
    }
    return true;
}

}  // namespace llmap::annot
