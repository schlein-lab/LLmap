#include "annot/classifier.h"
#include "annot/classifier_internal.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <unordered_map>

namespace llmap::annot {

using internal::Json;
using internal::JsonReader;
using internal::SlurpFile;
using internal::MatchOne;
using internal::DecodeMappingHints;
using internal::DecodeSignature;

Classifier::Classifier(std::vector<ClassifierRule> rules)
    : rules_(std::move(rules)) {
    std::sort(rules_.begin(), rules_.end(),
        [](const ClassifierRule& a, const ClassifierRule& b) {
            return a.priority > b.priority;
        });
}

std::unique_ptr<Classifier> Classifier::Load(const std::filesystem::path& path) {
    auto contents = SlurpFile(path);
    if (!contents) {
        std::cerr << "[annot] failed to open " << path.string() << "\n";
        return nullptr;
    }
    JsonReader rdr(*contents);
    auto root = rdr.Parse();
    if (!root || root->type != Json::Type::Object) {
        std::cerr << "[annot] failed to parse JSON in " << path.string() << "\n";
        return nullptr;
    }

    const Json* rules_node = root->Get("rules");
    if (!rules_node || rules_node->type != Json::Type::Array) {
        std::cerr << "[annot] no 'rules' array in " << path.string() << "\n";
        return nullptr;
    }

    std::vector<ClassifierRule> rules;
    for (const auto& r : rules_node->arr_v) {
        if (r.type != Json::Type::Object) continue;
        ClassifierRule rule;
        if (auto* n = r.Get("region_name"))
            rule.region_name = n->AsString();
        if (auto* n = r.Get("priority"))
            rule.priority = static_cast<int>(n->AsNumber());
        if (auto* sig = r.Get("feature_signature"))
            rule.predicates = DecodeSignature(*sig);
        else if (auto* preds = r.Get("predicates"))
            rule.predicates = DecodeSignature(*preds);
        if (auto* hints = r.Get("mapping_hints"))
            DecodeMappingHints(*hints, rule.mapping_hints);
        rules.push_back(std::move(rule));
    }

    // Merge mapping_hints from sibling regions/*.json files when the rule
    // itself didn't carry hints inline.
    auto regions_dir = path.parent_path() / "regions";
    if (std::filesystem::exists(regions_dir) &&
        std::filesystem::is_directory(regions_dir)) {

        std::unordered_map<std::string, ParamOverride> hints_by_region;
        for (const auto& entry : std::filesystem::directory_iterator(regions_dir)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".json") continue;
            auto file_contents = SlurpFile(entry.path());
            if (!file_contents) continue;
            JsonReader rrdr(*file_contents);
            auto region_root = rrdr.Parse();
            if (!region_root || region_root->type != Json::Type::Object) continue;
            const Json* nm = region_root->Get("name");
            const Json* h = region_root->Get("mapping_hints");
            if (!nm || nm->type != Json::Type::String) continue;
            if (!h) continue;
            ParamOverride p;
            DecodeMappingHints(*h, p);
            hints_by_region[nm->AsString()] = p;
        }

        for (auto& rule : rules) {
            auto it = hints_by_region.find(rule.region_name);
            if (it == hints_by_region.end()) continue;
            ParamOverride merged = rule.mapping_hints;
            merged.OverlayOver(it->second);
            rule.mapping_hints = merged;
        }
    }

    return std::make_unique<Classifier>(std::move(rules));
}

const ClassifierRule* Classifier::Match(const WindowFeatures& feats) const {
    for (const auto& rule : rules_) {
        bool all_match = true;
        for (const auto& p : rule.predicates) {
            if (!MatchOne(feats, p)) { all_match = false; break; }
        }
        if (all_match) return &rule;
    }
    return nullptr;
}

std::vector<AnnotationInterval> Classifier::Classify(
    const std::vector<WindowFeatures>& features) const {

    std::vector<AnnotationInterval> out;
    if (features.empty()) return out;

    AnnotationInterval cur;
    bool have_cur = false;

    auto flush = [&]() {
        if (have_cur) out.push_back(std::move(cur));
        have_cur = false;
    };

    for (const auto& f : features) {
        const auto* r = Match(f);
        std::string name = r ? r->region_name : "unique_single_copy";
        ParamOverride hints = r ? r->mapping_hints : ParamOverride{};

        if (have_cur && cur.ref_id == f.ref_id && cur.end == f.start &&
            cur.region_name == name) {
            cur.end = f.end;
            continue;
        }
        flush();
        cur = AnnotationInterval{};
        cur.ref_id = f.ref_id;
        cur.start = f.start;
        cur.end = f.end;
        cur.region_name = std::move(name);
        cur.source = r ? "taxonomy" : "default";
        cur.layer = r ? AnnotationLayer::Taxonomy : AnnotationLayer::Default;
        cur.params = std::move(hints);
        have_cur = true;
    }
    flush();
    return out;
}

}  // namespace llmap::annot
