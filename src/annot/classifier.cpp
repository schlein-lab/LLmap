#include "annot/classifier.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

namespace llmap::annot {

namespace {

// --- tiny JSON reader -----------------------------------------------------
//
// Supports: objects, arrays, strings (no escape handling beyond \"), numbers
// (decimal), bools, null. Whitespace tolerated. Comments not supported (we
// stick to standard JSON for the rule files).

struct Json {
    enum class Type { Null, Bool, Number, String, Array, Object };
    Type type = Type::Null;
    bool bool_v = false;
    double num_v = 0.0;
    std::string str_v;
    std::vector<Json> arr_v;
    std::vector<std::pair<std::string, Json>> obj_v;

    const Json* Get(const std::string& key) const {
        if (type != Type::Object) return nullptr;
        for (const auto& kv : obj_v) if (kv.first == key) return &kv.second;
        return nullptr;
    }
    bool AsBool(bool defv = false) const {
        if (type == Type::Bool) return bool_v;
        return defv;
    }
    double AsNumber(double defv = 0.0) const {
        if (type == Type::Number) return num_v;
        return defv;
    }
    std::string AsString(const std::string& defv = "") const {
        if (type == Type::String) return str_v;
        return defv;
    }
};

class JsonReader {
public:
    explicit JsonReader(std::string src) : src_(std::move(src)), pos_(0) {}
    std::optional<Json> Parse() {
        SkipWs();
        Json out;
        if (!ParseValue(out)) return std::nullopt;
        SkipWs();
        return out;
    }

private:
    std::string src_;
    size_t pos_;

    void SkipWs() {
        while (pos_ < src_.size() &&
               (std::isspace(static_cast<unsigned char>(src_[pos_])))) ++pos_;
    }
    bool Peek(char c) {
        SkipWs();
        return pos_ < src_.size() && src_[pos_] == c;
    }
    bool Eat(char c) {
        SkipWs();
        if (pos_ < src_.size() && src_[pos_] == c) { ++pos_; return true; }
        return false;
    }
    bool ParseValue(Json& out) {
        SkipWs();
        if (pos_ >= src_.size()) return false;
        char c = src_[pos_];
        if (c == '{') return ParseObject(out);
        if (c == '[') return ParseArray(out);
        if (c == '"') return ParseString(out);
        if (c == 't' || c == 'f') return ParseBool(out);
        if (c == 'n') return ParseNull(out);
        if (c == '-' || std::isdigit(static_cast<unsigned char>(c)))
            return ParseNumber(out);
        return false;
    }
    bool ParseString(Json& out) {
        if (!Eat('"')) return false;
        std::string s;
        while (pos_ < src_.size() && src_[pos_] != '"') {
            if (src_[pos_] == '\\' && pos_ + 1 < src_.size()) {
                char n = src_[pos_ + 1];
                switch (n) {
                    case '"': s.push_back('"'); break;
                    case '\\': s.push_back('\\'); break;
                    case 'n': s.push_back('\n'); break;
                    case 't': s.push_back('\t'); break;
                    case 'r': s.push_back('\r'); break;
                    default: s.push_back(n);
                }
                pos_ += 2;
            } else {
                s.push_back(src_[pos_]);
                ++pos_;
            }
        }
        if (!Eat('"')) return false;
        out.type = Json::Type::String;
        out.str_v = std::move(s);
        return true;
    }
    bool ParseNumber(Json& out) {
        size_t start = pos_;
        if (src_[pos_] == '-') ++pos_;
        while (pos_ < src_.size() &&
               (std::isdigit(static_cast<unsigned char>(src_[pos_])) ||
                src_[pos_] == '.' || src_[pos_] == 'e' || src_[pos_] == 'E' ||
                src_[pos_] == '+' || src_[pos_] == '-')) ++pos_;
        try {
            out.num_v = std::stod(src_.substr(start, pos_ - start));
        } catch (...) { return false; }
        out.type = Json::Type::Number;
        return true;
    }
    bool ParseBool(Json& out) {
        if (src_.compare(pos_, 4, "true") == 0) {
            pos_ += 4;
            out.type = Json::Type::Bool;
            out.bool_v = true;
            return true;
        }
        if (src_.compare(pos_, 5, "false") == 0) {
            pos_ += 5;
            out.type = Json::Type::Bool;
            out.bool_v = false;
            return true;
        }
        return false;
    }
    bool ParseNull(Json& out) {
        if (src_.compare(pos_, 4, "null") == 0) {
            pos_ += 4;
            out.type = Json::Type::Null;
            return true;
        }
        return false;
    }
    bool ParseArray(Json& out) {
        if (!Eat('[')) return false;
        out.type = Json::Type::Array;
        SkipWs();
        if (Eat(']')) return true;
        while (true) {
            Json elem;
            if (!ParseValue(elem)) return false;
            out.arr_v.push_back(std::move(elem));
            SkipWs();
            if (Eat(',')) continue;
            if (Eat(']')) return true;
            return false;
        }
    }
    bool ParseObject(Json& out) {
        if (!Eat('{')) return false;
        out.type = Json::Type::Object;
        SkipWs();
        if (Eat('}')) return true;
        while (true) {
            Json key;
            SkipWs();
            if (!ParseString(key)) return false;
            SkipWs();
            if (!Eat(':')) return false;
            Json val;
            if (!ParseValue(val)) return false;
            out.obj_v.emplace_back(std::move(key.str_v), std::move(val));
            SkipWs();
            if (Eat(',')) continue;
            if (Eat('}')) return true;
            return false;
        }
    }
};

// --- predicate evaluation -------------------------------------------------

bool MatchOne(const WindowFeatures& f, const FeaturePredicate& p) {
    auto getf = [&](const std::string& name, float& out) -> bool {
        if      (name == "shannon_5mer") { out = f.shannon_5mer; return true; }
        else if (name == "gc_content")   { out = f.gc_content; return true; }
        else if (name == "palindrome_density") { out = f.palindrome_density; return true; }
        else if (name == "orf_density")  { out = f.orf_density; return true; }
        else if (name == "tandem_period") {
            if (f.tandem_period < 0) return false;
            out = static_cast<float>(f.tandem_period); return true;
        }
        else if (name == "kmer_multiplicity_p95") {
            out = static_cast<float>(f.kmer_multiplicity_p95); return true;
        }
        return false;
    };

    float fv;
    switch (p.op) {
        case FeaturePredicate::Op::GE:
            if (!getf(p.feature, fv)) return false;
            return fv >= p.bound_lo;
        case FeaturePredicate::Op::LE:
            if (!getf(p.feature, fv)) return false;
            return fv <= p.bound_hi;
        case FeaturePredicate::Op::EQ:
            if (!getf(p.feature, fv)) return false;
            return std::abs(fv - p.bound_lo) < 1e-6f;
        case FeaturePredicate::Op::RangeIn:
            if (!getf(p.feature, fv)) return false;
            return fv >= p.bound_lo && fv <= p.bound_hi;
        case FeaturePredicate::Op::MultiplicityMin:
            return f.kmer_multiplicity_p95 >= static_cast<uint32_t>(p.int_bound);
    }
    return false;
}

// Apply mapping_hints JSON object to ParamOverride.
void DecodeMappingHints(const Json& obj, ParamOverride& p) {
    if (obj.type != Json::Type::Object) return;
    for (const auto& kv : obj.obj_v) {
        const std::string& k = kv.first;
        const Json& v = kv.second;
        if      (k == "k" && v.type == Json::Type::Number)        p.k = static_cast<uint8_t>(v.num_v);
        else if (k == "w" && v.type == Json::Type::Number)        p.w = static_cast<uint8_t>(v.num_v);
        else if (k == "max_occ" && v.type == Json::Type::Number)  p.max_occ = static_cast<uint32_t>(v.num_v);
        else if (k == "max_occ_multiplier" && v.type == Json::Type::Number)
            p.max_occ = static_cast<uint32_t>(5000 * v.num_v);  // scale against typical 5000
        else if (k == "lambda_scale" && v.type == Json::Type::Number)
            p.lambda_scale = static_cast<float>(v.num_v);
        else if (k == "identity_threshold" && v.type == Json::Type::Number)
            p.identity_threshold = static_cast<float>(v.num_v);
        else if (k == "anchor_weight_scale" && v.type == Json::Type::Number)
            p.anchor_weight_scale = static_cast<float>(v.num_v);
        else if (k == "report_multi_position")  p.report_multi_position = v.AsBool();
        else if (k == "require_psv_disambiguation") p.require_psv_disambig = v.AsBool();
        else if (k == "allow_high_mismatch") p.allow_high_mismatch = v.AsBool();
        else if (k == "require_llm_at_runtime") p.require_llm_at_runtime = v.AsBool();
    }
}

// Parse a single feature_signature predicate of either of two syntaxes:
//   (a) human-style:  {"feature": "x", "op": "between"|"eq"|"lt"|"gt"|"le"|"ge"|"in", "value": ...}
//   (b) compact:      {"x": {"range": [lo, hi]}} or {"x": {"min": v}} / {"max": v}
//
// Returns std::nullopt for unrecognised forms; caller silently drops them.
std::optional<FeaturePredicate> DecodeOnePredicateA(const Json& node) {
    auto* feat_n = node.Get("feature");
    auto* op_n = node.Get("op");
    auto* val_n = node.Get("value");
    if (!feat_n || !op_n || feat_n->type != Json::Type::String ||
        op_n->type != Json::Type::String) {
        return std::nullopt;
    }
    FeaturePredicate p;
    p.feature = feat_n->str_v;
    const std::string& op = op_n->str_v;

    if (op == "between" && val_n && val_n->type == Json::Type::Array &&
        val_n->arr_v.size() == 2) {
        p.op = FeaturePredicate::Op::RangeIn;
        p.bound_lo = static_cast<float>(val_n->arr_v[0].AsNumber());
        p.bound_hi = static_cast<float>(val_n->arr_v[1].AsNumber());
        return p;
    }
    if ((op == "ge" || op == "gte" || op == "geq" || op == "gt") && val_n &&
        val_n->type == Json::Type::Number) {
        p.op = FeaturePredicate::Op::GE;
        p.bound_lo = static_cast<float>(val_n->num_v);
        return p;
    }
    if ((op == "le" || op == "lte" || op == "leq" || op == "lt") && val_n &&
        val_n->type == Json::Type::Number) {
        p.op = FeaturePredicate::Op::LE;
        p.bound_hi = static_cast<float>(val_n->num_v);
        return p;
    }
    if (op == "eq" && val_n && val_n->type == Json::Type::Number) {
        p.op = FeaturePredicate::Op::EQ;
        p.bound_lo = static_cast<float>(val_n->num_v);
        return p;
    }
    if (op == "min" && val_n && val_n->type == Json::Type::Number) {
        if (p.feature == "kmer_multiplicity_p95") {
            p.op = FeaturePredicate::Op::MultiplicityMin;
            p.int_bound = static_cast<int>(val_n->num_v);
        } else {
            p.op = FeaturePredicate::Op::GE;
            p.bound_lo = static_cast<float>(val_n->num_v);
        }
        return p;
    }
    // Note: "in" with a string-set value is dropped here because we don't yet
    // have consensus_match feature plumbing; the corresponding rule simply
    // matches via the other predicates of the conjunction.
    return std::nullopt;
}

// Parse a feature_signature node in either array-of-predicates form
// (style a) or compact object form (style b).
std::vector<FeaturePredicate> DecodeSignature(const Json& sig) {
    std::vector<FeaturePredicate> out;
    if (sig.type == Json::Type::Array) {
        for (const auto& p : sig.arr_v) {
            auto pred = DecodeOnePredicateA(p);
            if (pred) out.push_back(std::move(*pred));
        }
        return out;
    }
    if (sig.type != Json::Type::Object) return out;
    for (const auto& kv : sig.obj_v) {
        const std::string& feat = kv.first;
        const Json& v = kv.second;
        if (v.type != Json::Type::Object) continue;

        const Json* range_node = v.Get("range");
        const Json* min_node = v.Get("min");
        const Json* max_node = v.Get("max");

        if (range_node && range_node->type == Json::Type::Array &&
            range_node->arr_v.size() == 2 &&
            range_node->arr_v[0].type == Json::Type::Number &&
            range_node->arr_v[1].type == Json::Type::Number) {
            FeaturePredicate p;
            p.feature = feat;
            p.op = FeaturePredicate::Op::RangeIn;
            p.bound_lo = static_cast<float>(range_node->arr_v[0].num_v);
            p.bound_hi = static_cast<float>(range_node->arr_v[1].num_v);
            out.push_back(std::move(p));
        }
        if (min_node && min_node->type == Json::Type::Number) {
            FeaturePredicate p;
            p.feature = feat;
            if (feat == "kmer_multiplicity_p95") {
                p.op = FeaturePredicate::Op::MultiplicityMin;
                p.int_bound = static_cast<int>(min_node->num_v);
            } else {
                p.op = FeaturePredicate::Op::GE;
                p.bound_lo = static_cast<float>(min_node->num_v);
            }
            out.push_back(std::move(p));
        }
        if (max_node && max_node->type == Json::Type::Number) {
            FeaturePredicate p;
            p.feature = feat;
            p.op = FeaturePredicate::Op::LE;
            p.bound_hi = static_cast<float>(max_node->num_v);
            out.push_back(std::move(p));
        }
    }
    return out;
}

}  // namespace

Classifier::Classifier(std::vector<ClassifierRule> rules)
    : rules_(std::move(rules)) {
    std::sort(rules_.begin(), rules_.end(),
        [](const ClassifierRule& a, const ClassifierRule& b) {
            return a.priority > b.priority;
        });
}

std::unique_ptr<Classifier> Classifier::Load(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "[annot] failed to open " << path.string() << "\n";
        return nullptr;
    }
    std::stringstream ss;
    ss << f.rdbuf();
    JsonReader rdr(ss.str());
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
