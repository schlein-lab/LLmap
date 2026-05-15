#include "annot/classifier_internal.h"

#include <cmath>

namespace llmap::annot::internal {

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

void DecodeMappingHints(const Json& obj, ParamOverride& p) {
    if (obj.type != Json::Type::Object) return;
    for (const auto& kv : obj.obj_v) {
        const std::string& k = kv.first;
        const Json& v = kv.second;
        if      (k == "k" && v.type == Json::Type::Number)        p.k = static_cast<uint8_t>(v.num_v);
        else if (k == "w" && v.type == Json::Type::Number)        p.w = static_cast<uint8_t>(v.num_v);
        else if (k == "max_occ" && v.type == Json::Type::Number)  p.max_occ = static_cast<uint32_t>(v.num_v);
        else if (k == "max_occ_multiplier" && v.type == Json::Type::Number)
            p.max_occ = static_cast<uint32_t>(5000 * v.num_v);
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

namespace {

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
    return std::nullopt;
}

}  // namespace

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

}  // namespace llmap::annot::internal
