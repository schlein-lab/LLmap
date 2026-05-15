// LLmap -- rule-based classifier: features -> region label + ParamOverride.
//
// Rules are loaded from classifier_rules.json (one per organism). Each rule
// is a conjunction of feature predicates; if all predicates match, the rule
// fires with a given region_name + mapping_hints (which become the
// ParamOverride for that window).
//
// Higher-priority rules win on conflict.
//
// The rules file format is hand-readable JSON, written via a minimal
// recursive-descent parser (no external dependency).

#pragma once

#include "annot/annot_types.h"
#include "annot/feature_extractor.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace llmap::annot {

struct FeaturePredicate {
    enum class Op { GE, LE, EQ, RangeIn, MultiplicityMin };
    std::string feature;
    Op op;
    float bound_lo = 0.0f;
    float bound_hi = 0.0f;
    int int_bound = 0;
    std::string str_bound;
};

struct ClassifierRule {
    std::string region_name;
    int priority = 0;
    std::vector<FeaturePredicate> predicates;
    ParamOverride mapping_hints;
};

class Classifier {
public:
    // Load rules from a classifier_rules.json file. Returns nullptr on parse
    // failure; failure message goes to stderr.
    static std::unique_ptr<Classifier> Load(const std::filesystem::path& path);

    // Construct with explicit rules (for unit tests).
    explicit Classifier(std::vector<ClassifierRule> rules);

    // Apply rules to a window's features. Returns the highest-priority rule
    // that matches; nullptr if none matched.
    const ClassifierRule* Match(const WindowFeatures& feats) const;

    // Convenience: classify every window and build an AnnotationInterval set
    // ready for AnnotationStore. Adjacent windows with the same region_name
    // are merged.
    std::vector<AnnotationInterval> Classify(
        const std::vector<WindowFeatures>& features) const;

    size_t NumRules() const { return rules_.size(); }

private:
    std::vector<ClassifierRule> rules_;
};

}  // namespace llmap::annot
