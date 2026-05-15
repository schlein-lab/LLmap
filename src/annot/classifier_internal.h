// Internal header for classifier implementation - shared types and utilities.
#pragma once

#include "annot/classifier.h"
#include "annot/annot_types.h"
#include "annot/feature_extractor.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace llmap::annot {
namespace internal {

// --- tiny JSON reader (supports: objects, arrays, strings, numbers, bools, null)
struct Json {
    enum class Type { Null, Bool, Number, String, Array, Object };
    Type type = Type::Null;
    bool bool_v = false;
    double num_v = 0.0;
    std::string str_v;
    std::vector<Json> arr_v;
    std::vector<std::pair<std::string, Json>> obj_v;

    const Json* Get(const std::string& key) const;
    bool AsBool(bool defv = false) const;
    double AsNumber(double defv = 0.0) const;
    std::string AsString(const std::string& defv = "") const;
};

class JsonReader {
public:
    explicit JsonReader(std::string src);
    std::optional<Json> Parse();

private:
    std::string src_;
    size_t pos_;

    void SkipWs();
    bool Peek(char c);
    bool Eat(char c);
    bool ParseValue(Json& out);
    bool ParseString(Json& out);
    bool ParseNumber(Json& out);
    bool ParseBool(Json& out);
    bool ParseNull(Json& out);
    bool ParseArray(Json& out);
    bool ParseObject(Json& out);
};

// Read whole file into string, or return std::nullopt.
std::optional<std::string> SlurpFile(const std::filesystem::path& p);

// --- predicate evaluation
bool MatchOne(const WindowFeatures& f, const FeaturePredicate& p);

// Apply mapping_hints JSON object to ParamOverride.
void DecodeMappingHints(const Json& obj, ParamOverride& p);

// Parse a feature_signature node in either array-of-predicates form or compact object form.
std::vector<FeaturePredicate> DecodeSignature(const Json& sig);

}  // namespace internal
}  // namespace llmap::annot
