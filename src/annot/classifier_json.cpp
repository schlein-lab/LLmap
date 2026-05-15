#include "annot/classifier_internal.h"

#include <cctype>
#include <fstream>
#include <sstream>

namespace llmap::annot::internal {

const Json* Json::Get(const std::string& key) const {
    if (type != Type::Object) return nullptr;
    for (const auto& kv : obj_v) if (kv.first == key) return &kv.second;
    return nullptr;
}

bool Json::AsBool(bool defv) const {
    if (type == Type::Bool) return bool_v;
    return defv;
}

double Json::AsNumber(double defv) const {
    if (type == Type::Number) return num_v;
    return defv;
}

std::string Json::AsString(const std::string& defv) const {
    if (type == Type::String) return str_v;
    return defv;
}

JsonReader::JsonReader(std::string src) : src_(std::move(src)), pos_(0) {}

std::optional<Json> JsonReader::Parse() {
    SkipWs();
    Json out;
    if (!ParseValue(out)) return std::nullopt;
    SkipWs();
    return out;
}

void JsonReader::SkipWs() {
    while (pos_ < src_.size() &&
           (std::isspace(static_cast<unsigned char>(src_[pos_])))) ++pos_;
}

bool JsonReader::Peek(char c) {
    SkipWs();
    return pos_ < src_.size() && src_[pos_] == c;
}

bool JsonReader::Eat(char c) {
    SkipWs();
    if (pos_ < src_.size() && src_[pos_] == c) { ++pos_; return true; }
    return false;
}

bool JsonReader::ParseValue(Json& out) {
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

bool JsonReader::ParseString(Json& out) {
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

bool JsonReader::ParseNumber(Json& out) {
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

bool JsonReader::ParseBool(Json& out) {
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

bool JsonReader::ParseNull(Json& out) {
    if (src_.compare(pos_, 4, "null") == 0) {
        pos_ += 4;
        out.type = Json::Type::Null;
        return true;
    }
    return false;
}

bool JsonReader::ParseArray(Json& out) {
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

bool JsonReader::ParseObject(Json& out) {
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

std::optional<std::string> SlurpFile(const std::filesystem::path& p) {
    std::ifstream f(p);
    if (!f) return std::nullopt;
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace llmap::annot::internal
