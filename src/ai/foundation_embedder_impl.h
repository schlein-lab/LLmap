#pragma once

#include "ai/foundation_embedder.h"

#include <memory>
#include <string>
#include <vector>

// Conditionally include ONNX Runtime headers based on availability
#ifdef LLMAP_HAS_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

namespace llmap::ai {

// PIMPL implementation class that holds ONNX Runtime state
class OnnxRuntimeSession {
public:
#ifdef LLMAP_HAS_ONNXRUNTIME
    Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "llmap"};
    std::unique_ptr<Ort::Session> session;
    Ort::AllocatorWithDefaultOptions allocator;
    std::vector<const char*> input_names;
    std::vector<const char*> output_names;
    std::vector<std::string> input_name_strings;
    std::vector<std::string> output_name_strings;
#endif
    bool initialized = false;
    std::string error_message;
};

// DNA nucleotide to token ID mapping (simplified tokenization)
// Real Caduceus-Ph uses BPE, but for now we use character-level
constexpr int64_t kPadToken = 0;
constexpr int64_t kUnkToken = 1;
constexpr int64_t kClsToken = 2;
constexpr int64_t kSepToken = 3;
constexpr int64_t kAToken = 4;
constexpr int64_t kCToken = 5;
constexpr int64_t kGToken = 6;
constexpr int64_t kTToken = 7;
constexpr int64_t kNToken = 8;

int64_t NucleotideToToken(char c);
bool IsValidNucleotide(char c);

}  // namespace llmap::ai
