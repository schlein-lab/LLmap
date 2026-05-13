#pragma once

// Internal header for BucketEmbedder implementation details
// Not part of the public API - only used by bucket_embedder_*.cpp files

#ifdef LLMAP_HAS_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

#include <memory>
#include <string>
#include <vector>

namespace llmap::ai {

// PIMPL implementation class that holds ONNX Runtime state
class OnnxRuntimeSession {
public:
#ifdef LLMAP_HAS_ONNXRUNTIME
    Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "llmap_bucket"};
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

}  // namespace llmap::ai
