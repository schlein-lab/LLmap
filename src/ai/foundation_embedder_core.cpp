#include "ai/foundation_embedder.h"
#include "ai/foundation_embedder_impl.h"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

namespace llmap::ai {

FoundationEmbedder::FoundationEmbedder(const EmbedderConfig& config)
    : config_(config),
      active_provider_(ExecutionProvider::CPU),
      session_(std::make_unique<OnnxRuntimeSession>()) {}

FoundationEmbedder::~FoundationEmbedder() = default;

FoundationEmbedder::FoundationEmbedder(FoundationEmbedder&&) noexcept = default;
FoundationEmbedder& FoundationEmbedder::operator=(FoundationEmbedder&&) noexcept = default;

std::unique_ptr<FoundationEmbedder> FoundationEmbedder::Create(const EmbedderConfig& config) {
    auto embedder = std::unique_ptr<FoundationEmbedder>(new FoundationEmbedder(config));
    if (!embedder->Initialize()) {
        return nullptr;
    }
    return embedder;
}

bool FoundationEmbedder::Initialize() {
#ifdef LLMAP_HAS_ONNXRUNTIME
    try {
        Ort::SessionOptions session_options;
        session_options.SetIntraOpNumThreads(
            config_.intra_op_threads > 0 ? config_.intra_op_threads : 0);
        session_options.SetInterOpNumThreads(
            config_.inter_op_threads > 0 ? config_.inter_op_threads : 0);

        if (config_.enable_memory_pattern) {
            session_options.EnableMemPattern();
        }

        session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

        // Try requested execution provider, fall back to CPU if unavailable
        active_provider_ = config_.provider;

        if (config_.provider == ExecutionProvider::CUDA) {
            try {
                OrtCUDAProviderOptions cuda_options;
                cuda_options.device_id = config_.device_id;
                session_options.AppendExecutionProvider_CUDA(cuda_options);
            } catch (const Ort::Exception& e) {
                active_provider_ = ExecutionProvider::CPU;
            }
        } else if (config_.provider == ExecutionProvider::TensorRT) {
            try {
                OrtTensorRTProviderOptions trt_options;
                trt_options.device_id = config_.device_id;
                session_options.AppendExecutionProvider_TensorRT(trt_options);
            } catch (const Ort::Exception& e) {
                active_provider_ = ExecutionProvider::CPU;
            }
        }

        // Check if model file exists
        if (!std::filesystem::exists(config_.model_path)) {
            session_->error_message = "Model file not found: " + config_.model_path.string();
            session_->initialized = false;
            return false;
        }

        session_->session = std::make_unique<Ort::Session>(
            session_->env, config_.model_path.c_str(), session_options);

        // Get input/output names
        size_t num_inputs = session_->session->GetInputCount();
        for (size_t i = 0; i < num_inputs; ++i) {
            auto name = session_->session->GetInputNameAllocated(i, session_->allocator);
            session_->input_name_strings.push_back(name.get());
        }
        for (const auto& name : session_->input_name_strings) {
            session_->input_names.push_back(name.c_str());
        }

        size_t num_outputs = session_->session->GetOutputCount();
        for (size_t i = 0; i < num_outputs; ++i) {
            auto name = session_->session->GetOutputNameAllocated(i, session_->allocator);
            session_->output_name_strings.push_back(name.get());
        }
        for (const auto& name : session_->output_name_strings) {
            session_->output_names.push_back(name.c_str());
        }

        session_->initialized = true;
        return true;
    } catch (const Ort::Exception& e) {
        session_->error_message = e.what();
        session_->initialized = false;
        return false;
    } catch (const std::exception& e) {
        session_->error_message = e.what();
        session_->initialized = false;
        return false;
    }
#else
    // No ONNX Runtime available — embedder is not functional but compiles
    session_->error_message = "ONNX Runtime not available at compile time";
    session_->initialized = false;
    return false;
#endif
}

bool FoundationEmbedder::IsGpuEnabled() const {
    return active_provider_ == ExecutionProvider::CUDA ||
           active_provider_ == ExecutionProvider::TensorRT;
}

std::string FoundationEmbedder::ModelName() const {
    return config_.model_path.stem().string();
}

std::string FoundationEmbedder::ModelVersion() const {
    return "1.0";  // Placeholder — would extract from ONNX model metadata
}

std::string FoundationEmbedder::LastError() const {
    if (session_) {
        return session_->error_message;
    }
    return "Session not initialized";
}

// Global utility functions

bool IsOnnxRuntimeAvailable() {
#ifdef LLMAP_HAS_ONNXRUNTIME
    return true;
#else
    return false;
#endif
}

bool IsCudaEPAvailable() {
#ifdef LLMAP_HAS_ONNXRUNTIME
    auto providers = Ort::GetAvailableProviders();
    return std::find(providers.begin(), providers.end(), "CUDAExecutionProvider") != providers.end();
#else
    return false;
#endif
}

bool IsTensorRTEPAvailable() {
#ifdef LLMAP_HAS_ONNXRUNTIME
    auto providers = Ort::GetAvailableProviders();
    return std::find(providers.begin(), providers.end(), "TensorrtExecutionProvider") != providers.end();
#else
    return false;
#endif
}

std::vector<std::string> ListAvailableProviders() {
#ifdef LLMAP_HAS_ONNXRUNTIME
    return Ort::GetAvailableProviders();
#else
    return {"CPUExecutionProvider (stub)"};
#endif
}

}  // namespace llmap::ai
