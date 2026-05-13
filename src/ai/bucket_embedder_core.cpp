#include "ai/bucket_embedder.h"
#include "ai/bucket_embedder_impl.h"

#include <cmath>
#include <stdexcept>

namespace llmap::ai {

BucketEmbedder::BucketEmbedder(const BucketEmbedderConfig& config)
    : config_(config),
      active_provider_(ExecutionProvider::CPU),
      session_(std::make_unique<OnnxRuntimeSession>()) {}

BucketEmbedder::~BucketEmbedder() = default;

BucketEmbedder::BucketEmbedder(BucketEmbedder&&) noexcept = default;
BucketEmbedder& BucketEmbedder::operator=(BucketEmbedder&&) noexcept = default;

std::unique_ptr<BucketEmbedder> BucketEmbedder::Create(const BucketEmbedderConfig& config) {
    auto embedder = std::unique_ptr<BucketEmbedder>(new BucketEmbedder(config));
    if (!embedder->Initialize()) {
        return nullptr;
    }
    return embedder;
}

bool BucketEmbedder::Initialize() {
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

        if (!std::filesystem::exists(config_.model_path)) {
            session_->error_message = "Model file not found: " + config_.model_path.string();
            session_->initialized = false;
            return false;
        }

        session_->session = std::make_unique<Ort::Session>(
            session_->env, config_.model_path.c_str(), session_options);

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
    session_->error_message = "ONNX Runtime not available at compile time";
    session_->initialized = false;
    return false;
#endif
}

bool BucketEmbedder::IsGpuEnabled() const {
    return active_provider_ == ExecutionProvider::CUDA ||
           active_provider_ == ExecutionProvider::TensorRT;
}

std::string BucketEmbedder::ModelName() const {
    return config_.model_path.stem().string();
}

std::string BucketEmbedder::ModelVersion() const {
    return "1.0";
}

std::string BucketEmbedder::LastError() const {
    if (session_) {
        return session_->error_message;
    }
    return "Session not initialized";
}

// Utility functions

float CosineSimilarity(std::span<const float> a, std::span<const float> b) {
    if (a.size() != b.size() || a.empty()) {
        return 0.0f;
    }

    float dot = 0.0f;
    float norm_a = 0.0f;
    float norm_b = 0.0f;

    for (size_t i = 0; i < a.size(); ++i) {
        dot += a[i] * b[i];
        norm_a += a[i] * a[i];
        norm_b += b[i] * b[i];
    }

    float denom = std::sqrt(norm_a) * std::sqrt(norm_b);
    if (denom < 1e-12f) {
        return 0.0f;
    }

    return dot / denom;
}

void ComputeCosineSimilarityMatrix(
    std::span<const float> embeddings,
    size_t num_vectors,
    size_t embedding_dim,
    std::span<float> output) {

    if (output.size() < num_vectors * num_vectors) {
        return;
    }

    for (size_t i = 0; i < num_vectors; ++i) {
        for (size_t j = 0; j < num_vectors; ++j) {
            std::span<const float> vec_i(
                embeddings.data() + i * embedding_dim, embedding_dim);
            std::span<const float> vec_j(
                embeddings.data() + j * embedding_dim, embedding_dim);
            output[i * num_vectors + j] = CosineSimilarity(vec_i, vec_j);
        }
    }
}

}  // namespace llmap::ai
