#include "ai/foundation_embedder.h"

#include <algorithm>
#include <chrono>
#include <stdexcept>

// Conditionally include ONNX Runtime headers based on availability
#ifdef LLMAP_HAS_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

namespace llmap::ai {

namespace {

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

int64_t NucleotideToToken(char c) {
    switch (c) {
        case 'A': case 'a': return kAToken;
        case 'C': case 'c': return kCToken;
        case 'G': case 'g': return kGToken;
        case 'T': case 't': return kTToken;
        case 'N': case 'n': return kNToken;
        default: return kUnkToken;
    }
}

bool IsValidNucleotide(char c) {
    switch (c) {
        case 'A': case 'a':
        case 'C': case 'c':
        case 'G': case 'g':
        case 'T': case 't':
        case 'N': case 'n':
            return true;
        default:
            return false;
    }
}

}  // namespace

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
    }
#else
    // No ONNX Runtime available — embedder is not functional but compiles
    session_->error_message = "ONNX Runtime not available at compile time";
    session_->initialized = false;
    return false;
#endif
}

std::vector<int64_t> FoundationEmbedder::Tokenize(std::string_view sequence) const {
    std::vector<int64_t> tokens;
    tokens.reserve(config_.max_sequence_length + 2);  // +2 for CLS and SEP

    // Add CLS token
    tokens.push_back(kClsToken);

    // Tokenize sequence
    size_t seq_len = std::min(sequence.length(), config_.max_sequence_length);
    for (size_t i = 0; i < seq_len; ++i) {
        tokens.push_back(NucleotideToToken(sequence[i]));
    }

    // Add SEP token
    tokens.push_back(kSepToken);

    // Pad to max length
    while (tokens.size() < config_.max_sequence_length + 2) {
        tokens.push_back(kPadToken);
    }

    return tokens;
}

void FoundationEmbedder::TokenizeBatch(
    std::span<const std::string_view> sequences,
    std::vector<int64_t>& token_ids,
    std::vector<int64_t>& attention_mask) const {

    size_t batch_size = sequences.size();
    size_t seq_len = config_.max_sequence_length + 2;  // +2 for CLS and SEP

    token_ids.resize(batch_size * seq_len);
    attention_mask.resize(batch_size * seq_len);

    for (size_t b = 0; b < batch_size; ++b) {
        size_t offset = b * seq_len;
        const auto& seq = sequences[b];

        // CLS token
        token_ids[offset] = kClsToken;
        attention_mask[offset] = 1;

        // Sequence tokens
        size_t actual_len = std::min(seq.length(), config_.max_sequence_length);
        for (size_t i = 0; i < actual_len; ++i) {
            token_ids[offset + 1 + i] = NucleotideToToken(seq[i]);
            attention_mask[offset + 1 + i] = 1;
        }

        // SEP token
        token_ids[offset + 1 + actual_len] = kSepToken;
        attention_mask[offset + 1 + actual_len] = 1;

        // Padding
        for (size_t i = actual_len + 2; i < seq_len; ++i) {
            token_ids[offset + i] = kPadToken;
            attention_mask[offset + i] = 0;
        }
    }
}

std::optional<EmbeddingResult> FoundationEmbedder::Embed(std::string_view sequence) const {
    if (!session_->initialized) {
        return std::nullopt;
    }

    // Validate sequence
    for (char c : sequence) {
        if (!IsValidNucleotide(c)) {
            return std::nullopt;
        }
    }

    if (sequence.length() > config_.max_sequence_length) {
        sequence = sequence.substr(0, config_.max_sequence_length);
    }

    std::array<std::string_view, 1> seqs = {sequence};
    auto batch_result = EmbedBatch(seqs);
    if (!batch_result) {
        return std::nullopt;
    }

    return EmbeddingResult{
        .embedding = std::move(batch_result->embeddings[0]),
        .inference_time_us = batch_result->total_inference_time_us,
    };
}

std::optional<BatchEmbeddingResult> FoundationEmbedder::EmbedBatch(
    std::span<const std::string_view> sequences) const {

    if (!session_->initialized || sequences.empty()) {
        return std::nullopt;
    }

#ifdef LLMAP_HAS_ONNXRUNTIME
    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Tokenize batch
        std::vector<int64_t> token_ids;
        std::vector<int64_t> attention_mask;
        TokenizeBatch(sequences, token_ids, attention_mask);

        // Create input tensors
        size_t batch_size = sequences.size();
        size_t seq_len = config_.max_sequence_length + 2;

        std::array<int64_t, 2> input_shape = {
            static_cast<int64_t>(batch_size),
            static_cast<int64_t>(seq_len)
        };

        Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(
            OrtArenaAllocator, OrtMemTypeDefault);

        std::vector<Ort::Value> input_tensors;
        input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
            memory_info, token_ids.data(), token_ids.size(),
            input_shape.data(), input_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
            memory_info, attention_mask.data(), attention_mask.size(),
            input_shape.data(), input_shape.size()));

        // Run inference
        auto output_tensors = session_->session->Run(
            Ort::RunOptions{nullptr},
            session_->input_names.data(),
            input_tensors.data(),
            input_tensors.size(),
            session_->output_names.data(),
            session_->output_names.size());

        auto end = std::chrono::high_resolution_clock::now();
        float total_time_us = std::chrono::duration<float, std::micro>(end - start).count();

        // Extract embeddings from output tensor
        // Expected shape: [batch_size, embedding_dim]
        auto& output_tensor = output_tensors[0];
        auto output_shape = output_tensor.GetTensorTypeAndShapeInfo().GetShape();
        const float* output_data = output_tensor.GetTensorData<float>();

        BatchEmbeddingResult result;
        result.total_inference_time_us = total_time_us;
        result.avg_time_per_read_us = total_time_us / static_cast<float>(batch_size);
        result.embeddings.reserve(batch_size);

        size_t embedding_dim = config_.embedding_dim;
        if (output_shape.size() >= 2) {
            embedding_dim = static_cast<size_t>(output_shape[1]);
        }

        for (size_t b = 0; b < batch_size; ++b) {
            std::vector<float> embedding(embedding_dim);
            std::copy_n(output_data + b * embedding_dim, embedding_dim, embedding.begin());
            result.embeddings.push_back(std::move(embedding));
        }

        return result;
    } catch (const Ort::Exception& e) {
        return std::nullopt;
    }
#else
    return std::nullopt;
#endif
}

bool FoundationEmbedder::EmbedBatchInto(
    std::span<const std::string_view> sequences,
    std::span<float> output) const {

    auto result = EmbedBatch(sequences);
    if (!result) {
        return false;
    }

    size_t expected_size = sequences.size() * config_.embedding_dim;
    if (output.size() < expected_size) {
        return false;
    }

    size_t offset = 0;
    for (const auto& emb : result->embeddings) {
        std::copy(emb.begin(), emb.end(), output.begin() + offset);
        offset += emb.size();
    }

    return true;
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
