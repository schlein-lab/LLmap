#include "ai/bucket_embedder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>

#ifdef LLMAP_HAS_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

namespace llmap::ai {

namespace {

// DNA nucleotide to token ID mapping
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
// (duplicate of foundation_embedder.cpp - internal implementation detail)
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

std::vector<int64_t> BucketEmbedder::Tokenize(std::string_view sequence) const {
    std::vector<int64_t> tokens;
    size_t max_len = config_.max_sequence_length;
    tokens.reserve(max_len + 2);

    tokens.push_back(kClsToken);

    size_t seq_len = std::min(sequence.length(), max_len);
    for (size_t i = 0; i < seq_len; ++i) {
        tokens.push_back(NucleotideToToken(sequence[i]));
    }

    tokens.push_back(kSepToken);

    while (tokens.size() < max_len + 2) {
        tokens.push_back(kPadToken);
    }

    return tokens;
}

void BucketEmbedder::TokenizeBatch(
    std::span<const std::string_view> sequences,
    std::vector<int64_t>& token_ids,
    std::vector<int64_t>& attention_mask) const {

    size_t batch_size = sequences.size();
    size_t seq_len = config_.max_sequence_length + 2;

    token_ids.resize(batch_size * seq_len);
    attention_mask.resize(batch_size * seq_len);

    for (size_t b = 0; b < batch_size; ++b) {
        size_t offset = b * seq_len;
        const auto& seq = sequences[b];

        token_ids[offset] = kClsToken;
        attention_mask[offset] = 1;

        size_t actual_len = std::min(seq.length(), config_.max_sequence_length);
        for (size_t i = 0; i < actual_len; ++i) {
            token_ids[offset + 1 + i] = NucleotideToToken(seq[i]);
            attention_mask[offset + 1 + i] = 1;
        }

        token_ids[offset + 1 + actual_len] = kSepToken;
        attention_mask[offset + 1 + actual_len] = 1;

        for (size_t i = actual_len + 2; i < seq_len; ++i) {
            token_ids[offset + i] = kPadToken;
            attention_mask[offset + i] = 0;
        }
    }
}

void BucketEmbedder::NormalizeEmbedding(std::span<float> embedding) const {
    float norm_sq = 0.0f;
    for (float v : embedding) {
        norm_sq += v * v;
    }
    if (norm_sq > 1e-12f) {
        float inv_norm = 1.0f / std::sqrt(norm_sq);
        for (float& v : embedding) {
            v *= inv_norm;
        }
    }
}

std::vector<float> BucketEmbedder::PoolChunkEmbeddings(
    const std::vector<std::vector<float>>& chunk_embeddings) const {

    if (chunk_embeddings.empty()) {
        return std::vector<float>(config_.embedding_dim, 0.0f);
    }

    if (chunk_embeddings.size() == 1) {
        return chunk_embeddings[0];
    }

    size_t dim = chunk_embeddings[0].size();
    std::vector<float> pooled(dim, 0.0f);

    for (const auto& chunk : chunk_embeddings) {
        for (size_t i = 0; i < dim; ++i) {
            pooled[i] += chunk[i];
        }
    }

    float inv_n = 1.0f / static_cast<float>(chunk_embeddings.size());
    for (float& v : pooled) {
        v *= inv_n;
    }

    if (config_.normalize_output) {
        NormalizeEmbedding(pooled);
    }

    return pooled;
}

std::optional<std::vector<float>> BucketEmbedder::EmbedLongSequence(
    std::string_view sequence) const {

    if (!session_->initialized) {
        return std::nullopt;
    }

    size_t chunk_size = config_.chunk_size;
    size_t overlap = config_.chunk_overlap;
    size_t step = chunk_size - overlap;

    std::vector<std::string_view> chunks;
    for (size_t start = 0; start < sequence.length(); start += step) {
        size_t len = std::min(chunk_size, sequence.length() - start);
        chunks.push_back(sequence.substr(start, len));
        if (start + len >= sequence.length()) break;
    }

    if (chunks.empty()) {
        return std::nullopt;
    }

    auto chunk_results = EmbedBuckets(chunks);
    if (!chunk_results) {
        return std::nullopt;
    }

    return PoolChunkEmbeddings(*chunk_results);
}

std::optional<std::vector<float>> BucketEmbedder::EmbedBucket(
    std::string_view sequence) const {

    if (!session_->initialized) {
        return std::nullopt;
    }

    for (char c : sequence) {
        if (!IsValidNucleotide(c)) {
            return std::nullopt;
        }
    }

    if (sequence.length() > config_.max_sequence_length) {
        return EmbedLongSequence(sequence);
    }

    std::array<std::string_view, 1> seqs = {sequence};
    auto batch_result = EmbedBuckets(seqs);
    if (!batch_result || batch_result->empty()) {
        return std::nullopt;
    }

    return std::move((*batch_result)[0]);
}

std::optional<std::vector<std::vector<float>>> BucketEmbedder::EmbedBuckets(
    std::span<const std::string_view> sequences) const {

    if (!session_->initialized || sequences.empty()) {
        return std::nullopt;
    }

#ifdef LLMAP_HAS_ONNXRUNTIME
    try {
        std::vector<int64_t> token_ids;
        std::vector<int64_t> attention_mask;
        TokenizeBatch(sequences, token_ids, attention_mask);

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

        auto output_tensors = session_->session->Run(
            Ort::RunOptions{nullptr},
            session_->input_names.data(),
            input_tensors.data(),
            input_tensors.size(),
            session_->output_names.data(),
            session_->output_names.size());

        auto& output_tensor = output_tensors[0];
        auto output_shape = output_tensor.GetTensorTypeAndShapeInfo().GetShape();
        const float* output_data = output_tensor.GetTensorData<float>();

        size_t embedding_dim = config_.embedding_dim;
        if (output_shape.size() >= 2) {
            embedding_dim = static_cast<size_t>(output_shape[1]);
        }

        std::vector<std::vector<float>> results;
        results.reserve(batch_size);

        for (size_t b = 0; b < batch_size; ++b) {
            std::vector<float> embedding(embedding_dim);
            std::copy_n(output_data + b * embedding_dim, embedding_dim, embedding.begin());

            if (config_.normalize_output) {
                NormalizeEmbedding(embedding);
            }

            results.push_back(std::move(embedding));
        }

        return results;
    } catch (const Ort::Exception& e) {
        return std::nullopt;
    }
#else
    return std::nullopt;
#endif
}

bool BucketEmbedder::EmbedBucketsInto(
    std::span<const std::string_view> sequences,
    std::span<float> output,
    BucketEmbeddingStats* stats) const {

    auto start = std::chrono::high_resolution_clock::now();

    auto result = EmbedBuckets(sequences);
    if (!result) {
        return false;
    }

    size_t expected_size = sequences.size() * config_.embedding_dim;
    if (output.size() < expected_size) {
        return false;
    }

    size_t offset = 0;
    for (const auto& emb : *result) {
        std::copy(emb.begin(), emb.end(), output.begin() + offset);
        offset += emb.size();
    }

    auto end = std::chrono::high_resolution_clock::now();

    if (stats) {
        stats->total_buckets = sequences.size();
        stats->total_nucleotides = 0;
        for (const auto& seq : sequences) {
            stats->total_nucleotides += seq.size();
        }
        stats->total_chunks = 0;  // TODO: track actual chunks for long sequences
        stats->total_time_ms = std::chrono::duration<float, std::milli>(end - start).count();
        stats->avg_time_per_bucket_ms = stats->total_time_ms / static_cast<float>(sequences.size());

        if (stats->total_time_ms > 0.0f) {
            stats->throughput_buckets_per_sec =
                static_cast<float>(stats->total_buckets) / (stats->total_time_ms / 1000.0f);
            stats->throughput_mb_per_sec =
                static_cast<float>(stats->total_nucleotides) / (stats->total_time_ms / 1000.0f) / 1e6f;
        }
    }

    return true;
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
