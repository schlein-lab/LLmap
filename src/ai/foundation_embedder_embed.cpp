#include "ai/foundation_embedder.h"
#include "ai/foundation_embedder_impl.h"

#include <algorithm>
#include <chrono>

namespace llmap::ai {

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

}  // namespace llmap::ai
