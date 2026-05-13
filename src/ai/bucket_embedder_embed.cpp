#include "ai/bucket_embedder.h"
#include "ai/bucket_embedder_impl.h"

#include <algorithm>
#include <chrono>

namespace llmap::ai {

namespace {

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

std::optional<std::vector<float>> BucketEmbedder::EmbedLongSequence(
    std::string_view sequence) const {

    if (!session_ || !session_->initialized) {
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

    if (!session_ || !session_->initialized) {
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

    if (!session_ || !session_->initialized || sequences.empty()) {
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
        stats->total_chunks = 0;
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

}  // namespace llmap::ai
