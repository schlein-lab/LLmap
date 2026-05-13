#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "ai/foundation_embedder.h"

namespace llmap::ai {

// Bucket embedder configuration
// Tuned for longer reference sequences (bucket-level, 10kb-5MB windows)
struct BucketEmbedderConfig {
    std::filesystem::path model_path;       // Path to .onnx file (Evo-1.5B distilled)
    ExecutionProvider provider = ExecutionProvider::CPU;
    int device_id = 0;                       // GPU device ID when using CUDA/TensorRT
    size_t batch_size = 64;                  // Max buckets per batch (larger than reads)
    size_t embedding_dim = 256;              // Output embedding dimension
    size_t max_sequence_length = 8192;       // Max nucleotides per bucket (8kb default)
    size_t chunk_size = 1024;                // Sub-chunk for very long sequences
    size_t chunk_overlap = 128;              // Overlap between chunks for continuity
    bool enable_memory_pattern = true;       // ONNX Runtime memory optimization
    int intra_op_threads = 0;                // 0 = auto-detect
    int inter_op_threads = 0;                // 0 = auto-detect
    bool normalize_output = true;            // L2-normalize embeddings
};

// Statistics from embedding a set of buckets
struct BucketEmbeddingStats {
    size_t total_buckets = 0;
    size_t total_nucleotides = 0;
    size_t total_chunks = 0;                 // For long sequences that required chunking
    float total_time_ms = 0.0f;
    float avg_time_per_bucket_ms = 0.0f;
    float throughput_buckets_per_sec = 0.0f;
    float throughput_mb_per_sec = 0.0f;      // Megabases per second
};

// Forward declarations
class OnnxRuntimeSession;

// Bucket embedder for reference genome windows
// Wraps Evo-1.5B distilled model (or similar long-context DNA model)
// Used at index-time to pre-compute bucket embeddings for π_AI(b|r)
class BucketEmbedder {
public:
    // Factory: create embedder from config
    // Returns nullptr if initialization fails
    static std::unique_ptr<BucketEmbedder> Create(const BucketEmbedderConfig& config);

    ~BucketEmbedder();

    // Non-copyable, movable
    BucketEmbedder(const BucketEmbedder&) = delete;
    BucketEmbedder& operator=(const BucketEmbedder&) = delete;
    BucketEmbedder(BucketEmbedder&&) noexcept;
    BucketEmbedder& operator=(BucketEmbedder&&) noexcept;

    // Embed a single bucket sequence
    // For long sequences (> max_sequence_length), chunks and pools
    // Returns nullopt on error
    std::optional<std::vector<float>> EmbedBucket(std::string_view sequence) const;

    // Embed multiple bucket sequences
    // Returns matrix [num_buckets x embedding_dim]
    std::optional<std::vector<std::vector<float>>> EmbedBuckets(
        std::span<const std::string_view> sequences) const;

    // Embed multiple buckets with pre-allocated output
    // Output shape: [sequences.size(), embedding_dim]
    // Returns false on error, also fills stats if provided
    bool EmbedBucketsInto(
        std::span<const std::string_view> sequences,
        std::span<float> output,
        BucketEmbeddingStats* stats = nullptr) const;

    // Query configuration
    const BucketEmbedderConfig& Config() const { return config_; }
    size_t EmbeddingDim() const { return config_.embedding_dim; }
    size_t MaxSequenceLength() const { return config_.max_sequence_length; }
    bool IsGpuEnabled() const;
    ExecutionProvider ActiveProvider() const { return active_provider_; }

    // Model info
    std::string ModelName() const;
    std::string ModelVersion() const;

    // Status
    bool IsReady() const { return session_ != nullptr; }
    std::string LastError() const;

private:
    explicit BucketEmbedder(const BucketEmbedderConfig& config);

    bool Initialize();

    // Tokenize DNA sequence
    std::vector<int64_t> Tokenize(std::string_view sequence) const;

    // Batch tokenization with padding
    void TokenizeBatch(
        std::span<const std::string_view> sequences,
        std::vector<int64_t>& token_ids,
        std::vector<int64_t>& attention_mask) const;

    // For long sequences: chunk, embed each, then pool
    std::optional<std::vector<float>> EmbedLongSequence(std::string_view sequence) const;

    // Pool multiple chunk embeddings into one
    // Strategy: mean pooling with optional attention weighting
    std::vector<float> PoolChunkEmbeddings(
        const std::vector<std::vector<float>>& chunk_embeddings) const;

    // L2-normalize embedding in place
    void NormalizeEmbedding(std::span<float> embedding) const;

    BucketEmbedderConfig config_;
    ExecutionProvider active_provider_;
    std::unique_ptr<OnnxRuntimeSession> session_;
};

// Utility: compute cosine similarity between two embeddings
float CosineSimilarity(std::span<const float> a, std::span<const float> b);

// Utility: compute pairwise cosine similarity matrix
// Output: [n x n] matrix where output[i*n + j] = cosine_sim(embeddings[i], embeddings[j])
void ComputeCosineSimilarityMatrix(
    std::span<const float> embeddings,
    size_t num_vectors,
    size_t embedding_dim,
    std::span<float> output);

}  // namespace llmap::ai
