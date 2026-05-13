#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::ai {

// Execution provider for ONNX Runtime
enum class ExecutionProvider {
    CPU,       // Always available
    CUDA,      // Requires CUDA toolkit + ONNX CUDA EP
    TensorRT,  // Requires TensorRT + ONNX TRT EP
};

// Configuration for the embedder
struct EmbedderConfig {
    std::filesystem::path model_path;       // Path to .onnx file
    ExecutionProvider provider = ExecutionProvider::CPU;
    int device_id = 0;                       // GPU device ID when using CUDA/TensorRT
    size_t batch_size = 1024;                // Max reads per batch
    size_t embedding_dim = 256;              // Output embedding dimension
    size_t max_sequence_length = 512;        // Max nucleotides per read
    bool enable_memory_pattern = true;       // ONNX Runtime memory optimization
    int intra_op_threads = 0;                // 0 = auto-detect
    int inter_op_threads = 0;                // 0 = auto-detect
};

// Single embedding result
struct EmbeddingResult {
    std::vector<float> embedding;  // Shape: [embedding_dim]
    float inference_time_us;       // Microseconds for this inference
};

// Batch embedding result
struct BatchEmbeddingResult {
    std::vector<std::vector<float>> embeddings;  // Shape: [batch_size, embedding_dim]
    float total_inference_time_us;                // Total inference time
    float avg_time_per_read_us;                   // Average per read
};

// Forward declarations for ONNX Runtime (PIMPL pattern to hide ORT headers)
class OnnxRuntimeSession;

// Foundation model embedder for DNA sequences
// Wraps Caduceus-Ph distilled model (or fallback DNABERT-2 on CPU)
class FoundationEmbedder {
public:
    // Factory: create embedder from config
    // Returns nullptr if initialization fails (e.g., missing ONNX Runtime, bad model path)
    static std::unique_ptr<FoundationEmbedder> Create(const EmbedderConfig& config);

    ~FoundationEmbedder();

    // Non-copyable, movable
    FoundationEmbedder(const FoundationEmbedder&) = delete;
    FoundationEmbedder& operator=(const FoundationEmbedder&) = delete;
    FoundationEmbedder(FoundationEmbedder&&) noexcept;
    FoundationEmbedder& operator=(FoundationEmbedder&&) noexcept;

    // Embed a single DNA sequence
    // Input: nucleotide string (A, C, G, T, N allowed)
    // Returns nullopt on error (e.g., sequence too long, invalid characters)
    std::optional<EmbeddingResult> Embed(std::string_view sequence) const;

    // Embed a batch of sequences
    // Sequences are padded/truncated to max_sequence_length
    std::optional<BatchEmbeddingResult> EmbedBatch(
        std::span<const std::string_view> sequences) const;

    // Embed batch with pre-allocated output buffer
    // Output shape: [sequences.size(), embedding_dim]
    // Returns false on error
    bool EmbedBatchInto(
        std::span<const std::string_view> sequences,
        std::span<float> output) const;

    // Query configuration
    const EmbedderConfig& Config() const { return config_; }
    size_t EmbeddingDim() const { return config_.embedding_dim; }
    size_t MaxSequenceLength() const { return config_.max_sequence_length; }
    bool IsGpuEnabled() const;
    ExecutionProvider ActiveProvider() const { return active_provider_; }

    // Model info
    std::string ModelName() const;
    std::string ModelVersion() const;

    // Status check
    bool IsReady() const { return session_ != nullptr; }

    // Get last error message (for debugging initialization failures)
    std::string LastError() const;

private:
    explicit FoundationEmbedder(const EmbedderConfig& config);

    // Initialize ONNX Runtime session
    bool Initialize();

    // Tokenize DNA sequence to model input format
    // Returns token IDs suitable for the model
    std::vector<int64_t> Tokenize(std::string_view sequence) const;

    // Batch tokenization with padding
    void TokenizeBatch(
        std::span<const std::string_view> sequences,
        std::vector<int64_t>& token_ids,
        std::vector<int64_t>& attention_mask) const;

    EmbedderConfig config_;
    ExecutionProvider active_provider_;  // May differ from config if fallback occurred
    std::unique_ptr<OnnxRuntimeSession> session_;
};

// Utility: check if ONNX Runtime is available at runtime
bool IsOnnxRuntimeAvailable();

// Utility: check if CUDA execution provider is available
bool IsCudaEPAvailable();

// Utility: check if TensorRT execution provider is available
bool IsTensorRTEPAvailable();

// Utility: list available execution providers
std::vector<std::string> ListAvailableProviders();

}  // namespace llmap::ai
