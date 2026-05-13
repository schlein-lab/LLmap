#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace llmap::self_interference {

// Search result for a single query vector
struct KnnResult {
    std::vector<int64_t> indices;     // k nearest neighbor indices
    std::vector<float> distances;     // corresponding distances (L2 or inner product)
};

// Batch search results
struct BatchKnnResult {
    std::vector<int64_t> indices;     // Shape: [num_queries * k]
    std::vector<float> distances;     // Shape: [num_queries * k]
    size_t num_queries = 0;
    size_t k = 0;
    float search_time_ms = 0.0f;
};

// FAISS index type
enum class FaissIndexType {
    FlatL2,           // Exact search with L2 distance (brute force)
    FlatIP,           // Exact search with inner product
    IVFFlat,          // Inverted file index with L2 (approximate, faster)
    IVFPQ,            // Inverted file + product quantization (smaller, faster, less accurate)
};

// FAISS execution provider
enum class FaissProvider {
    CPU,              // CPU-only FAISS
    GPU,              // FAISS-GPU (requires CUDA)
};

// Configuration for the FAISS index
struct FaissIndexConfig {
    FaissIndexType index_type = FaissIndexType::IVFFlat;
    FaissProvider provider = FaissProvider::CPU;
    int gpu_device_id = 0;              // GPU device when using GPU provider

    size_t embedding_dim = 256;         // Dimension of input vectors

    // IVF-specific parameters
    size_t nlist = 1024;                // Number of Voronoi cells for IVF
    size_t nprobe = 32;                 // Number of cells to visit during search (accuracy vs speed)

    // PQ-specific parameters (for IVFPQ)
    size_t m_subquantizers = 8;         // Number of sub-quantizers
    size_t nbits_per_idx = 8;           // Bits per sub-quantizer index

    // Training parameters
    size_t training_sample_size = 0;    // 0 = use all vectors; otherwise subsample for training

    // Search parameters
    size_t default_k = 50;              // Default number of neighbors

    bool normalize_vectors = true;      // L2-normalize before indexing (for cosine similarity)
};

// Statistics about the index
struct FaissIndexStats {
    size_t num_vectors = 0;
    size_t embedding_dim = 0;
    size_t memory_usage_bytes = 0;
    bool is_trained = false;
    float training_time_ms = 0.0f;
    float add_time_ms = 0.0f;
};

// Forward declaration for PIMPL
class FaissIndexImpl;

// FAISS index wrapper for ANN search
// Primary use case: Stage 1 self-interference (read-vs-read similarity)
class FaissIndex {
public:
    // Factory: create index from config
    // Returns nullptr if initialization fails (e.g., FAISS-GPU not available)
    static std::unique_ptr<FaissIndex> Create(const FaissIndexConfig& config);

    ~FaissIndex();

    // Non-copyable, movable
    FaissIndex(const FaissIndex&) = delete;
    FaissIndex& operator=(const FaissIndex&) = delete;
    FaissIndex(FaissIndex&&) noexcept;
    FaissIndex& operator=(FaissIndex&&) noexcept;

    // ========== Index construction ==========

    // Train the index (required for IVF-based indices before adding)
    // vectors: [n x dim] row-major float array
    bool Train(std::span<const float> vectors, size_t num_vectors);

    // Add vectors to the index
    // vectors: [n x dim] row-major float array
    // Returns number of vectors successfully added
    size_t Add(std::span<const float> vectors, size_t num_vectors);

    // Add vectors with explicit IDs
    // ids must have num_vectors elements
    size_t AddWithIds(
        std::span<const float> vectors,
        std::span<const int64_t> ids,
        size_t num_vectors);

    // Train and add in one call (convenience)
    bool TrainAndAdd(std::span<const float> vectors, size_t num_vectors);

    // ========== Search ==========

    // Search for k nearest neighbors of a single query
    std::optional<KnnResult> Search(std::span<const float> query, size_t k) const;

    // Batch search: find k nearest neighbors for multiple queries
    // queries: [n x dim] row-major float array
    std::optional<BatchKnnResult> SearchBatch(
        std::span<const float> queries,
        size_t num_queries,
        size_t k) const;

    // Search with pre-allocated output buffers
    // indices: [num_queries * k], distances: [num_queries * k]
    // Returns false on error
    bool SearchBatchInto(
        std::span<const float> queries,
        size_t num_queries,
        size_t k,
        std::span<int64_t> indices,
        std::span<float> distances) const;

    // ========== Configuration ==========

    // Query configuration
    const FaissIndexConfig& Config() const { return config_; }
    size_t EmbeddingDim() const { return config_.embedding_dim; }
    size_t DefaultK() const { return config_.default_k; }

    // Set number of cells to probe (affects accuracy vs speed tradeoff)
    void SetNprobe(size_t nprobe);
    size_t GetNprobe() const;

    // ========== Status ==========

    bool IsReady() const;
    bool IsTrained() const;
    size_t NumVectors() const;
    FaissIndexStats GetStats() const;

    // Check if GPU is being used
    bool IsGpuEnabled() const;
    FaissProvider ActiveProvider() const { return active_provider_; }

    // Last error message for debugging
    std::string LastError() const;

    // ========== Serialization ==========

    // Save index to file
    bool Save(const std::string& path) const;

    // Load index from file (static factory)
    static std::unique_ptr<FaissIndex> Load(
        const std::string& path,
        const FaissIndexConfig& config);

private:
    explicit FaissIndex(const FaissIndexConfig& config);

    bool Initialize();

    // L2-normalize vectors in place
    void NormalizeVectors(std::span<float> vectors, size_t num_vectors) const;

    // Normalize a copy of input vectors (for const methods)
    std::vector<float> NormalizeVectorsCopy(
        std::span<const float> vectors,
        size_t num_vectors) const;

    FaissIndexConfig config_;
    FaissProvider active_provider_;
    std::unique_ptr<FaissIndexImpl> impl_;
};

// ========== Utilities ==========

// Check if FAISS is available at runtime
bool IsFaissAvailable();

// Check if FAISS-GPU is available
bool IsFaissGpuAvailable();

// Get FAISS version string
std::string GetFaissVersion();

// Compute recall@k between approximate and exact search results
// Returns fraction of true k-NN found by approximate search
float ComputeRecallAtK(
    std::span<const int64_t> approx_indices,
    std::span<const int64_t> exact_indices,
    size_t num_queries,
    size_t k);

}  // namespace llmap::self_interference
