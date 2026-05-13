#include "self_interference/faiss_wrapper.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>

#ifdef LLMAP_HAS_FAISS
#include <faiss/Index.h>
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/IndexIVFPQ.h>
#include <faiss/index_io.h>
#ifdef LLMAP_HAS_FAISS_GPU
#include <faiss/gpu/GpuIndexFlat.h>
#include <faiss/gpu/GpuIndexIVFFlat.h>
#include <faiss/gpu/GpuResources.h>
#include <faiss/gpu/StandardGpuResources.h>
#endif
#endif

namespace llmap::self_interference {

// PIMPL implementation holding the actual FAISS index
class FaissIndexImpl {
public:
#ifdef LLMAP_HAS_FAISS
    std::unique_ptr<faiss::Index> index;
#ifdef LLMAP_HAS_FAISS_GPU
    std::unique_ptr<faiss::gpu::StandardGpuResources> gpu_resources;
#endif
#endif
    FaissIndexStats stats;
    std::string last_error;
};

FaissIndex::FaissIndex(const FaissIndexConfig& config)
    : config_(config)
    , active_provider_(FaissProvider::CPU)
    , impl_(std::make_unique<FaissIndexImpl>())
{
}

FaissIndex::~FaissIndex() = default;

FaissIndex::FaissIndex(FaissIndex&&) noexcept = default;
FaissIndex& FaissIndex::operator=(FaissIndex&&) noexcept = default;

std::unique_ptr<FaissIndex> FaissIndex::Create(const FaissIndexConfig& config) {
    auto index = std::unique_ptr<FaissIndex>(new FaissIndex(config));
    if (!index->Initialize()) {
        return nullptr;
    }
    return index;
}

bool FaissIndex::Initialize() {
#ifdef LLMAP_HAS_FAISS
    try {
        const auto dim = static_cast<int>(config_.embedding_dim);

        // Try GPU first if requested
        bool use_gpu = false;
#ifdef LLMAP_HAS_FAISS_GPU
        if (config_.provider == FaissProvider::GPU) {
            impl_->gpu_resources = std::make_unique<faiss::gpu::StandardGpuResources>();
            use_gpu = true;
            active_provider_ = FaissProvider::GPU;
        }
#else
        if (config_.provider == FaissProvider::GPU) {
            impl_->last_error = "FAISS-GPU requested but not compiled with GPU support";
            active_provider_ = FaissProvider::CPU;
        }
#endif

        // Create the index based on type
        switch (config_.index_type) {
            case FaissIndexType::FlatL2: {
#ifdef LLMAP_HAS_FAISS_GPU
                if (use_gpu) {
                    impl_->index = std::make_unique<faiss::gpu::GpuIndexFlatL2>(
                        impl_->gpu_resources.get(),
                        dim);
                } else
#endif
                {
                    impl_->index = std::make_unique<faiss::IndexFlatL2>(dim);
                }
                impl_->stats.is_trained = true;  // Flat indices are always trained
                break;
            }

            case FaissIndexType::FlatIP: {
#ifdef LLMAP_HAS_FAISS_GPU
                if (use_gpu) {
                    impl_->index = std::make_unique<faiss::gpu::GpuIndexFlatIP>(
                        impl_->gpu_resources.get(),
                        dim);
                } else
#endif
                {
                    impl_->index = std::make_unique<faiss::IndexFlatIP>(dim);
                }
                impl_->stats.is_trained = true;
                break;
            }

            case FaissIndexType::IVFFlat: {
                // IVF requires a quantizer (Flat index for cell centroids)
                auto quantizer = std::make_unique<faiss::IndexFlatL2>(dim);

#ifdef LLMAP_HAS_FAISS_GPU
                if (use_gpu) {
                    faiss::gpu::GpuIndexIVFFlatConfig gpu_config;
                    gpu_config.device = config_.gpu_device_id;
                    impl_->index = std::make_unique<faiss::gpu::GpuIndexIVFFlat>(
                        impl_->gpu_resources.get(),
                        quantizer.release(),  // GpuIndexIVFFlat takes ownership
                        dim,
                        static_cast<int>(config_.nlist),
                        faiss::METRIC_L2,
                        gpu_config);
                } else
#endif
                {
                    // CPU: IndexIVFFlat does NOT take ownership of quantizer
                    auto ivf = std::make_unique<faiss::IndexIVFFlat>(
                        quantizer.get(),
                        dim,
                        config_.nlist,
                        faiss::METRIC_L2);
                    ivf->own_fields = true;  // Take ownership
                    quantizer.release();
                    impl_->index = std::move(ivf);
                }
                impl_->stats.is_trained = false;  // IVF needs training
                break;
            }

            case FaissIndexType::IVFPQ: {
                auto quantizer = std::make_unique<faiss::IndexFlatL2>(dim);
                // CPU only for IVFPQ in this implementation
                if (use_gpu) {
                    impl_->last_error = "IVFPQ on GPU not yet supported; falling back to CPU";
                    active_provider_ = FaissProvider::CPU;
                }
                auto ivfpq = std::make_unique<faiss::IndexIVFPQ>(
                    quantizer.get(),
                    dim,
                    config_.nlist,
                    config_.m_subquantizers,
                    config_.nbits_per_idx);
                ivfpq->own_fields = true;
                quantizer.release();
                impl_->index = std::move(ivfpq);
                impl_->stats.is_trained = false;
                break;
            }
        }

        if (!impl_->index) {
            impl_->last_error = "Failed to create FAISS index";
            return false;
        }

        impl_->stats.embedding_dim = config_.embedding_dim;
        return true;

    } catch (const std::exception& e) {
        impl_->last_error = std::string("FAISS initialization failed: ") + e.what();
        return false;
    }
#else
    impl_->last_error = "LLmap compiled without FAISS support";
    return false;
#endif
}

bool FaissIndex::Train(std::span<const float> vectors, size_t num_vectors) {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) {
        impl_->last_error = "Index not initialized";
        return false;
    }

    if (vectors.size() < num_vectors * config_.embedding_dim) {
        impl_->last_error = "Vector buffer too small for specified num_vectors";
        return false;
    }

    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Normalize if configured
        std::vector<float> normalized;
        const float* data_ptr = vectors.data();
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(vectors, num_vectors);
            data_ptr = normalized.data();
        }

        // Subsample if training_sample_size is set
        if (config_.training_sample_size > 0 && config_.training_sample_size < num_vectors) {
            std::vector<float> sample(config_.training_sample_size * config_.embedding_dim);
            // Simple strided sampling
            const size_t stride = num_vectors / config_.training_sample_size;
            for (size_t i = 0; i < config_.training_sample_size; ++i) {
                const size_t src_idx = i * stride;
                std::memcpy(
                    sample.data() + i * config_.embedding_dim,
                    data_ptr + src_idx * config_.embedding_dim,
                    config_.embedding_dim * sizeof(float));
            }
            impl_->index->train(
                static_cast<faiss::idx_t>(config_.training_sample_size),
                sample.data());
        } else {
            impl_->index->train(
                static_cast<faiss::idx_t>(num_vectors),
                data_ptr);
        }

        auto end = std::chrono::high_resolution_clock::now();
        impl_->stats.training_time_ms = std::chrono::duration<float, std::milli>(end - start).count();
        impl_->stats.is_trained = true;
        return true;

    } catch (const std::exception& e) {
        impl_->last_error = std::string("Training failed: ") + e.what();
        return false;
    }
#else
    impl_->last_error = "LLmap compiled without FAISS support";
    return false;
#endif
}

size_t FaissIndex::Add(std::span<const float> vectors, size_t num_vectors) {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) {
        impl_->last_error = "Index not initialized";
        return 0;
    }

    if (!impl_->stats.is_trained) {
        impl_->last_error = "Index not trained; call Train() first for IVF indices";
        return 0;
    }

    if (vectors.size() < num_vectors * config_.embedding_dim) {
        impl_->last_error = "Vector buffer too small";
        return 0;
    }

    try {
        auto start = std::chrono::high_resolution_clock::now();

        const float* data_ptr = vectors.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(vectors, num_vectors);
            data_ptr = normalized.data();
        }

        impl_->index->add(static_cast<faiss::idx_t>(num_vectors), data_ptr);

        auto end = std::chrono::high_resolution_clock::now();
        impl_->stats.add_time_ms += std::chrono::duration<float, std::milli>(end - start).count();
        impl_->stats.num_vectors = impl_->index->ntotal;
        return num_vectors;

    } catch (const std::exception& e) {
        impl_->last_error = std::string("Add failed: ") + e.what();
        return 0;
    }
#else
    impl_->last_error = "LLmap compiled without FAISS support";
    return 0;
#endif
}

size_t FaissIndex::AddWithIds(
    std::span<const float> vectors,
    std::span<const int64_t> ids,
    size_t num_vectors)
{
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) {
        impl_->last_error = "Index not initialized";
        return 0;
    }

    if (!impl_->stats.is_trained) {
        impl_->last_error = "Index not trained";
        return 0;
    }

    if (ids.size() < num_vectors) {
        impl_->last_error = "IDs buffer too small";
        return 0;
    }

    try {
        auto start = std::chrono::high_resolution_clock::now();

        const float* data_ptr = vectors.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(vectors, num_vectors);
            data_ptr = normalized.data();
        }

        impl_->index->add_with_ids(
            static_cast<faiss::idx_t>(num_vectors),
            data_ptr,
            ids.data());

        auto end = std::chrono::high_resolution_clock::now();
        impl_->stats.add_time_ms += std::chrono::duration<float, std::milli>(end - start).count();
        impl_->stats.num_vectors = impl_->index->ntotal;
        return num_vectors;

    } catch (const std::exception& e) {
        impl_->last_error = std::string("AddWithIds failed: ") + e.what();
        return 0;
    }
#else
    impl_->last_error = "LLmap compiled without FAISS support";
    return 0;
#endif
}

bool FaissIndex::TrainAndAdd(std::span<const float> vectors, size_t num_vectors) {
    if (!Train(vectors, num_vectors)) {
        return false;
    }
    return Add(vectors, num_vectors) == num_vectors;
}

std::optional<KnnResult> FaissIndex::Search(std::span<const float> query, size_t k) const {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index || impl_->stats.num_vectors == 0) {
        return std::nullopt;
    }

    if (query.size() < config_.embedding_dim) {
        return std::nullopt;
    }

    try {
        KnnResult result;
        result.indices.resize(k);
        result.distances.resize(k);

        const float* data_ptr = query.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(query, 1);
            data_ptr = normalized.data();
        }

        impl_->index->search(
            1,
            data_ptr,
            static_cast<faiss::idx_t>(k),
            result.distances.data(),
            result.indices.data());

        return result;

    } catch (const std::exception&) {
        return std::nullopt;
    }
#else
    return std::nullopt;
#endif
}

std::optional<BatchKnnResult> FaissIndex::SearchBatch(
    std::span<const float> queries,
    size_t num_queries,
    size_t k) const
{
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index || impl_->stats.num_vectors == 0) {
        impl_->last_error = "Index empty or not initialized";
        return std::nullopt;
    }

    if (queries.size() < num_queries * config_.embedding_dim) {
        impl_->last_error = "Query buffer too small";
        return std::nullopt;
    }

    try {
        auto start = std::chrono::high_resolution_clock::now();

        BatchKnnResult result;
        result.num_queries = num_queries;
        result.k = k;
        result.indices.resize(num_queries * k);
        result.distances.resize(num_queries * k);

        const float* data_ptr = queries.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(queries, num_queries);
            data_ptr = normalized.data();
        }

        impl_->index->search(
            static_cast<faiss::idx_t>(num_queries),
            data_ptr,
            static_cast<faiss::idx_t>(k),
            result.distances.data(),
            result.indices.data());

        auto end = std::chrono::high_resolution_clock::now();
        result.search_time_ms = std::chrono::duration<float, std::milli>(end - start).count();

        return result;

    } catch (const std::exception& e) {
        impl_->last_error = std::string("SearchBatch failed: ") + e.what();
        return std::nullopt;
    }
#else
    return std::nullopt;
#endif
}

bool FaissIndex::SearchBatchInto(
    std::span<const float> queries,
    size_t num_queries,
    size_t k,
    std::span<int64_t> indices,
    std::span<float> distances) const
{
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index || impl_->stats.num_vectors == 0) {
        return false;
    }

    if (indices.size() < num_queries * k || distances.size() < num_queries * k) {
        return false;
    }

    try {
        const float* data_ptr = queries.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(queries, num_queries);
            data_ptr = normalized.data();
        }

        impl_->index->search(
            static_cast<faiss::idx_t>(num_queries),
            data_ptr,
            static_cast<faiss::idx_t>(k),
            distances.data(),
            indices.data());

        return true;

    } catch (const std::exception&) {
        return false;
    }
#else
    return false;
#endif
}

void FaissIndex::SetNprobe(size_t nprobe) {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) return;

    // nprobe is an IVF-specific parameter
    if (auto* ivf = dynamic_cast<faiss::IndexIVF*>(impl_->index.get())) {
        ivf->nprobe = nprobe;
    }
#ifdef LLMAP_HAS_FAISS_GPU
    if (auto* gpu_ivf = dynamic_cast<faiss::gpu::GpuIndexIVF*>(impl_->index.get())) {
        gpu_ivf->setNumProbes(static_cast<int>(nprobe));
    }
#endif
#endif
}

size_t FaissIndex::GetNprobe() const {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) return 0;

    if (auto* ivf = dynamic_cast<faiss::IndexIVF*>(impl_->index.get())) {
        return ivf->nprobe;
    }
#ifdef LLMAP_HAS_FAISS_GPU
    if (auto* gpu_ivf = dynamic_cast<faiss::gpu::GpuIndexIVF*>(impl_->index.get())) {
        return gpu_ivf->getNumProbes();
    }
#endif
#endif
    return config_.nprobe;
}

bool FaissIndex::IsReady() const {
#ifdef LLMAP_HAS_FAISS
    return impl_->index != nullptr && impl_->stats.is_trained;
#else
    return false;
#endif
}

bool FaissIndex::IsTrained() const {
    return impl_->stats.is_trained;
}

size_t FaissIndex::NumVectors() const {
    return impl_->stats.num_vectors;
}

FaissIndexStats FaissIndex::GetStats() const {
    return impl_->stats;
}

bool FaissIndex::IsGpuEnabled() const {
    return active_provider_ == FaissProvider::GPU;
}

std::string FaissIndex::LastError() const {
    return impl_->last_error;
}

bool FaissIndex::Save(const std::string& path) const {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) {
        return false;
    }
    try {
        faiss::write_index(impl_->index.get(), path.c_str());
        return true;
    } catch (const std::exception&) {
        return false;
    }
#else
    return false;
#endif
}

std::unique_ptr<FaissIndex> FaissIndex::Load(
    const std::string& path,
    const FaissIndexConfig& config)
{
#ifdef LLMAP_HAS_FAISS
    try {
        auto wrapper = std::unique_ptr<FaissIndex>(new FaissIndex(config));
        wrapper->impl_->index.reset(faiss::read_index(path.c_str()));

        if (!wrapper->impl_->index) {
            return nullptr;
        }

        wrapper->impl_->stats.num_vectors = wrapper->impl_->index->ntotal;
        wrapper->impl_->stats.embedding_dim = config.embedding_dim;
        wrapper->impl_->stats.is_trained = wrapper->impl_->index->is_trained;

        return wrapper;
    } catch (const std::exception&) {
        return nullptr;
    }
#else
    return nullptr;
#endif
}

void FaissIndex::NormalizeVectors(std::span<float> vectors, size_t num_vectors) const {
    const size_t dim = config_.embedding_dim;
    for (size_t i = 0; i < num_vectors; ++i) {
        float* vec = vectors.data() + i * dim;
        float norm = 0.0f;
        for (size_t j = 0; j < dim; ++j) {
            norm += vec[j] * vec[j];
        }
        norm = std::sqrt(norm);
        if (norm > 1e-8f) {
            const float inv_norm = 1.0f / norm;
            for (size_t j = 0; j < dim; ++j) {
                vec[j] *= inv_norm;
            }
        }
    }
}

std::vector<float> FaissIndex::NormalizeVectorsCopy(
    std::span<const float> vectors,
    size_t num_vectors) const
{
    std::vector<float> normalized(vectors.begin(),
        vectors.begin() + num_vectors * config_.embedding_dim);
    NormalizeVectors(normalized, num_vectors);
    return normalized;
}

// ========== Global utilities ==========

bool IsFaissAvailable() {
#ifdef LLMAP_HAS_FAISS
    return true;
#else
    return false;
#endif
}

bool IsFaissGpuAvailable() {
#ifdef LLMAP_HAS_FAISS_GPU
    return true;
#else
    return false;
#endif
}

std::string GetFaissVersion() {
#ifdef LLMAP_HAS_FAISS
    // FAISS doesn't expose a version string easily; return compile-time info
    return "FAISS (compiled)";
#else
    return "FAISS (not available)";
#endif
}

float ComputeRecallAtK(
    std::span<const int64_t> approx_indices,
    std::span<const int64_t> exact_indices,
    size_t num_queries,
    size_t k)
{
    if (num_queries == 0 || k == 0) {
        return 0.0f;
    }

    if (approx_indices.size() < num_queries * k ||
        exact_indices.size() < num_queries * k) {
        return 0.0f;
    }

    size_t total_found = 0;
    for (size_t q = 0; q < num_queries; ++q) {
        const auto* approx = approx_indices.data() + q * k;
        const auto* exact = exact_indices.data() + q * k;

        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                if (approx[i] == exact[j] && approx[i] >= 0) {
                    ++total_found;
                    break;
                }
            }
        }
    }

    return static_cast<float>(total_found) / static_cast<float>(num_queries * k);
}

}  // namespace llmap::self_interference
