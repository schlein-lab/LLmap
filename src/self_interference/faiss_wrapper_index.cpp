#include "self_interference/faiss_wrapper_impl.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>

#ifdef LLMAP_HAS_FAISS
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/IndexIVFPQ.h>
#include <faiss/index_io.h>
#ifdef LLMAP_HAS_FAISS_GPU
#include <faiss/gpu/GpuIndexFlat.h>
#include <faiss/gpu/GpuIndexIVFFlat.h>
#include <faiss/gpu/GpuResources.h>
#endif
#endif

namespace llmap::self_interference {

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

        switch (config_.index_type) {
            case FaissIndexType::FlatL2:
                impl_->index = CreateFlatL2Index(dim, use_gpu);
                impl_->stats.is_trained = true;
                break;

            case FaissIndexType::FlatIP:
                impl_->index = CreateFlatIPIndex(dim, use_gpu);
                impl_->stats.is_trained = true;
                break;

            case FaissIndexType::IVFFlat:
                impl_->index = CreateIVFFlatIndex(dim, use_gpu);
                impl_->stats.is_trained = false;
                break;

            case FaissIndexType::IVFPQ:
                impl_->index = CreateIVFPQIndex(dim, use_gpu);
                impl_->stats.is_trained = false;
                break;
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

#ifdef LLMAP_HAS_FAISS
std::unique_ptr<faiss::Index> FaissIndex::CreateFlatL2Index(int dim, bool use_gpu) {
#ifdef LLMAP_HAS_FAISS_GPU
    if (use_gpu) {
        return std::make_unique<faiss::gpu::GpuIndexFlatL2>(
            impl_->gpu_resources.get(), dim);
    }
#endif
    (void)use_gpu;
    return std::make_unique<faiss::IndexFlatL2>(dim);
}

std::unique_ptr<faiss::Index> FaissIndex::CreateFlatIPIndex(int dim, bool use_gpu) {
#ifdef LLMAP_HAS_FAISS_GPU
    if (use_gpu) {
        return std::make_unique<faiss::gpu::GpuIndexFlatIP>(
            impl_->gpu_resources.get(), dim);
    }
#endif
    (void)use_gpu;
    return std::make_unique<faiss::IndexFlatIP>(dim);
}

std::unique_ptr<faiss::Index> FaissIndex::CreateIVFFlatIndex(int dim, bool use_gpu) {
    auto quantizer = std::make_unique<faiss::IndexFlatL2>(dim);

#ifdef LLMAP_HAS_FAISS_GPU
    if (use_gpu) {
        faiss::gpu::GpuIndexIVFFlatConfig gpu_config;
        gpu_config.device = config_.gpu_device_id;
        return std::make_unique<faiss::gpu::GpuIndexIVFFlat>(
            impl_->gpu_resources.get(),
            quantizer.release(),
            dim,
            static_cast<int>(config_.nlist),
            faiss::METRIC_L2,
            gpu_config);
    }
#endif
    (void)use_gpu;
    auto ivf = std::make_unique<faiss::IndexIVFFlat>(
        quantizer.get(), dim, config_.nlist, faiss::METRIC_L2);
    ivf->own_fields = true;
    quantizer.release();
    return ivf;
}

std::unique_ptr<faiss::Index> FaissIndex::CreateIVFPQIndex(int dim, bool use_gpu) {
    if (use_gpu) {
        impl_->last_error = "IVFPQ on GPU not yet supported; falling back to CPU";
        active_provider_ = FaissProvider::CPU;
    }
    auto quantizer = std::make_unique<faiss::IndexFlatL2>(dim);
    auto ivfpq = std::make_unique<faiss::IndexIVFPQ>(
        quantizer.get(), dim, config_.nlist,
        config_.m_subquantizers, config_.nbits_per_idx);
    ivfpq->own_fields = true;
    quantizer.release();
    return ivfpq;
}
#endif

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

        std::vector<float> normalized;
        const float* data_ptr = vectors.data();
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(vectors, num_vectors);
            data_ptr = normalized.data();
        }

        if (config_.training_sample_size > 0 && config_.training_sample_size < num_vectors) {
            std::vector<float> sample(config_.training_sample_size * config_.embedding_dim);
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

}  // namespace llmap::self_interference
