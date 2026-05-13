// Embedder throughput benchmark
// Run with: ./bench_embedder_throughput [model_path] [num_sequences] [sequence_length]
//
// This is a simple standalone benchmark, not using Google Benchmark framework
// to keep dependencies minimal for Phase 1.

#include "ai/foundation_embedder.h"
#include "ai/bucket_embedder.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

std::string GenerateRandomSequence(size_t length, std::mt19937& gen) {
    static const char nucleotides[] = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);

    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq += nucleotides[dist(gen)];
    }
    return seq;
}

void PrintUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " <model_path> [num_sequences] [sequence_length]\n"
              << "  model_path:      Path to ONNX model file\n"
              << "  num_sequences:   Number of sequences to embed (default: 1000)\n"
              << "  sequence_length: Length of each sequence (default: 150)\n";
}

void PrintSeparator() {
    std::cout << std::string(70, '-') << "\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }

    std::filesystem::path model_path = argv[1];
    size_t num_sequences = (argc > 2) ? std::atoi(argv[2]) : 1000;
    size_t seq_length = (argc > 3) ? std::atoi(argv[3]) : 150;

    if (!std::filesystem::exists(model_path)) {
        std::cerr << "Error: Model file not found: " << model_path << "\n";
        return 1;
    }

    std::cout << "\n";
    PrintSeparator();
    std::cout << "LLmap Embedder Throughput Benchmark\n";
    PrintSeparator();
    std::cout << "Model:            " << model_path.filename().string() << "\n"
              << "Sequences:        " << num_sequences << "\n"
              << "Sequence length:  " << seq_length << " bp\n"
              << "Total nucleotides:" << (num_sequences * seq_length) << " bp\n";
    PrintSeparator();
    std::cout << "\n";

    // Check ONNX Runtime availability
    if (!llmap::ai::IsOnnxRuntimeAvailable()) {
        std::cerr << "Error: ONNX Runtime not available\n";
        return 1;
    }

    std::cout << "ONNX Runtime: Available\n";
    std::cout << "Providers: ";
    for (const auto& p : llmap::ai::ListAvailableProviders()) {
        std::cout << p << " ";
    }
    std::cout << "\n\n";

    // Generate test sequences
    std::cout << "Generating " << num_sequences << " random sequences...\n";
    std::mt19937 gen(42);
    std::vector<std::string> sequences_storage;
    sequences_storage.reserve(num_sequences);
    for (size_t i = 0; i < num_sequences; ++i) {
        sequences_storage.push_back(GenerateRandomSequence(seq_length, gen));
    }

    std::vector<std::string_view> sequences;
    sequences.reserve(num_sequences);
    for (const auto& s : sequences_storage) {
        sequences.push_back(s);
    }
    std::cout << "Done.\n\n";

    // FoundationEmbedder benchmark
    {
        std::cout << "=== FoundationEmbedder (read embedder) ===\n";

        llmap::ai::EmbedderConfig config;
        config.model_path = model_path;
        config.provider = llmap::ai::ExecutionProvider::CPU;
        config.embedding_dim = 256;
        config.max_sequence_length = 512;
        config.batch_size = 64;

        auto embedder = llmap::ai::FoundationEmbedder::Create(config);
        if (!embedder) {
            std::cerr << "Failed to create FoundationEmbedder\n";
            return 1;
        }

        std::cout << "Active provider: "
                  << (embedder->ActiveProvider() == llmap::ai::ExecutionProvider::CUDA ? "CUDA" :
                      embedder->ActiveProvider() == llmap::ai::ExecutionProvider::TensorRT ? "TensorRT" : "CPU")
                  << "\n";

        // Warm up
        std::cout << "Warming up...\n";
        for (int i = 0; i < 10; ++i) {
            embedder->Embed(sequences[0]);
        }

        // Single sequence benchmark
        std::cout << "\n--- Single sequence embedding ---\n";
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < num_sequences; ++i) {
            auto result = embedder->Embed(sequences[i]);
            if (!result) {
                std::cerr << "Embed failed at sequence " << i << "\n";
                return 1;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();

        double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        double reads_per_sec = num_sequences / (elapsed_ms / 1000.0);
        double us_per_read = (elapsed_ms * 1000.0) / num_sequences;
        double mb_per_sec = (num_sequences * seq_length) / (elapsed_ms / 1000.0) / 1e6;

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Total time:       " << elapsed_ms << " ms\n"
                  << "Per-read time:    " << us_per_read << " µs\n"
                  << "Throughput:       " << reads_per_sec << " reads/sec\n"
                  << "                  " << mb_per_sec << " Mb/sec\n";

        // Batch benchmark
        std::cout << "\n--- Batch embedding (batch_size=" << config.batch_size << ") ---\n";
        start = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < num_sequences; i += config.batch_size) {
            size_t batch_end = std::min(i + config.batch_size, num_sequences);
            std::span<const std::string_view> batch(sequences.data() + i, batch_end - i);
            auto result = embedder->EmbedBatch(batch);
            if (!result) {
                std::cerr << "EmbedBatch failed at batch starting " << i << "\n";
                return 1;
            }
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        reads_per_sec = num_sequences / (elapsed_ms / 1000.0);
        us_per_read = (elapsed_ms * 1000.0) / num_sequences;
        mb_per_sec = (num_sequences * seq_length) / (elapsed_ms / 1000.0) / 1e6;

        std::cout << "Total time:       " << elapsed_ms << " ms\n"
                  << "Per-read time:    " << us_per_read << " µs\n"
                  << "Throughput:       " << reads_per_sec << " reads/sec\n"
                  << "                  " << mb_per_sec << " Mb/sec\n";
    }

    std::cout << "\n";

    // BucketEmbedder benchmark (with longer sequences)
    {
        std::cout << "=== BucketEmbedder (reference bucket embedder) ===\n";

        // Generate longer bucket-like sequences
        size_t bucket_len = 1000;  // 1kb buckets
        size_t num_buckets = num_sequences / 10;  // fewer but larger

        std::vector<std::string> bucket_storage;
        bucket_storage.reserve(num_buckets);
        std::mt19937 gen2(123);
        for (size_t i = 0; i < num_buckets; ++i) {
            bucket_storage.push_back(GenerateRandomSequence(bucket_len, gen2));
        }

        std::vector<std::string_view> buckets;
        buckets.reserve(num_buckets);
        for (const auto& s : bucket_storage) {
            buckets.push_back(s);
        }

        std::cout << "Bucket count:     " << num_buckets << "\n"
                  << "Bucket length:    " << bucket_len << " bp\n"
                  << "Total bp:         " << (num_buckets * bucket_len) << "\n\n";

        llmap::ai::BucketEmbedderConfig config;
        config.model_path = model_path;
        config.provider = llmap::ai::ExecutionProvider::CPU;
        config.embedding_dim = 256;
        config.max_sequence_length = 2048;
        config.batch_size = 32;

        auto embedder = llmap::ai::BucketEmbedder::Create(config);
        if (!embedder) {
            std::cerr << "Failed to create BucketEmbedder\n";
            return 1;
        }

        std::cout << "Active provider: "
                  << (embedder->ActiveProvider() == llmap::ai::ExecutionProvider::CUDA ? "CUDA" :
                      embedder->ActiveProvider() == llmap::ai::ExecutionProvider::TensorRT ? "TensorRT" : "CPU")
                  << "\n";

        // Warm up
        std::cout << "Warming up...\n";
        for (int i = 0; i < 5; ++i) {
            embedder->EmbedBucket(buckets[0]);
        }

        // Batch benchmark with stats
        std::cout << "\n--- Batch embedding with stats ---\n";

        std::vector<float> output(num_buckets * config.embedding_dim);
        llmap::ai::BucketEmbeddingStats stats;

        auto start = std::chrono::high_resolution_clock::now();
        bool success = embedder->EmbedBucketsInto(buckets, output, &stats);
        auto end = std::chrono::high_resolution_clock::now();

        if (!success) {
            std::cerr << "EmbedBucketsInto failed\n";
            return 1;
        }

        double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Total time:       " << elapsed_ms << " ms\n"
                  << "Per-bucket time:  " << stats.avg_time_per_bucket_ms << " ms\n"
                  << "Throughput:       " << stats.throughput_buckets_per_sec << " buckets/sec\n"
                  << "                  " << stats.throughput_mb_per_sec << " Mb/sec\n";
    }

    std::cout << "\n";
    PrintSeparator();
    std::cout << "Benchmark complete.\n";
    PrintSeparator();

    return 0;
}
