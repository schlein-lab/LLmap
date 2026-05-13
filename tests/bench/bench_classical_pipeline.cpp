// Classical pipeline benchmark
// Run with: ./bench_classical_pipeline [num_reads] [read_length] [ref_length]
//
// Benchmarks the seed-chain-extend alignment path:
//   - Minimizer extraction
//   - Index lookup and chaining
//   - WFA2 extension alignment

#include "classical/classical_pipeline.h"
#include "classical/minimizer_index.h"
#include "classical/chain.h"
#include "classical/wfa2_aligner.h"
#include "core/profiler.h"

#include <chrono>
#include <cstdlib>
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

std::string MutateSequence(const std::string& seq, double error_rate, std::mt19937& gen) {
    static const char nucleotides[] = "ACGT";
    std::uniform_real_distribution<> prob(0.0, 1.0);
    std::uniform_int_distribution<> base(0, 3);
    std::uniform_int_distribution<> error_type(0, 2);  // sub, ins, del

    std::string result;
    result.reserve(seq.size());

    for (char c : seq) {
        if (prob(gen) < error_rate) {
            int type = error_type(gen);
            if (type == 0) {
                result += nucleotides[base(gen)];
            } else if (type == 1) {
                result += c;
                result += nucleotides[base(gen)];
            }
            // type == 2: deletion, skip the character
        } else {
            result += c;
        }
    }
    return result;
}

void PrintSeparator() {
    std::cout << std::string(70, '-') << "\n";
}

void BenchmarkMinimizerExtraction(
    const std::vector<std::string>& reads,
    const llmap::classical::MinimizerConfig& config) {
    std::cout << "\n=== Minimizer Extraction Benchmark ===\n";

    llmap::core::ManualTimer timer;
    size_t total_minimizers = 0;

    timer.Start();
    for (const auto& read : reads) {
        auto minimizers = llmap::classical::ExtractMinimizers(read, config);
        total_minimizers += minimizers.size();
    }
    timer.Stop();

    double reads_per_sec = reads.size() / (timer.ElapsedMs() / 1000.0);
    double us_per_read = timer.ElapsedUs() / reads.size();
    double minimizers_per_read = static_cast<double>(total_minimizers) / reads.size();

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total reads:        " << reads.size() << "\n"
              << "Total minimizers:   " << total_minimizers << "\n"
              << "Minimizers/read:    " << minimizers_per_read << "\n"
              << "Total time:         " << timer.ElapsedMs() << " ms\n"
              << "Per-read time:      " << us_per_read << " µs\n"
              << "Throughput:         " << reads_per_sec << " reads/sec\n";
}

void BenchmarkChaining(
    const std::vector<std::vector<llmap::classical::Anchor>>& anchor_sets,
    const llmap::classical::ChainConfig& config) {
    std::cout << "\n=== Chaining Benchmark ===\n";

    llmap::core::ManualTimer timer;
    size_t total_chains = 0;

    timer.Start();
    for (const auto& anchors : anchor_sets) {
        if (!anchors.empty()) {
            auto result = llmap::classical::ExtractChainsFromAnchors(anchors, 150, config);
            total_chains += result.chains.size();
        }
    }
    timer.Stop();

    double sets_per_sec = anchor_sets.size() / (timer.ElapsedMs() / 1000.0);
    double us_per_set = timer.ElapsedUs() / anchor_sets.size();
    double chains_per_set = static_cast<double>(total_chains) / anchor_sets.size();

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total anchor sets:  " << anchor_sets.size() << "\n"
              << "Total chains:       " << total_chains << "\n"
              << "Chains/set:         " << chains_per_set << "\n"
              << "Total time:         " << timer.ElapsedMs() << " ms\n"
              << "Per-set time:       " << us_per_set << " µs\n"
              << "Throughput:         " << sets_per_sec << " sets/sec\n";
}

void BenchmarkWFA2Extension(
    const std::vector<std::pair<std::string, std::string>>& pairs,
    const llmap::classical::WFA2Config& config) {
    std::cout << "\n=== WFA2 Extension Benchmark ===\n";

    llmap::classical::WFA2Aligner aligner(config);
    llmap::core::ManualTimer timer;
    size_t successful = 0;

    timer.Start();
    for (const auto& [query, ref] : pairs) {
        auto result = aligner.Align(query, ref);
        if (result.has_value()) {
            successful++;
        }
    }
    timer.Stop();

    double pairs_per_sec = pairs.size() / (timer.ElapsedMs() / 1000.0);
    double us_per_pair = timer.ElapsedUs() / pairs.size();
    double success_rate = 100.0 * successful / pairs.size();

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total pairs:        " << pairs.size() << "\n"
              << "Successful:         " << successful << " (" << success_rate << "%)\n"
              << "Total time:         " << timer.ElapsedMs() << " ms\n"
              << "Per-pair time:      " << us_per_pair << " µs\n"
              << "Throughput:         " << pairs_per_sec << " pairs/sec\n";
}

void BenchmarkFullPipeline(
    const std::vector<std::string>& reads,
    const std::vector<std::string>& read_names,
    llmap::classical::ClassicalPipeline& pipeline) {
    std::cout << "\n=== Full Pipeline Benchmark ===\n";

    llmap::core::ManualTimer timer;

    timer.Start();
    auto results = pipeline.AlignReads(read_names, reads);
    timer.Stop();

    const auto& stats = pipeline.Stats();

    double reads_per_sec = reads.size() / (timer.ElapsedMs() / 1000.0);
    double us_per_read = timer.ElapsedUs() / reads.size();
    double mapping_rate = 100.0 * stats.reads_aligned / stats.total_reads;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total reads:        " << stats.total_reads << "\n"
              << "Aligned:            " << stats.reads_aligned << " (" << mapping_rate << "%)\n"
              << "Unmapped:           " << stats.reads_unmapped << "\n"
              << "Total time:         " << timer.ElapsedMs() << " ms\n"
              << "Per-read time:      " << us_per_read << " µs\n"
              << "Throughput:         " << reads_per_sec << " reads/sec\n"
              << "\nPhase breakdown:\n"
              << "  Seeding:          " << stats.seeding_time_ms << " ms\n"
              << "  Chaining:         " << stats.chaining_time_ms << " ms\n"
              << "  Extension:        " << stats.extension_time_ms << " ms\n"
              << "\nAverage identity:   " << (stats.avg_identity * 100) << "%\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    size_t num_reads = (argc > 1) ? std::atoi(argv[1]) : 10000;
    size_t read_length = (argc > 2) ? std::atoi(argv[2]) : 150;
    size_t ref_length = (argc > 3) ? std::atoi(argv[3]) : 100000;

    std::cout << "\n";
    PrintSeparator();
    std::cout << "LLmap Classical Pipeline Benchmark\n";
    PrintSeparator();
    std::cout << "Reads:             " << num_reads << "\n"
              << "Read length:       " << read_length << " bp\n"
              << "Reference length:  " << ref_length << " bp\n";
    PrintSeparator();
    std::cout << "\n";

    // Generate reference and reads
    std::cout << "Generating synthetic data...\n";
    std::mt19937 gen(42);

    std::string ref_seq = GenerateRandomSequence(ref_length, gen);
    std::vector<std::string> ref_names = {"chr1"};
    std::vector<std::string> ref_seqs = {ref_seq};

    std::uniform_int_distribution<size_t> pos_dist(0, ref_length - read_length - 1);
    std::vector<std::string> reads;
    std::vector<std::string> read_names;
    reads.reserve(num_reads);
    read_names.reserve(num_reads);

    for (size_t i = 0; i < num_reads; ++i) {
        size_t start = pos_dist(gen);
        std::string read = ref_seq.substr(start, read_length);
        read = MutateSequence(read, 0.01, gen);  // 1% error rate
        reads.push_back(read);
        read_names.push_back("read_" + std::to_string(i));
    }
    std::cout << "Done.\n";

    // Configuration
    llmap::classical::MinimizerConfig mini_config;
    mini_config.k = 15;
    mini_config.w = 10;

    llmap::classical::ChainConfig chain_config;
    chain_config.max_gap_ref = 500;
    chain_config.min_chain_score = 30;

    llmap::classical::WFA2Config wfa_config;
    wfa_config.match_score = 0;
    wfa_config.mismatch_penalty = 4;
    wfa_config.gap_open = 6;
    wfa_config.gap_extend = 2;
    wfa_config.max_alignment_score = 1000;

    llmap::classical::ClassicalPipelineConfig pipeline_config;
    pipeline_config.minimizer_config = mini_config;
    pipeline_config.chain_config = chain_config;
    pipeline_config.extension_config = wfa_config;

    // Run individual benchmarks
    BenchmarkMinimizerExtraction(reads, mini_config);

    // Build index for chaining benchmark
    std::cout << "\nBuilding minimizer index...\n";
    llmap::classical::MinimizerIndex::Builder builder(mini_config);
    for (size_t i = 0; i < ref_seqs.size(); ++i) {
        builder.AddSequence(ref_names[i], ref_seqs[i]);
    }
    auto index = builder.Build();
    std::cout << "Index built with " << index->Size() << " minimizers.\n";

    // Generate anchor sets for chaining benchmark
    std::cout << "\nGenerating anchor sets for chaining...\n";
    std::vector<std::vector<llmap::classical::Anchor>> anchor_sets;
    anchor_sets.reserve(num_reads);
    for (const auto& read : reads) {
        auto hits = index->Query(read);
        std::vector<llmap::classical::Anchor> anchors;
        for (const auto& hit : hits) {
            anchors.push_back({
                hit.ref_id,
                hit.ref_pos,
                hit.query_pos,
                hit.same_strand
            });
        }
        anchor_sets.push_back(std::move(anchors));
    }

    BenchmarkChaining(anchor_sets, chain_config);

    // Generate pairs for WFA2 benchmark (subsample)
    std::cout << "\nPreparing sequence pairs for WFA2...\n";
    size_t wfa_sample = std::min(num_reads, size_t(1000));
    std::vector<std::pair<std::string, std::string>> pairs;
    pairs.reserve(wfa_sample);
    for (size_t i = 0; i < wfa_sample; ++i) {
        size_t start = pos_dist(gen);
        std::string ref_region = ref_seq.substr(start, read_length + 20);
        pairs.emplace_back(reads[i], ref_region);
    }

    BenchmarkWFA2Extension(pairs, wfa_config);

    // Full pipeline benchmark
    std::cout << "\nPreparing full pipeline...\n";
    llmap::classical::ClassicalPipeline pipeline(pipeline_config);
    pipeline.SetIndex(std::move(index));

    BenchmarkFullPipeline(reads, read_names, pipeline);

    std::cout << "\n";
    PrintSeparator();
    std::cout << "Benchmark complete.\n";
    PrintSeparator();

    return 0;
}
