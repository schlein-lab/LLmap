// LLmap — triobin: parental-haplotype binning of child reads.
//
// Builds parent-specific k-mer membership filters from each parent's sequences
// (reads or assembly, plain/gz FASTA/FASTQ) and classifies every child read as
// paternal / maternal / ambiguous by which parent's specific k-mers it carries.
// This is the standard trio-binning truth signal (cf. yak triobin), implemented
// as a first-class LLmap subcommand on LLmap's io + threading primitives.
//
// Memory: whole-genome parental k-mer sets are billions of k-mers, so we use a
// pair of Bloom filters (fixed memory, tunable false-positive rate) rather than
// exact hash sets. Parent-specific = present in one parent's filter and absent
// in the other's; a child read is binned by majority over its specific-k-mer
// hits. Bloom false positives add a little noise that the per-read majority vote
// absorbs.
//
// Streaming: parents/child may be given as a path, "-"/"/dev/stdin", or a
// process-substitution path (e.g. <(curl s3-url)), so parental reads can be
// streamed without staging to disk.

#include "cli/commands.h"
#include "io/fastq_reader.h"

#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::cli {
namespace {

// ---- compact Bloom filter over 64-bit k-mer hashes ------------------------
class BloomFilter {
public:
    explicit BloomFilter(std::uint64_t n_bits, int n_hash)
        : bits_(n_bits / 64 + 1, 0), m_(n_bits), k_(n_hash) {}

    void add(std::uint64_t h) {
        for (int i = 0; i < k_; ++i) {
            std::uint64_t p = mix(h, i) % m_;
            bits_[p >> 6] |= (1ULL << (p & 63));
        }
    }
    bool maybe(std::uint64_t h) const {
        for (int i = 0; i < k_; ++i) {
            std::uint64_t p = mix(h, i) % m_;
            if (!(bits_[p >> 6] & (1ULL << (p & 63)))) return false;
        }
        return true;
    }

private:
    static std::uint64_t mix(std::uint64_t h, int i) {
        // splitmix64 seeded by the hash function index.
        std::uint64_t z = h + 0x9E3779B97F4A7C15ULL * (std::uint64_t)(i + 1);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    std::vector<std::uint64_t> bits_;
    std::uint64_t m_;
    int k_;
};

inline int base2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;  // N / other → breaks the k-mer
    }
}

inline std::uint64_t hash64(std::uint64_t key, std::uint64_t mask) {
    // minimap2-style invertible integer hash.
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

// Visit every canonical k-mer hash in `seq`, calling fn(hash).
template <class F>
void for_each_canonical_kmer(std::string_view seq, int k, F&& fn) {
    const std::uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    const int shift = 2 * (k - 1);
    std::uint64_t fwd = 0, rev = 0;
    int len = 0;
    for (char c : seq) {
        int b = base2bit(c);
        if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
        fwd = ((fwd << 2) | (std::uint64_t)b) & mask;
        rev = (rev >> 2) | ((std::uint64_t)(3 - b) << shift);
        if (++len >= k) {
            std::uint64_t canon = fwd < rev ? fwd : rev;
            fn(hash64(canon, mask));
        }
    }
}

struct Args {
    std::string paternal, maternal, child, output;
    int k = 31;
    std::uint64_t bloom_gib = 4;     // per-parent Bloom size, GiB
    double margin = 2.0;             // require winner >= margin × loser
    int min_specific = 5;            // min specific k-mers to call a parent
    bool help = false;
};

void usage() {
    std::puts(
        "Usage: llmap triobin --paternal P --maternal M --child C -o OUT [opts]\n\n"
        "Bin child reads to parental haplotypes via parent-specific k-mers.\n"
        "Inputs are FASTA/FASTQ (plain or .gz); use '-'/'/dev/stdin' or a\n"
        "process-substitution path to stream (e.g. --paternal <(curl s3-url)).\n\n"
        "Required:\n"
        "  --paternal FILE   paternal reads/assembly\n"
        "  --maternal FILE   maternal reads/assembly\n"
        "  --child FILE      child reads to classify\n"
        "  -o, --output FILE TSV: read_id<TAB>bin<TAB>pat_hits<TAB>mat_hits\n\n"
        "Options:\n"
        "  -k INT            k-mer size [31]\n"
        "  --bloom-gib INT   per-parent Bloom filter size in GiB [4]\n"
        "  --margin FLOAT    winner must beat loser by this factor [2.0]\n"
        "  --min-specific N  min specific k-mers to assign a parent [5]\n");
}

bool parse(int argc, char** argv, Args& a) {
    for (int i = 0; i < argc; ++i) {
        std::string s = argv[i];
        auto val = [&](const char* name) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", name); return nullptr; }
            return argv[++i];
        };
        if (s == "--paternal") { auto v = val("--paternal"); if (!v) return false; a.paternal = v; }
        else if (s == "--maternal") { auto v = val("--maternal"); if (!v) return false; a.maternal = v; }
        else if (s == "--child") { auto v = val("--child"); if (!v) return false; a.child = v; }
        else if (s == "-o" || s == "--output") { auto v = val("--output"); if (!v) return false; a.output = v; }
        else if (s == "-k") { auto v = val("-k"); if (!v) return false; a.k = std::atoi(v); }
        else if (s == "--bloom-gib") { auto v = val("--bloom-gib"); if (!v) return false; a.bloom_gib = std::strtoull(v, nullptr, 10); }
        else if (s == "--margin") { auto v = val("--margin"); if (!v) return false; a.margin = std::atof(v); }
        else if (s == "--min-specific") { auto v = val("--min-specific"); if (!v) return false; a.min_specific = std::atoi(v); }
        else if (s == "-h" || s == "--help") { a.help = true; return true; }
        else { std::fprintf(stderr, "unknown arg: %s\n", s.c_str()); return false; }
    }
    return true;
}

// Add every read's k-mers (FASTQ, plain or .gz, streamable) into a Bloom
// filter. Returns false on open failure.
bool fill_bloom(const std::string& path, int k, BloomFilter& bf,
                std::uint64_t& n_kmers) {
    n_kmers = 0;
    auto fq = io::FastqReader::Open(path);
    if (!fq) return false;
    while (fq->HasMore()) {
        auto recs = fq->NextBatch(4096);
        if (recs.empty()) break;
        for (const auto& r : recs)
            for_each_canonical_kmer(r.sequence, k,
                [&](std::uint64_t h) { bf.add(h); ++n_kmers; });
    }
    return true;
}

}  // namespace

int run_triobin(int argc, char** argv) {
    Args a;
    if (!parse(argc, argv, a)) { usage(); return 1; }
    if (a.help) { usage(); return 0; }
    if (a.paternal.empty() || a.maternal.empty() || a.child.empty() || a.output.empty()) {
        usage(); return 1;
    }

    const std::uint64_t bits = a.bloom_gib * 8ULL * 1024 * 1024 * 1024;
    BloomFilter pat(bits, 4), mat(bits, 4);

    std::fprintf(stderr, "[triobin] building paternal k-mer filter (k=%d)\n", a.k);
    std::uint64_t np = 0, nm = 0;
    if (!fill_bloom(a.paternal, a.k, pat, np)) {
        std::fprintf(stderr, "[triobin] cannot read paternal: %s\n", a.paternal.c_str());
        return 1;
    }
    std::fprintf(stderr, "[triobin] paternal k-mers added: %llu\n", (unsigned long long)np);

    std::fprintf(stderr, "[triobin] building maternal k-mer filter\n");
    if (!fill_bloom(a.maternal, a.k, mat, nm)) {
        std::fprintf(stderr, "[triobin] cannot read maternal: %s\n", a.maternal.c_str());
        return 1;
    }
    std::fprintf(stderr, "[triobin] maternal k-mers added: %llu\n", (unsigned long long)nm);

    std::FILE* out = std::fopen(a.output.c_str(), "w");
    if (!out) { std::fprintf(stderr, "[triobin] cannot write %s\n", a.output.c_str()); return 1; }
    std::fprintf(out, "read_id\tbin\tpat_hits\tmat_hits\n");

    auto fq = io::FastqReader::Open(a.child);
    if (!fq) { std::fprintf(stderr, "[triobin] cannot read child: %s\n", a.child.c_str()); std::fclose(out); return 1; }

    std::uint64_t n_pat = 0, n_mat = 0, n_amb = 0, n_tot = 0;
    while (fq->HasMore()) {
        auto recs = fq->NextBatch(8192);
        if (recs.empty()) break;
        for (const auto& r : recs) {
            int p = 0, m = 0;
            for_each_canonical_kmer(r.sequence, a.k, [&](std::uint64_t h) {
                const bool in_p = pat.maybe(h);
                const bool in_m = mat.maybe(h);
                if (in_p && !in_m) ++p;        // paternal-specific
                else if (in_m && !in_p) ++m;   // maternal-specific
            });
            const char* bin;
            if (p >= a.min_specific && p >= a.margin * m) { bin = "paternal"; ++n_pat; }
            else if (m >= a.min_specific && m >= a.margin * p) { bin = "maternal"; ++n_mat; }
            else { bin = "ambiguous"; ++n_amb; }
            std::fprintf(out, "%s\t%s\t%d\t%d\n", r.id.c_str(), bin, p, m);
            ++n_tot;
        }
    }
    std::fclose(out);
    std::fprintf(stderr,
        "[triobin] done: %llu reads | paternal=%llu maternal=%llu ambiguous=%llu\n",
        (unsigned long long)n_tot, (unsigned long long)n_pat,
        (unsigned long long)n_mat, (unsigned long long)n_amb);
    return 0;
}

}  // namespace llmap::cli
