// LLmap — Input sniffer implementation.

#include "io/input_sniffer.h"

#include "core/logging.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::io {

namespace {

// Lowercase ASCII copy.
std::string ToLower(std::string_view s) {
    std::string out(s);
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

// Basename (strip directory).
std::string Basename(const std::string& path) {
    const auto slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

bool ContainsAny(std::string_view hay, std::initializer_list<std::string_view> needles) {
    for (const auto n : needles) {
        if (hay.find(n) != std::string_view::npos) return true;
    }
    return false;
}

bool EndsWith(std::string_view s, std::string_view suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Median of a length vector (mutates by sorting). 0 if empty.
std::uint64_t Median(std::vector<std::uint64_t>& v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const std::size_t mid = v.size() / 2;
    if (v.size() % 2 == 1) return v[mid];
    return (v[mid - 1] + v[mid]) / 2;
}

// N50 over a length vector: smallest length L such that the cumulative sum
// of all lengths >= L (taken from the largest down) reaches half the total.
std::uint64_t N50(std::vector<std::uint64_t> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end(), std::greater<>{});
    std::uint64_t total = 0;
    for (const auto x : v) total += x;
    const std::uint64_t half = total / 2;
    std::uint64_t acc = 0;
    for (const auto x : v) {
        acc += x;
        if (acc * 2 >= total || acc >= half) return x;
    }
    return v.back();
}

}  // namespace

const char* FileFormatName(FileFormat f) noexcept {
    switch (f) {
        case FileFormat::Unknown: return "unknown";
        case FileFormat::Fasta:   return "fasta";
        case FileFormat::Fastq:   return "fastq";
        case FileFormat::Sam:     return "sam";
        case FileFormat::Bam:     return "bam";
    }
    return "unknown";
}

FileFormat SniffFormat(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return FileFormat::Unknown;

    std::array<unsigned char, 4> magic{};
    in.read(reinterpret_cast<char*>(magic.data()),
            static_cast<std::streamsize>(magic.size()));
    const std::streamsize got = in.gcount();
    if (got <= 0) return FileFormat::Unknown;

    // gzip magic — BAM is BGZF (gzip). We can't cheaply inflate to confirm
    // BAM vs gz-FASTA/FASTQ here, so fall back to the extension: only a
    // ".bam" suffix is reported as Bam. gz-FASTA/FASTQ stay Unknown and
    // ResolveMode applies the basename heuristic + GenomeReads default.
    if (got >= 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
        const std::string lower = ToLower(path);
        if (EndsWith(lower, ".bam")) return FileFormat::Bam;
        return FileFormat::Unknown;
    }

    const char c0 = static_cast<char>(magic[0]);
    if (c0 == '>') return FileFormat::Fasta;

    if (c0 == '@') {
        // SAM header lines are '@' + two uppercase letters (HD/SQ/PG/RG/CO/PG)
        // followed by a tab. FASTQ records are '@' + read name (rarely two
        // uppercase letters + tab). Disambiguate on that shape.
        if (got >= 3) {
            const char c1 = static_cast<char>(magic[1]);
            const char c2 = static_cast<char>(magic[2]);
            const bool two_upper =
                std::isupper(static_cast<unsigned char>(c1)) &&
                std::isupper(static_cast<unsigned char>(c2));
            if (two_upper) return FileFormat::Sam;
        }
        return FileFormat::Fastq;
    }

    // Raw (uncompressed) BAM begins with "BAM\1".
    if (got >= 4 && magic[0] == 'B' && magic[1] == 'A' && magic[2] == 'M' &&
        magic[3] == 1) {
        return FileFormat::Bam;
    }

    return FileFormat::Unknown;
}

FastaStats ComputeFastaStats(const std::string& path, std::uint64_t sample_limit) {
    FastaStats stats;
    std::ifstream in(path, std::ios::binary);
    if (!in) return stats;

    std::vector<std::uint64_t> lengths;  // sampled sequence lengths
    std::uint64_t cur_len = 0;
    bool in_seq = false;

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            // close previous sequence
            if (in_seq && lengths.size() < sample_limit) {
                lengths.push_back(cur_len);
            }
            ++stats.n_seqs;
            cur_len = 0;
            in_seq = true;
        } else if (in_seq) {
            cur_len += line.size();
        }
    }
    // close the final sequence
    if (in_seq && lengths.size() < sample_limit) {
        lengths.push_back(cur_len);
    }

    stats.sampled = lengths.size();
    stats.sampled_truncated = stats.n_seqs > stats.sampled;
    stats.median_len = Median(lengths);          // sorts lengths in place
    stats.n50 = N50(lengths);
    return stats;
}

namespace {

// Inspect SAM header (@-prefixed leading lines) for splice-aware program
// tokens. Returns true if any is found.
bool SamHeaderIsSpliceAware(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] != '@') break;  // header ends at first alignment line
        const std::string lower = ToLower(line);
        if (ContainsAny(lower, {"isoseq", "lima", "splice", "star", "cdna",
                                "flnc"})) {
            return true;
        }
    }
    return false;
}

bool BasenameLooksTranscript(const std::string& path) {
    const std::string base = ToLower(Basename(path));
    return ContainsAny(base, {"flnc", "isoseq", "cdna", "rna"});
}

SniffResult Make(FileFormat fmt, core::TranscriptMode mode, std::string reason) {
    SniffResult r;
    r.format = fmt;
    r.mode = mode;
    r.reason = std::move(reason);
    return r;
}

}  // namespace

SniffResult ResolveMode(const std::string& primary_path,
                        core::TranscriptMode mode_override,
                        bool has_reads,
                        bool has_assembly) {
    using core::TranscriptMode;

    // 1. Manual override always wins.
    if (mode_override != TranscriptMode::Auto) {
        SniffResult r = Make(FileFormat::Unknown, mode_override, "override");
        r.reason = std::string("override: --mode ") +
                   core::TranscriptModeName(mode_override);
        return r;
    }

    // 2. Both reads + assembly supplied → combined mode.
    if (has_reads && has_assembly) {
        return Make(FileFormat::Unknown, TranscriptMode::ReadsVsAssembly,
                    "both --reads and --assembly supplied");
    }

    // 3. Sniff the primary input.
    const FileFormat fmt = SniffFormat(primary_path);

    switch (fmt) {
        case FileFormat::Bam:
        case FileFormat::Sam: {
            if (SamHeaderIsSpliceAware(primary_path)) {
                return Make(fmt, TranscriptMode::Transcript,
                            "SAM/BAM @PG splice-aware token");
            }
            // TODO(v2): pilot-pass long-N CIGAR fraction (doc §3) needs a real
            // alignment pass — deferred. Default to genome reads for now.
            return Make(fmt, TranscriptMode::GenomeReads,
                        "SAM/BAM, no splice-aware @PG token");
        }

        case FileFormat::Fastq: {
            if (BasenameLooksTranscript(primary_path)) {
                return Make(fmt, TranscriptMode::Transcript,
                            "FASTQ basename matches FLNC/isoseq/cdna/rna");
            }
            return Make(fmt, TranscriptMode::GenomeReads, "FASTQ, DNA basename");
        }

        case FileFormat::Fasta: {
            FastaStats st = ComputeFastaStats(primary_path);
            // Assembly: long contigs, few of them.
            if (st.median_len > 50'000 && st.n50 > 100'000 && st.n_seqs < 5000) {
                SniffResult r = Make(fmt, TranscriptMode::Assembly,
                                     "FASTA median>50kb, N50>100kb, n_seqs<5k");
                r.fasta_stats = st;
                return r;
            }
            // iso-seq FLNC-as-FASTA signature: many medium-length records.
            if (st.median_len >= 300 && st.median_len <= 15'000 &&
                st.n_seqs > 50'000) {
                SniffResult r = Make(fmt, TranscriptMode::Transcript,
                                     "FASTA median 0.3-15kb, n_seqs>50k "
                                     "(FLNC-as-FASTA)");
                r.fasta_stats = st;
                return r;
            }
            // Basename hint even when length stats are inconclusive.
            if (BasenameLooksTranscript(primary_path)) {
                SniffResult r = Make(fmt, TranscriptMode::Transcript,
                                     "FASTA basename matches FLNC/isoseq/cdna/rna");
                r.fasta_stats = st;
                return r;
            }
            SniffResult r = Make(fmt, TranscriptMode::GenomeReads,
                                 "FASTA, default (low-confidence)");
            r.fasta_stats = st;
            return r;
        }

        case FileFormat::Unknown:
        default: {
            // gz-FASTQ/FASTA or unrecognised: still honour a transcript-looking
            // basename, else default to genome reads.
            if (BasenameLooksTranscript(primary_path)) {
                return Make(fmt, TranscriptMode::Transcript,
                            "unknown/gzip format, transcript basename");
            }
            return Make(fmt, TranscriptMode::GenomeReads,
                        "unknown format, default GenomeReads");
        }
    }
}

}  // namespace llmap::io
