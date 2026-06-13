// LLmap — Contamination / provenance spectrum aggregator implementation.

#include "provenance/contamination_spectrum.h"

#include <cstdio>
#include <sstream>

namespace llmap::provenance {

namespace {
// The 7 bioconfounder bits in declaration order, for tallying + reporting.
constexpr BioConfounder kBioFlags[] = {
    BioConfounder::Vdj,            BioConfounder::Shm,
    BioConfounder::ClassSwitch,   BioConfounder::GeneConversion,
    BioConfounder::MtHeteroplasmy, BioConfounder::Mosaicism,
    BioConfounder::Imprinting,    BioConfounder::SvSpanning,
};
}  // namespace

void ContaminationSpectrum::Add(const ReadProvenance& rp) noexcept {
    const auto i = static_cast<std::size_t>(rp.origin);
    if (i >= kN) return;
    ++n_reads_[i];
    bases_[i] += rp.aligned_bases;
    post_sum_[i] += rp.posterior;
    ++total_;
    // Layer-3 overlay: count each set bioconfounder bit (does NOT touch total_).
    for (const BioConfounder f : kBioFlags) {
        if (rp.bioconfounder & static_cast<std::uint16_t>(f)) {
            const auto b = static_cast<std::size_t>(f);
            // map bit value (1<<k) to a dense index via the bit position
            std::size_t idx = 0, v = b;
            while (v > 1) { v >>= 1; ++idx; }
            if (idx < kBio) ++bio_[idx];
        }
    }
}

ClassStat ContaminationSpectrum::Stat(ProvenanceClass c) const noexcept {
    const auto i = static_cast<std::size_t>(c);
    ClassStat s;
    if (i >= kN) return s;
    s.n_reads = n_reads_[i];
    s.bases = bases_[i];
    s.fraction = total_ ? static_cast<double>(n_reads_[i]) / static_cast<double>(total_) : 0.0;
    s.mean_posterior = n_reads_[i] ? post_sum_[i] / static_cast<double>(n_reads_[i]) : 0.0;
    return s;
}

std::uint64_t ContaminationSpectrum::BioConfounderReads(BioConfounder f) const noexcept {
    std::size_t idx = 0, v = static_cast<std::size_t>(f);
    if (v == 0) return 0;
    while (v > 1) { v >>= 1; ++idx; }
    return idx < kBio ? bio_[idx] : 0;
}

std::string ContaminationSpectrum::ToString() const {
    std::ostringstream os;
    // Layer-1 partition (Σ == total).
    os << "# layer1_origin_partition (sum == n_input)\n";
    os << "class\tn_reads\tfraction\tbases\tmean_posterior\n";
    for (std::size_t i = 0; i < kN; ++i) {
        if (n_reads_[i] == 0) continue;
        const auto c = static_cast<ProvenanceClass>(i);
        const ClassStat s = Stat(c);
        char fbuf[64], pbuf[64];
        std::snprintf(fbuf, sizeof fbuf, "%.6f", s.fraction);
        std::snprintf(pbuf, sizeof pbuf, "%.4f", s.mean_posterior);
        os << ProvenanceClassTag(c) << '\t' << s.n_reads << '\t' << fbuf << '\t'
           << s.bases << '\t' << pbuf << '\n';
    }
    // Layer-3 overlay (orthogonal; may overlap Host — NOT part of the partition).
    bool any_bio = false;
    for (std::size_t k = 0; k < kBio; ++k) if (bio_[k]) any_bio = true;
    if (any_bio) {
        os << "# layer3_bioconfounder_overlay (co-occurs with origin, not summed)\n";
        os << "flag\tn_reads\n";
        for (const BioConfounder f : kBioFlags) {
            const std::uint64_t n = BioConfounderReads(f);
            if (n) os << BioConfounderTag(f) << '\t' << n << '\n';
        }
    }
    return os.str();
}

bool ContaminationSpectrum::WriteTsv(const std::string& path) const {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    const std::string s = ToString();
    const bool ok = std::fwrite(s.data(), 1, s.size(), f) == s.size();
    std::fclose(f);
    return ok;
}

}  // namespace llmap::provenance
