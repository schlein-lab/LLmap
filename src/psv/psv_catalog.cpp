// LLmap — PSV catalog implementation.

#include "psv/psv_catalog.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>

namespace llmap::psv {

void PsvCatalog::AddSite(PsvSite site) {
    const auto idx = sites_.size();
    id_index_[site.psv_id] = idx;
    pos_index_[site.chrom][site.position] = idx;
    sites_.push_back(std::move(site));
    indexed_ = false;
}

std::vector<const PsvSite*> PsvCatalog::GetSitesInRegion(
    std::string_view chrom,
    std::uint64_t start,
    std::uint64_t end) const {

    std::vector<const PsvSite*> result;
    const std::string chrom_str{chrom};

    auto chrom_it = sorted_positions_.find(chrom_str);
    if (chrom_it == sorted_positions_.end()) {
        return result;
    }

    const auto& positions = chrom_it->second;
    auto it_start = std::lower_bound(positions.begin(), positions.end(), start);
    auto it_end = std::upper_bound(positions.begin(), positions.end(), end);

    auto pos_it = pos_index_.find(chrom_str);
    if (pos_it == pos_index_.end()) {
        return result;
    }

    for (auto it = it_start; it != it_end; ++it) {
        auto idx_it = pos_it->second.find(*it);
        if (idx_it != pos_it->second.end()) {
            result.push_back(&sites_[idx_it->second]);
        }
    }

    return result;
}

const PsvSite* PsvCatalog::GetSiteById(std::uint64_t psv_id) const {
    auto it = id_index_.find(psv_id);
    if (it == id_index_.end()) {
        return nullptr;
    }
    return &sites_[it->second];
}

const PsvSite* PsvCatalog::GetSiteAtPosition(
    std::string_view chrom,
    std::uint64_t position) const {

    const std::string chrom_str{chrom};
    auto chrom_it = pos_index_.find(chrom_str);
    if (chrom_it == pos_index_.end()) {
        return nullptr;
    }

    auto pos_it = chrom_it->second.find(position);
    if (pos_it == chrom_it->second.end()) {
        return nullptr;
    }

    return &sites_[pos_it->second];
}

std::vector<std::string> PsvCatalog::GetParalogs() const {
    std::set<std::string> paralogs;
    for (const auto& site : sites_) {
        for (const auto& [paralog_id, _] : site.paralog_alleles) {
            paralogs.insert(paralog_id);
        }
    }
    return {paralogs.begin(), paralogs.end()};
}

std::vector<std::string> PsvCatalog::GetChromosomes() const {
    std::set<std::string> chroms;
    for (const auto& site : sites_) {
        chroms.insert(site.chrom);
    }
    return {chroms.begin(), chroms.end()};
}

void PsvCatalog::Clear() noexcept {
    sites_.clear();
    id_index_.clear();
    pos_index_.clear();
    sorted_positions_.clear();
    indexed_ = false;
}

void PsvCatalog::BuildIndex() {
    sorted_positions_.clear();

    for (const auto& [chrom, pos_map] : pos_index_) {
        auto& positions = sorted_positions_[chrom];
        positions.reserve(pos_map.size());
        for (const auto& [pos, _] : pos_map) {
            positions.push_back(pos);
        }
        std::sort(positions.begin(), positions.end());
    }

    indexed_ = true;
}

std::optional<PsvCatalog> LoadPsvCatalogFromBed(
    const std::filesystem::path& path) {

    std::ifstream file(path);
    if (!file) {
        return std::nullopt;
    }

    PsvCatalog catalog;
    std::string line;
    std::uint64_t psv_id = 0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string chrom;
        std::uint64_t pos;
        char ref;
        std::string alleles_str;
        float informativeness = 1.0f;

        if (!(iss >> chrom >> pos >> ref >> alleles_str)) {
            continue;
        }
        iss >> informativeness;

        PsvSite site;
        site.psv_id = psv_id++;
        site.chrom = chrom;
        site.position = pos;
        site.ref_allele = ref;
        site.informativeness = informativeness;

        std::istringstream alleles_iss(alleles_str);
        std::string pair;
        while (std::getline(alleles_iss, pair, ',')) {
            auto colon = pair.find(':');
            if (colon != std::string::npos && colon + 1 < pair.size()) {
                std::string paralog = pair.substr(0, colon);
                char allele = pair[colon + 1];
                site.paralog_alleles[paralog] = allele;
            }
        }

        if (!site.paralog_alleles.empty()) {
            catalog.AddSite(std::move(site));
        }
    }

    catalog.BuildIndex();
    return catalog;
}

std::optional<PsvCatalog> LoadPsvCatalogFromVcf(
    const std::filesystem::path& path) {

    std::ifstream file(path);
    if (!file) {
        return std::nullopt;
    }

    PsvCatalog catalog;
    std::string line;
    std::uint64_t psv_id = 0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string chrom;
        std::uint64_t pos;
        std::string id, ref, alt, qual, filter, info;

        if (!(iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            continue;
        }

        if (ref.size() != 1 || alt.size() != 1) {
            continue;
        }

        PsvSite site;
        site.psv_id = psv_id++;
        site.chrom = chrom;
        site.position = pos;
        site.ref_allele = ref[0];

        auto psv_info_pos = info.find("PSV=");
        if (psv_info_pos != std::string::npos) {
            auto start = psv_info_pos + 4;
            auto end = info.find(';', start);
            std::string psv_val = info.substr(start, end - start);

            std::istringstream psv_iss(psv_val);
            std::string pair;
            while (std::getline(psv_iss, pair, ',')) {
                auto colon = pair.find(':');
                if (colon != std::string::npos && colon + 1 < pair.size()) {
                    std::string paralog = pair.substr(0, colon);
                    char allele = pair[colon + 1];
                    site.paralog_alleles[paralog] = allele;
                }
            }
        }

        site.informativeness = ComputeInformativeness(site);

        if (!site.paralog_alleles.empty()) {
            catalog.AddSite(std::move(site));
        }
    }

    catalog.BuildIndex();
    return catalog;
}

bool SavePsvCatalogToBed(
    const PsvCatalog& catalog,
    const std::filesystem::path& path) {

    std::ofstream file(path);
    if (!file) {
        return false;
    }

    file << "# PSV catalog: chrom\tpos\tref\tparalog_alleles\tinformativeness\n";

    for (const auto& site : catalog.Sites()) {
        file << site.chrom << '\t' << site.position << '\t' << site.ref_allele << '\t';

        bool first = true;
        for (const auto& [paralog, allele] : site.paralog_alleles) {
            if (!first) file << ',';
            file << paralog << ':' << allele;
            first = false;
        }

        file << '\t' << site.informativeness << '\n';
    }

    return true;
}

float ComputeInformativeness(const PsvSite& site) {
    if (site.paralog_alleles.empty()) {
        return 0.0f;
    }

    std::unordered_map<char, int> allele_counts;
    for (const auto& [_, allele] : site.paralog_alleles) {
        allele_counts[allele]++;
    }

    if (allele_counts.size() == 1) {
        return 0.0f;
    }

    const auto n = static_cast<float>(site.paralog_alleles.size());
    float entropy = 0.0f;
    for (const auto& [_, count] : allele_counts) {
        const float p = static_cast<float>(count) / n;
        if (p > 0.0f) {
            entropy -= p * std::log2(p);
        }
    }

    const float max_entropy = std::log2(n);
    return max_entropy > 0.0f ? entropy / max_entropy : 0.0f;
}

void PsvCatalogBuilder::AddParalogSequence(
    std::string paralog_id,
    std::string_view sequence) {

    paralogs_[paralog_id] = std::string{sequence};
}

void PsvCatalogBuilder::SetReferenceSequence(std::string_view sequence) {
    reference_ = std::string{sequence};
}

PsvCatalog PsvCatalogBuilder::Build(
    std::string_view chrom,
    std::uint64_t start_position,
    float min_informativeness) {

    PsvCatalog catalog;

    if (reference_.empty() || paralogs_.empty()) {
        return catalog;
    }

    std::size_t min_len = reference_.size();
    for (const auto& [_, seq] : paralogs_) {
        min_len = std::min(min_len, seq.size());
    }

    std::uint64_t psv_id = 0;
    for (std::size_t i = 0; i < min_len; ++i) {
        char ref_base = reference_[i];
        if (!IsValidAllele(ref_base)) {
            continue;
        }

        std::unordered_map<std::string, char> alleles;
        bool has_variation = false;

        for (const auto& [paralog_id, seq] : paralogs_) {
            char allele = seq[i];
            if (!IsValidAllele(allele)) {
                continue;
            }
            alleles[paralog_id] = allele;
            if (allele != ref_base) {
                has_variation = true;
            }
        }

        if (has_variation && !alleles.empty()) {
            PsvSite site;
            site.psv_id = psv_id++;
            site.chrom = std::string{chrom};
            site.position = start_position + i;
            site.ref_allele = ref_base;
            site.paralog_alleles = std::move(alleles);
            site.informativeness = ComputeInformativeness(site);

            if (site.informativeness >= min_informativeness) {
                catalog.AddSite(std::move(site));
            }
        }
    }

    catalog.BuildIndex();
    return catalog;
}

}  // namespace llmap::psv
