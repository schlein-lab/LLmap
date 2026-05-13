#include "classical/minimizer_index.h"
#include "classical/minimizer_index_impl.h"

#include <fstream>

namespace llmap::classical {

bool MinimizerIndex::Save(const std::string& path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;

    // Magic + version
    uint32_t magic = 0x4D494E49;  // "MINI"
    uint32_t version = 1;
    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));

    // Config
    out.write(reinterpret_cast<const char*>(&config_), sizeof(config_));

    // Sequences
    uint32_t num_seqs = static_cast<uint32_t>(impl_->sequences.size());
    out.write(reinterpret_cast<const char*>(&num_seqs), sizeof(num_seqs));
    for (const auto& seq : impl_->sequences) {
        uint32_t name_len = static_cast<uint32_t>(seq.name.size());
        out.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
        out.write(seq.name.data(), name_len);
        out.write(reinterpret_cast<const char*>(&seq.length), sizeof(seq.length));
    }

    // Index entries
    uint64_t num_entries = impl_->index.size();
    out.write(reinterpret_cast<const char*>(&num_entries), sizeof(num_entries));
    for (const auto& [hash, positions] : impl_->index) {
        out.write(reinterpret_cast<const char*>(&hash), sizeof(hash));
        uint32_t num_pos = static_cast<uint32_t>(positions.size());
        out.write(reinterpret_cast<const char*>(&num_pos), sizeof(num_pos));
        for (const auto& [ref_id, pos] : positions) {
            out.write(reinterpret_cast<const char*>(&ref_id), sizeof(ref_id));
            out.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
        }
    }

    // Stats
    out.write(reinterpret_cast<const char*>(&impl_->stats), sizeof(impl_->stats));

    return out.good();
}

std::unique_ptr<MinimizerIndex> MinimizerIndex::Load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return nullptr;

    // Magic + version
    uint32_t magic, version;
    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (magic != 0x4D494E49 || version != 1) return nullptr;

    auto index = std::make_unique<MinimizerIndex>();

    // Config
    in.read(reinterpret_cast<char*>(&index->config_), sizeof(index->config_));

    // Sequences
    uint32_t num_seqs;
    in.read(reinterpret_cast<char*>(&num_seqs), sizeof(num_seqs));
    index->impl_->sequences.resize(num_seqs);
    for (auto& seq : index->impl_->sequences) {
        uint32_t name_len;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
        seq.name.resize(name_len);
        in.read(seq.name.data(), name_len);
        in.read(reinterpret_cast<char*>(&seq.length), sizeof(seq.length));
    }

    // Index entries
    uint64_t num_entries;
    in.read(reinterpret_cast<char*>(&num_entries), sizeof(num_entries));
    for (uint64_t i = 0; i < num_entries; ++i) {
        uint64_t hash;
        in.read(reinterpret_cast<char*>(&hash), sizeof(hash));
        uint32_t num_pos;
        in.read(reinterpret_cast<char*>(&num_pos), sizeof(num_pos));
        auto& positions = index->impl_->index[hash];
        positions.resize(num_pos);
        for (auto& [ref_id, pos] : positions) {
            in.read(reinterpret_cast<char*>(&ref_id), sizeof(ref_id));
            in.read(reinterpret_cast<char*>(&pos), sizeof(pos));
        }
    }

    // Stats
    in.read(reinterpret_cast<char*>(&index->impl_->stats), sizeof(index->impl_->stats));

    if (!in.good()) return nullptr;
    return index;
}

}  // namespace llmap::classical
