// LLmap — Pangenome (GAF) bridge implementation.

#include "output/pangenome_bridge.h"

#include <fstream>
#include <sstream>

namespace llmap::output {

std::string EmitGafRow(const GafRow& r) {
    std::ostringstream os;
    os  << r.query_name        << '\t'
        << r.query_length      << '\t'
        << r.query_start       << '\t'
        << r.query_end         << '\t'
        << r.strand            << '\t'
        << r.path              << '\t'
        << r.path_length       << '\t'
        << r.path_start        << '\t'
        << r.path_end          << '\t'
        << r.residue_matches   << '\t'
        << r.block_length      << '\t'
        << static_cast<int>(r.mapping_quality);
    return os.str();
}

bool WriteGafFile(const std::filesystem::path& path,
                  std::span<const GafRow> rows) {
    std::ofstream out(path);
    if (!out) return false;
    for (const auto& r : rows) {
        out << EmitGafRow(r) << '\n';
    }
    return static_cast<bool>(out);
}

}  // namespace llmap::output
