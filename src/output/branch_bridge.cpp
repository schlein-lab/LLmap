// LLmap — BRANCH-compatible GAF writer.

#include "output/branch_bridge.h"

#include <fstream>
#include <sstream>

namespace llmap::output {

std::string EmitBranchGafRow(const BranchGafRow& r) {
    std::ostringstream os;
    os  << r.read_id            << '\t'
        << (r.cluster_id == 0 ? std::string{"-"}
                                : std::to_string(r.cluster_id)) << '\t'
        << r.haplotype          << '\t'
        << r.path               << '\t'
        << r.score              << '\t'
        << r.n_supporting_reads << '\t'
        << r.provenance;
    return os.str();
}

bool WriteBranchGafFile(const std::filesystem::path& path,
                         std::span<const BranchGafRow> rows) {
    std::ofstream out(path);
    if (!out) return false;
    for (const auto& r : rows) {
        out << EmitBranchGafRow(r) << '\n';
    }
    return static_cast<bool>(out);
}

}  // namespace llmap::output
