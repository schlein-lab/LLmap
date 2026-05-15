// LLmap -- on-disk annotation store (.regann files).
//
// Format: BED-like, one interval per line.
//   ref_name<TAB>start<TAB>end<TAB>region_name<TAB>source<TAB>layer<TAB>params
// where params is a comma-separated list of key=value pairs covering any of
// the ParamOverride fields.
//
// Example:
//   chr14   105500000   107300000   paralog_family_immunoglobulin
//           taxonomy   1   max_occ=500,lambda_scale=2.0,report_multi_position=1
//
// (the line is on one line in the file; broken here for documentation only)
//
// We keep this human-readable so users can grep / hand-edit if needed.

#pragma once

#include "annot/annot_types.h"
#include "annot/interval_tree.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace llmap::annot {

class AnnotationStore {
public:
    AnnotationStore() = default;

    // Build from in-memory intervals plus the contig name list (so we can
    // map ref_id <-> name on serialisation).
    static std::unique_ptr<AnnotationStore> Create(
        std::vector<AnnotationInterval> intervals,
        std::vector<std::string> contig_names);

    // Load from a .regann file. Returns nullptr on parse failure.
    static std::unique_ptr<AnnotationStore> Load(
        const std::filesystem::path& path,
        const std::vector<std::string>& contig_names);

    // Save to a .regann file.
    bool Save(const std::filesystem::path& path) const;

    const AnnotationIndex& Index() const { return index_; }
    size_t Size() const { return index_.Size(); }
    bool Empty() const { return index_.Empty(); }

    // Convenience: composed params at a position with caller-supplied defaults.
    ParamOverride ParamsAt(uint32_t ref_id, uint32_t pos,
                           const ParamOverride& defaults = {}) const {
        return index_.ParamsAt(ref_id, pos, defaults);
    }

    const std::vector<std::string>& ContigNames() const { return contig_names_; }

private:
    AnnotationIndex index_;
    std::vector<std::string> contig_names_;
};

}  // namespace llmap::annot
