# Catalog Loader Worktree Changes

Branch: `feat/nahr-flanks-3318-linter`
Scope: `src/catalog/`, `tests/catalog/`, root `CMakeLists.txt` (2-line subdir hookup)
Status: implemented + built + tested locally (Debug, no CUDA/FAISS/ONNX)

## Files added

### `src/catalog/` (module)
- `segdup_catalog.h`           — public API, 221 LOC
- `segdup_catalog.cpp`         — catalog container, indices, lookup, 225 LOC
- `segdup_catalog_parse.cpp`   — curated T1 JSON parser (nlohmann/json), 215 LOC
- `bulk_loader.cpp`            — T2 JSONL streaming loader, 126 LOC
- `CMakeLists.txt`             — `llmap_catalog` STATIC target, find_package(nlohmann_json) with FetchContent fallback, 46 LOC

### `tests/catalog/` (gtest suites)
- `test_curated_loader.cpp`    — 7 tests, loads the curated dir + IGHG entry assertions, 116 LOC
- `test_coord_lookup.cpp`      — 8 tests, single-point + overlap + boundary lookups, 127 LOC
- `test_bulk_loader.cpp`       — 3 tests (1 skip if canon JSONL absent), 110 LOC
- `CMakeLists.txt`             — 3 gtest_discover_tests targets, 36 LOC

## Root file edited

`CMakeLists.txt` — 5 lines inserted:
```diff
+# SegDup catalog loader — independent of src/CMakeLists.txt to keep the
+# catalog module self-contained (own nlohmann_json dep).
+add_subdirectory(src/catalog)
...
+    add_subdirectory(tests/catalog)
```
No other source file in the repo was touched.

## Public API surface

```cpp
namespace llmap::catalog {

struct GenomicCoords { ... };          // assembly, chrom, [start,end), strand
struct DiagnosticSnp { ... };          // T1 discriminating-SNP record
struct MappingPrimary { ... };
struct MappingFallbackStage { ... };
struct SegDupCatalogEntry {
    enum class Tier { T1_Curated, T2_Bulk };
    std::string locus_id, human_name, version, schema_version;
    std::string structural_architecture, haplotype_class, nahr_status;
    std::vector<std::string> mechanism, clinical_function;
    std::unordered_map<std::string,GenomicCoords> coords_by_assembly;
    std::vector<DiagnosticSnp> discriminating_snps;
    MappingPrimary mapping_primary;
    std::vector<MappingFallbackStage> fallback_chain;
    Tier tier;
    bool contains(assembly, chrom, pos) const;
    bool overlaps(assembly, chrom, start, end) const;
};

class SegDupCatalog {
  LoadStatus load_curated_file(path);
  LoadStatus load_curated_dir(dir);
  LoadStatus load_bulk_jsonl(path);
  void       add_entry(SegDupCatalogEntry);
  std::optional<SegDupCatalogEntry>
             lookup_by_coords(assembly, chrom, pos) const;     // most-specific
  std::vector<SegDupCatalogEntry>
             lookup_overlapping(assembly, chrom, start, end) const;
  std::vector<SegDupCatalogEntry>
             lookup_by_locus_id(id) const;
  std::size_t count_architecture(architecture) const;
};

std::optional<SegDupCatalogEntry>
parse_curated_json(text, source_path, err);
std::optional<SegDupCatalogEntry>
parse_bulk_jsonl_line(text, source_path, err);

}  // namespace llmap::catalog
```

## Design notes

1. **Most-specific resolution.** Three curated entries (`IGHG_canondup_nahr_block`,
   `IGHG4_chimdup_tandem`, `IGHG4_chimdup_canonical_arch`) legitimately overlap at
   chr14:105625900. `lookup_by_coords` returns the smallest-span entry; the
   `lookup_overlapping` API returns the full set. Tests cover both behaviours.
2. **nlohmann/json.** Module-local dependency via `find_package(nlohmann_json
   3.10 QUIET)` first, falling back to FetchContent at v3.11.3. No edit to
   `src/CMakeLists.txt`, no `third_party/` writes.
3. **Schema tolerance.** Files in `catalog/curated/` advertise `schema_version
   0.1`; the parser tolerates any 0.x because v0.2 adds fields only. Unknown
   keys are ignored. Missing optional fields fall back to neutral defaults.
4. **LOC budget.** Largest .cpp = 225 (segdup_catalog.cpp), well under the
   400-LOC soft cap.
5. **Build isolation.** `add_subdirectory(src/catalog)` from the root keeps
   the module self-contained and does not require changes to
   `src/CMakeLists.txt` (which other parallel agents are touching).

## Build + test result (local, Debug)

Configured in fresh `build_catalog/`:
```
cmake .. -DLLMAP_ENABLE_CUDA=OFF -DLLMAP_ENABLE_FAISS=OFF \
         -DLLMAP_ENABLE_FOUNDATION=OFF -DLLMAP_USE_NATIVE_ARCH=OFF \
         -DCMAKE_BUILD_TYPE=Debug \
         -DCMAKE_PREFIX_PATH=/home/<user>/miniconda3
make llmap_catalog test_curated_loader test_coord_lookup test_bulk_loader
```

`ctest --output-on-failure` from `build_catalog/tests/catalog/`:
- 17 passed
- 1 skipped (`BulkLoaderTest.LoadsCanonicalBulkIfPresent` — no
  `catalog/bulk/v2026.Q2.grch38.jsonl` shipped yet)
- 0 failed

Total `make` time on the lib + 3 tests ≈ 6 s.

## Known follow-ups (not part of this PR)

- Ship `catalog/bulk/v2026.Q2.grch38.jsonl` — currently empty dir. The
  bulk loader is fully tested with synthetic JSONL; the canonical file
  test is gated by `GTEST_SKIP()` until the file exists.
- The mapper / classifier should consume `lookup_overlapping` rather
  than `lookup_by_coords` when deciding fallback strategy, because
  curated entries intentionally nest.
- An interval-tree backing for `by_chrom_` becomes worthwhile once T2
  bulk crosses ~10⁵ records; current linear sweep with start-sorted
  buckets is adequate for the catalog sizes shipped today.

## Out-of-scope (untouched)

No edits to `src/classify/`, `src/mapper/`, any other `src/*`, or any
existing test file. The benchmark-data deletions visible in `git status`
predate this work session.
