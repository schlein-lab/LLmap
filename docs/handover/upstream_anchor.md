# WORKTREE_CHANGES — IGH-mode upstream-anchor classifier

Branch: `feat/nahr-flanks-3318-linter`
Scope: `src/classify/`, `tests/classify/` (+ two single-line subdir
hookups in `src/CMakeLists.txt` and `tests/CMakeLists.txt`)
Author: classify-agent
Date: 2026-06-02

## TL;DR

(a) **API** — `llmap::classify::igh::Classify()` takes an `ISegDupCatalog`,
a `ReadInput` (sequence + mapping pos + optional anchor offset), and an
optional `ClassifyOptions::hypothesis`. Returns
`{haplotype_class_call, promoter_motif_match, flanking_kmer_mismatches,
promoter_test_ran, flank_window_used}`. Motif test is binary
canonical/chimdup with `Mixed`/`None` edge cases.

(b) **Test results** — 16/16 tests pass under
`build/tests/classify/test_upstream_anchor` (10 utility, 6 end-to-end
including the three target cases: canonical-promoter / chimdup-promoter /
CH1-only).

(c) **LOC** — well under the 400-LOC soft cap:
   - `src/classify/segdup_catalog_iface.h` — 124
   - `src/classify/upstream_anchor.h`       — 157
   - `src/classify/upstream_anchor.cpp`     — 277
   - `tests/classify/test_upstream_anchor.cpp` — 352

(d) **Catalog-loader interface** — defined as the pure-virtual
`ISegDupCatalog` in `src/classify/segdup_catalog_iface.h`. The
classifier holds only a `const ISegDupCatalog&`. The real loader (parallel
work in `src/catalog/`) will satisfy this interface by exposing
`LookupByPosition(target_id, pos)` and `LookupByHaplotypeClass(tag)`,
each returning `std::optional<CatalogEntryView>`. The view exposes
`locus_id`, `haplotype_class`, `mapping_primary.{kmer_size,max_mismatch,
include_flanking_bp,include_flanking_anchor}`, and the full
`promoter_signature` block (window, anchor, canonical_motif,
duplicate_motif). Tests use an in-memory `MockCatalog` that implements the
interface — same code path as the real loader will hit.

## Files added

```
src/classify/CMakeLists.txt
src/classify/segdup_catalog_iface.h
src/classify/upstream_anchor.h
src/classify/upstream_anchor.cpp
tests/classify/CMakeLists.txt
tests/classify/test_upstream_anchor.cpp
```

## Files modified (single-line subdirectory hooks)

- `src/CMakeLists.txt` — `add_subdirectory(classify)` (one stanza, sits
  next to the existing `add_subdirectory(mapper)` added by the fallback
  agent).
- `tests/CMakeLists.txt` — `add_subdirectory(classify)` at end-of-file.

No other files outside `src/classify/` or `tests/classify/` were touched.

## Build + run

```
cd build
cmake -DLLMAP_ENABLE_TESTS=ON \
      -DLLMAP_ENABLE_CUDA=OFF \
      -DLLMAP_ENABLE_FOUNDATION=OFF \
      -DLLMAP_ENABLE_FAISS=OFF \
      -DLLMAP_ENABLE_CLAUDE=OFF ..
cmake --build . --target test_upstream_anchor
./tests/classify/test_upstream_anchor
# → 16 tests pass, 0 fail
cmake --build . --target llmap     # full binary still links
```

## Behaviour notes / contract

1. **No hard-coded motifs.** Canonical / duplicate sequences are pulled
   verbatim from
   `diagnostic_features.promoter_signature.canonical_motif` and
   `.duplicate_motif`. The classifier honours the literal `...` ellipsis
   in the catalog motif strings as an "ordered-subsequence" wildcard
   (prefix then suffix in 5'->3' order; spacer length unconstrained
   within the catalog window). This matches the curated
   `IGHG4_chimdup_tandem.json` syntax.

2. **Window selection.** Defaults to
   `[-include_flanking_bp, 0]` relative to the catalog anchor. If the
   entry has a `promoter_signature.window_relative_to_anchor` (e.g.
   `"-100..0"`), that wins. Reads with anchor at offset 0 (no upstream
   bases available) set `flank_window_used=false`.

3. **K-mer mismatch score.** Counted over the read's upstream window vs
   the canonical motif's literal content (ellipsis bytes stripped). Per
   k-mer contribution is clamped by `mapping_primary.max_mismatch` so a
   hot-spot region cannot dominate. Score is 0 for an exact canonical
   match, positive otherwise.

4. **Hypothesis path.** When the upstream caller knows the candidate
   haplotype class (e.g. PSV pre-classification), `opts.hypothesis`
   bypasses the position lookup. Useful when the read maps to a position
   that is itself canonical/chimdup-ambiguous.

5. **Concurrency.** `Classify()` is stateless; safe to call from
   multiple threads as long as the catalog is itself thread-safe for
   reads (immutable catalogs are the expected case).

## Hand-off to the catalog-loader agent

The catalog loader in `src/catalog/` should:

1. Implement `llmap::classify::igh::ISegDupCatalog` (either directly on
   the loader class or via a small adapter).
2. Populate `CatalogEntryView` fields verbatim from the curated JSON
   schema; the field names are 1:1 with the schema's terminal keys.
3. Index entries by `(target_id, [start,end))` for `LookupByPosition`
   and by `haplotype_class` for `LookupByHaplotypeClass`.

Once the real loader is wired into the CLI, the `llmap_classify` target
needs no changes — the classifier already depends only on the abstract
interface.

## Out of scope (explicit non-goals for this change)

- No CLI subcommand (this is a library-level component only).
- No alignment-record integration (the classifier returns a
  `ClassificationResult`; downstream BAM/Parquet wiring lives in the
  output module and stays untouched here).
- No catalog JSON parser (the loader agent owns that).
- No reverse-strand handling beyond a `is_reverse` field on
  `ReadInput`; the actual reverse-complementation contract will be
  defined when the loader provides the reference flanking sequence.
