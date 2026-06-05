# LLmap Mode-6 — Taxonomic Binning of Un-Human Reads

**Status:** v0 scaffold on branch `mode6-taxbin`. Engine + CLI + registration land;
in-process Stage-1/cascade integration is the next step.

## Motivation

The genome-wide hidden-CNV pipeline produces a residual: reads (or read parts)
that map to **neither GRCh38 nor CHM13**. This residual is two things mixed:
1. **collapsed human** — sequence the references show once but the individual
   carries ≥2× (the "dreadful dozen" deeper layer), and
2. **non-human** — EBV (ubiquitous in EBV-transformed LCLs), mycoplasma, viral,
   and other contamination, plus genuinely novel sequence.

Classifying (2) with Kraken2/BLAST means a fragile external dependency and a
huge `nt` database. Mode-6 does it **LLmap-natively** (LLmap-only mandate),
reusing machinery LLmap already has, and emitting **honest, explicit
uncertainty** rather than a single brittle best-hit.

## Pipeline

```
un-human reads
   │
   ├─[Stage 1] self_interference::AllpairPipeline
   │     embed → FAISS k-NN similarity graph → Leiden clusters → SelfWaveCollapse
   │     ⇒ per-read cluster_id   (reads that look alike are grouped BEFORE any reference)
   │
   ├─[Cascade] for each species in panel (specific → broad):
   │     classical::ClassicalPipeline(species_index).AlignReads(reads)
   │     ⇒ per-read confidence = primary alignment identity ∈ [0,1]
   │     reads that bind strongly short-circuit; only the unbound residual is
   │     carried to the next, broader reference ("keep mapping until it finds something")
   │     ⇒ per-read × per-species confidence matrix
   │
   └─[Engine] classify::Taxbin
         1. per-read raw likelihood over (species…, NOVEL):
              raw[s] = conf[s];  raw[NOVEL] = max(0, bind_threshold − max_s conf[s]);  normalize
         2. cluster aggregate = Σ member raw vectors, normalized
              → the cluster's collective mass decides its consensus species
         3. per-read posterior = (1−w)·own + w·cluster_aggregate, normalized
         4. hard call = argmax posterior  (NOVEL is just the last index)
              flag `by_cluster` when the cluster changed the call,
              flag `novel` when nothing in the panel binds
   ⇒ per-read: full likelihood vector + call + margin   (deterministic, posthoc)
   ⇒ per-cluster: consensus species + purity
```

This is exactly the user's design: *likelihood bindings to several species, then
posthoc deterministically (a) output the uncertainty and (b) collapse into one
species based on the cluster of read-to-read comparisons.*

## What exists vs. what is new

| Component | Source | Status |
|-----------|--------|--------|
| read-to-read similarity + Leiden clustering | `self_interference::{SimilarityGraph,LeidenClustering}` | exists |
| intra-cluster EM collapse | `self_interference::SelfWaveCollapse` | exists |
| reference alignment / confidence | `classical::ClassicalPipeline` | exists |
| multi-reference EM (optional richer backend) | `reference_collapse::Stage2Pipeline` | exists |
| **panel cascade + likelihood + cluster-consensus collapse** | `classify::Taxbin` | **new (this branch)** |
| `llmap taxbin` CLI | `cli/cmd_taxbin.cpp` | **new (this branch)** |

## CLI (v0)

```
llmap taxbin --conf conf.tsv [--clusters clusters.tsv] -o calls.tsv \
             [--clusters-out clusters_summary.tsv] \
             [--bind-threshold 0.80] [--cluster-weight 0.50] [--min-margin 0.10]
```

`conf.tsv`: header `read_id  sp1  sp2 …`, one row per read with confidences ∈[0,1].
`clusters.tsv`: `read_id  cluster_id` (from `llmap allpair`).

Output columns: `read_id, cluster_id, call, top_prob, margin, novel, by_cluster,
L_<species>…, L_NOVEL`.

## Next steps (for the autonomous build to complete)

1. **In-process cascade**: drive `classical::ClassicalPipeline` per panel entry
   inside `cmd_taxbin`, short-circuiting bound reads, to populate `TaxbinInput`
   directly (drop the `--conf` TSV path; keep it as a debug ingress).
2. **In-process Stage 1**: call `AllpairPipeline` for `cluster_ids` (needs the
   ONNX embedding model). Degrade gracefully to per-read = own-cluster when no
   model is supplied (engine already handles empty `cluster_ids`).
3. **Panel manager**: `--panel panel.tsv` (`label  index_path  cascade_rank`);
   build species reference indices (LLmap already indexes mouse/rat/apes/
   zebrafish/dros/celegans/yeast — extend with EBV/mycoplasma/viral panels).
4. **NOVEL handoff**: emit the `novel`/`by_cluster` reads as FASTA for the
   genuinely-dark bucket.
5. Unit tests under `tests/classify/` for `Taxbin::RawLikelihood` and
   cluster-consensus collapse (engine is pure/deterministic — easy to test).
```
