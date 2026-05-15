# Extending LLmap's region-knowledge layer

LLmap ships with curated knowledge modules for human (`organisms/human/`) and
mouse (`organisms/mouse/`). Everything else — other vertebrates, plants,
fungi, bacteria, viruses, synthetic / fictive DNA — works the same way: you
write a module describing the regions your genome contains, and the runtime
mapper consumes it.

This document walks through writing a new module from scratch.

## Three layers of region knowledge

LLmap consumes region knowledge in three increasingly specific layers. Each
layer is optional. Mapping always works without any of them (classical
seed-chain-extend); each layer activated gains capability.

| Layer | What | When used | Fallback if absent |
|-------|------|-----------|---------------------|
| **1. Universal taxonomy** | feature → region type via deterministic rules (this file's main subject) | mapping-time, every read | classical scoring |
| **2. Per-locus database** (`specific_loci/`) | named real-world loci (each centromere, each known SD, each named cluster) with their coordinates, expected paralog structure, known PSVs | mapping-time when a read lands in a flagged locus | falls back to Layer 1 |
| **3. Active agent** (`--llm=auto`) | live LLM consultation at mapping checkpoints, with tool access (bash, web fetch, region lookup, PSV check, local grep) | only at chain-ambiguity, paralog uncertainty, SD-expansion, or unknown-region checkpoints | falls back to Layer 2 → Layer 1 |

The three modes for Layer 3:
- `--llm=auto` (default): use the agent when available, silently fall back otherwise
- `--llm=off`: never invoke the agent — fully deterministic, fully offline
- `--llm=required`: fail loudly if the agent isn't reachable (use this for reproducible review runs)

When the agent IS active, every read whose mapping was influenced by it
gets an `XL:Z:llm` BAM tag plus a human-readable `XR:Z:<reason>` and (where
relevant) `XW:Z:<wave>` multi-position output. The end-of-run summary
breaks down how many decisions came from the agent and what classes of
finding the agent contributed that classical alone could not have made.

---

## Step 1 — Decide whether to use the universal taxonomy

The universal taxonomy in `schema.json` covers concepts that apply to any
genome: `unique_single_copy`, `low_complexity`, `tandem_repeat`,
`paralog_family`, `dispersed_repeat`, `terminal_repeat`, `centromere_like`,
`coding`, `pseudogene`, plus an escape hatch `custom`.

Most projects will reuse these. You only need `custom` entries for genome
types where the universal categories don't fit — for example:

- **Viruses:** packaging signals, gag/pol/env genes, regulatory motifs.
- **Synthetic DNA:** plasmid origin of replication, selection marker,
  spacer regions, CRISPR target site.
- **Fictive DNA** (educational / testing data): anything you invent.

When in doubt, start with the universal taxonomy and add `custom` entries
only when you need them.

---

## Step 2 — Inventory your region types

List the region types that matter for your organism / question.

Example for plants (Arabidopsis):

| Region | Universal ID | Notes |
|--------|--------------|-------|
| Centromere (CEN178 satellite) | `centromere_like` | 178 bp tandem |
| 45S rDNA repeat | `tandem_repeat` | 8–10 kb units |
| Transposable element (Athila LTR) | `dispersed_repeat` | LTR retroelement |
| Knob heterochromatin | `low_complexity` | low entropy |
| Coding gene | `coding` | average ORF |
| Disease-resistance gene cluster (NB-LRR) | `paralog_family` | high intra-cluster homology |

Example for SARS-CoV-2:

| Region | Universal ID | Notes |
|--------|--------------|-------|
| 5' UTR + leader | `terminal_repeat` | TRS-leader at 5' |
| ORF1ab | `coding` | overlapping ORFs |
| Spike RBD | `coding` (custom subtype "high_variation") | rapid evolution |
| 3' UTR + poly-A | `terminal_repeat` | |

Example for synthetic plasmid pUC19:

| Region | Universal ID | Notes |
|--------|--------------|-------|
| AmpR gene | `coding` | β-lactamase, common selection marker |
| Origin of replication | `custom: plasmid_ori` | pMB1 ori |
| MCS | `custom: cloning_site` | designed restriction sites |

---

## Step 3 — Define feature signatures for each region

For every region in your inventory, write one file in
`organisms/<your_organism>/regions/<region_name>.json`:

```json
{
  "name": "athila_ltr_retroelement",
  "taxonomy_id": "dispersed_repeat",
  "description": "Athila family LTR retroelements, ~10 kb internal + ~1.5 kb LTRs flanking. Found in pericentromeric Arabidopsis heterochromatin.",
  "feature_signature": {
    "shannon_5mer": {"range": [2.5, 3.5]},
    "gc_content": {"range": [0.30, 0.45]},
    "kmer_multiplicity_p95": {"min": 50},
    "consensus_match": {"any_of": ["Athila", "AthGypsy"]}
  },
  "mapping_hints": {
    "expected_unique": false,
    "max_occ_multiplier": 50.0,
    "lambda_scale": 1.5,
    "report_multi_position": true
  },
  "sources": [
    "doi:10.1038/35048500",
    "doi:10.1126/science.286.5441.1023"
  ],
  "confidence": 0.85
}
```

`feature_signature` is a predicate evaluated by the runtime classifier. A
window is classified as this region iff **all** sub-predicates match. The
classifier ranks competing matches by confidence × specificity-of-match.

---

## Step 4 — (Optional) Run LLM-assisted curation

If you don't already have the feature ranges / consensus sequences in hand,
LLmap can curate them via:

```bash
llmap knowledge-fetch \
    --organism arabidopsis \
    --topic "Athila LTR retroelements: typical length, GC, consensus, copy count" \
    --output knowledge/organisms/arabidopsis/regions/athila.json
```

This issues a structured-output Claude API call (template at
`knowledge/prompt_templates/identify_region.md`), saves the curated entry,
and caches the API response in `~/.llmap/llm_cache/`.

Crucially: the fetch is a **one-time activity at module-build time**. After
that the JSON is committed to your module and mapping is fully deterministic.

For organisms with limited literature (e.g. exotic species), you can give
the LLM a fasta of representative sequences and ask it to derive features:

```bash
llmap knowledge-fetch \
    --organism mywild_lizard \
    --from-fasta sample_repeat_consensus.fa \
    --topic "Identify repeat structure, estimate copy count" \
    --output knowledge/organisms/mywild_lizard/regions/sample.json
```

For **fictive DNA** (test data with no biological meaning), simply hand-author
the JSON. The fictive module in `organisms/test_fictive/` is included as a
reference example.

---

## Step 5 — Build the classifier rules

Once your regions are defined, generate the deterministic classifier:

```bash
llmap classifier-compile \
    --regions knowledge/organisms/arabidopsis/regions \
    --output knowledge/organisms/arabidopsis/classifier_rules.json
```

This compiles all region feature_signatures into a single fast lookup
structure (decision tree). The classifier is now organism-aware,
deterministic, and dependency-free.

---

## Step 6 — Annotate your reference

Run the classifier over your reference genome:

```bash
llmap annotate-ref \
    --reference my_genome.fa \
    --organism arabidopsis \
    --output my_genome.regann
```

This produces a `.regann` companion file alongside your reference. It
contains one entry per 1 kb window with the assigned region type and the
feature vector.

---

## Step 7 — Use it for mapping

```bash
llmap align \
    --reference my_genome.fa \
    --region-annot my_genome.regann \
    --reads my_reads.fastq \
    --output mapped.bam
```

The chain DP now weights anchors per region, expects multi-position output
in `tandem_repeat` / `paralog_family` regions, and tightens scoring in
`coding` / `unique_single_copy` regions.

---

## Special cases

### Polyploids (wheat, strawberry, cancer aneuploids)

Polyploids have homologous chromosome copies with ~95–99% identity — the
same problem as paralog families but genome-wide. Set
`expected_unique: "ambiguous"` and `require_psv_disambiguation: true`
across the affected regions, and supply a `psv_catalog.tsv` listing the
homeolog-specific variants.

### Viruses

Tiny genome, often <30 kb. The 1 kb windowing in `annotate-ref` may be too
coarse — pass `--window-bp 100` to use 100 bp windows.

For RNA viruses with mutational hotspots (HIV, influenza), add `custom`
region entries with `mapping_hints.allow_high_mismatch: true` so the
identity threshold is relaxed in those regions only.

### Synthetic DNA

For oligo pools / designed sequences, you usually already know the answer:
the regions you care about are defined when you designed the construct.
Skip the LLM step entirely and hand-author `regions/*.json` directly. Map
with the resulting annotation file.

### Fictive DNA

Same as synthetic: hand-author. We ship `organisms/test_fictive/` as a
minimal reference module — it defines three made-up region types
(`zorg_repeat`, `klorf_coding`, `null_spacer`) and serves as a copy-paste
template.

### Mixed-organism or chimaeric samples (metagenomics, xenografts)

Build separate organism modules for each organism present, then run
`annotate-ref` multiple times against a concatenated reference. Each contig
gets the annotations of its source organism.

---

## What you should commit back to LLmap

If you build a high-quality module for a new organism, please consider
opening a PR. The structure under `knowledge/organisms/` is meant to be
community-extendable. Modules of reasonable coverage and well-cited
sources are valuable to other users of the same organism.

A module is "well-formed" when:
- every region entry has at least one `sources` reference
- `classifier_rules.json` compiles cleanly
- `annotate-ref` on a representative reference produces non-empty
  annotations covering >95% of the reference
- a smoke-test BAM produced with the module shows higher concordance
  with the gold-standard mapper than the same data mapped without
  `--region-annot`
