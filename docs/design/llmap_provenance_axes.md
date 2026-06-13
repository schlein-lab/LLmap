# LLmap Provenance / Contamination Mode — Design Spec

**Status:** Draft v0.1 — 2026-06-13
**Author intent (Schlein):** Most low-VAF "variants" in real data are NOT germline
— they are exogenous contamination, chemistry/damage artefacts, mapping
confusion, or cross-sample bleed. LLmap must not hide these: route each read into
a **separate provenance bucket** with **content flags**, so the contamination
spectrum is **quantifiable** per sample/gene/variant.
**Companion docs:** [llmap_mode_architecture.md](llmap_mode_architecture.md),
[segdup_catalog_spec.md](segdup_catalog_spec.md).
**Working notes:** `artifacts/2026-06-13_contamination_mode_design.md` (full taxonomy).

---

## 0. TL;DR

The unifying purpose: **kill the fake low-VAF variant.** Every source the operator
enumerated (EBV/LCL, microbiome/kitome, index-hopping, PhiX/spike-ins, 8-oxoG/
FFPE/deamination, RNA-editing, NUMTs, paralogs, processed pseudogenes, decoy/rDNA,
chimeras, CHIP/mosaicism, tumor-in-normal) masquerades as a low-frequency variant.

The mechanism is **native to LLmap's WaveCollapse core**: contamination detection =
**adding contaminant/artefact candidate buckets to the same EM.** A read already
collapses over host buckets (incl. paralog buckets); we extend the candidate set so
every read gets a **provenance posterior** over {host-germline, host-paralog, EBV,
microbiome, NUMT, pseudogene, chimeric, …}. A read preferring a contaminant bucket
**collapses there, tagged, never dropped** — lossless, flag-not-filter (same doctrine
as Transcript-Mode).

**Quantifiable falls out for free:** bucket occupancy = the sample's contamination
spectrum (EBV load, microbiome fraction, PhiX %, NUMT reads, paralog-ambiguous %).

---

## 1. The boundary — three classes, NOT one bucket model

Not everything is a *sequence-origin* bucket. Split by mechanism:

| Class | What it is | Realisation | Reuses |
|---|---|---|---|
| **A. Competing reference** (exogenous organism, NUMT, paralog, processed pseudogene, decoy/rDNA) | the read's true origin is a *different sequence* | **provenance bucket** (extra WaveCollapse candidate) | Mode-6 taxbin + contaminant panel; PSV + WaveCollapse |
| **B. Context-conditioned error** (8-oxoG, FFPE, end-deamination, A→I/C→U editing) | host sequence, a *base-level artefact* | **per-site evidence tag** (strand-bias + read-position + substitution profile) — NOT a bucket | `rnamod`, new substitution classifier |
| **C. Structural split** (ligation/PCR chimera, ONT fused/concatemer, RT template-switch) | one read = two molecules | **split-provenance** record | chimera detector (Block 7) |
| **D. Cross-sample / clonal** (index-hopping, CHIP/mosaicism, tumor-in-normal, LCL-acquired somatic) | needs cohort/paired context | **flag only**, downstream decides | — |

**Key invariant:** B and D are *not* origin → they get per-site/per-read evidence
tags, never a provenance bucket. Putting LCL-acquired somatic mutations or A→I edits
into a contaminant bucket would *strip real biology*. (Agent 1's sharpening.)

---

## 2. Layering — staggered always-on (cost-aware)

| Tier | Trigger | Contents | Cost |
|---|---|---|---|
| **always-on** (default; `--no-provenance` to disable) | every run | Class A mapping-byproducts (paralog/NUMT/pseudogene — LLmap computes these anyway) + Class B per-site evidence (strand-bias/substitution profile) | ~0 |
| **opt-in** (`--contaminant-panel FILE`) | panel supplied | Class A exogenous-organism bucket matches (EBV/Mycoplasma/kitome/PhiX/spike-in) vs reference panels | panel-mapping pass |
| **opt-in** (`--reconcile-cohort`) | cohort supplied | Class D cross-sample (index-hop, tumor-in-normal) | cohort join |

Output is **always** emitted (so the Σ-invariant is always checkable).

---

## 3. The probabilistic core (Agent 2 lane)

`ProvenanceClass` enum + per-read posterior. The WaveCollapse candidate set for a
read becomes `host_buckets ∪ contaminant_buckets`. After EM collapse:
- a read's **provenance** = argmax bucket-class with its posterior (the collapse
  confidence — *not* a hard call: the posterior is retained);
- a read stays **host** unless a contaminant bucket beats host by a **calibrated
  margin** (host-conservative threshold + per-class prior prevalence — else a
  real low-VAF somatic gets absorbed into an artefact bucket; cf. Agent 1's
  prior/normalisation caveat).

### Lossless Σ-invariant (extends `lossless_aggregator`)
Every input read lands in **exactly one** bucket (host-call OR one provenance
class). `Σ over classes = N_input`. The aggregator asserts this; a violation is a
dropped read = a bug. The variant caller sees only the clean host bucket; B/D
artefacts are **quantified separately**, never merged into germline.

### Output (Agent 2 lane)
- per-read BAM tags: `XB:Z:<class>` (provenance bucket), `XQ:f:<posterior>`,
  `XE:Z:<evidence>` (Class-B per-site flags).
- `<out>.provenance.bam` — flagged reads routed here, never the main BAM.
- **`<out>.contamination_spectrum.parquet`** — one row per class:
  `class, n_reads, fraction, mean_posterior, bases` → directly plottable
  ("0.8% PhiX, 1.2% index-hop, 2.1% NUMT-collapse").

---

## 4. Detector classes (Agent 1 lane)

- **Mapping-confusion** (NUMT/paralog/pseudogene): PSV + Transcript-Mode signals →
  Class-A buckets. NUMT = mtDNA-vs-nuclear PSV; pseudogene (GBAP1/GBA1, PMS2/PMS2CL,
  SMN1/2) = intron-presence + transcript PSV.
- **Damage/editing** (Class B): substitution-profile classifier over `rnamod` —
  8-oxoG (G→T strand-bias), FFPE (C→T/G→A), end-deamination (C→T read-end), A→I
  (A→G in Alu, ADAR), C→U (APOBEC). Per-site evidence tag, never a bucket.
- **Chimera** (Class C): Block 7 detector → split-provenance.

---

## 5. Implementation order

1. **`ProvenanceClass` enum + per-read provenance struct + Σ-invariant accounting** (Agent 2) — the spine; testable standalone.
2. **`contamination_spectrum` aggregator + Parquet/BAM-tag output** (Agent 2).
3. **Damage/editing per-site classifier** (Agent 1) — self-contained, no EM dependency.
4. **Mapping-confusion buckets** (Agent 1) — wires PSV/Transcript into Class-A.
5. **Contaminant-panel bucket pass** (`--contaminant-panel`) — opt-in, reuses taxbin.
6. **Cross-sample reconcile** (`--reconcile-cohort`) — Class D, last.

Steps 1-3 are no-regret and parallelisable (Agent 1 on 3, Agent 2 on 1-2).

---
## 6. Akzeptanz-Test (Operator-Direktive 2026-06-13)
Nach Fertigstellung des Modes: auf **1-2 Longread- + 1-2 Shortread-Genomen** laufen lassen
und die Provenienz-Bucket-Verteilung (`contamination_spectrum.parquet`) prüfen — welche
Special-Buckets real auftauchen (EBV/LCL bei 1000G-Linien? NUMT-Fraktion? PhiX? Paralog-Ambig?).
Erwartung als Sanity: Longread (HiFi/ONT) hat weniger Chimär/Damage als Shortread; LCL-Samples
(1000G/HapMap) zeigen EBV; low-input/FFPE zeigt dmg:*. Σ-Invariante muss auf jedem Genom halten.
Kandidaten-Daten auf the HPC cluster (zu bestätigen): HG002 PacBio-HiFi (longread), HG001-007 Illumina
deepvariant-BAMs (shortread), in `/beegfs/u/<user>/` bzw. `~<user>`.
