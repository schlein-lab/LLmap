# Bucket-Driven Assemble-then-Map (lossless)

Status: design, converged with Operator (2026-06-13). Awaiting final confirm of
line **A** (transient probe) before the `cluster_consensus` V1 rework + pipeline
wiring. No smoothing is wired anywhere yet.

## Problem (measured, not assumed)

On the chr14 iso-seq benchmark (3000 reads):
- **Extension-bound:** WFA extension = 99.8% of align time; `--classical-only` ==
  full → the EM is *not* the bottleneck (the GPU-EM idea was wrong).
- **Chain-soup:** ~80–270 chains/read; ~17 WFA extensions/read, ~94% identity-
  rejected.
- **Recall:** SPLICE_RECALL 0.134, mapping 257/952 vs minimap2.

Speed and recall share ONE root: too many candidate placements explored per read,
and the genuine exons drown in the soup. The min_chain_score-cull fix (committed)
recovered short exons but the structural gap remains.

## The loop (Operator)

```
cluster ──▶ virtual probe construct ──▶ bucket (region prior) ──▶ map probe to the
   probable region ──▶ fit? ──┬─ yes ─▶ propagate placement to members (lossless)
                              └─ no  ─▶ re-bucket (next probable region)
                                         └─ buckets exhausted ─▶ flag genuine-ambiguous
```

A read is never dropped: a misfit is *informative* (the bucket hypothesis was
wrong) and triggers escalation, not deletion.

## Losslessness invariant (the Operator constraint)

- **Stored output = `members[]`:** every read keeps its own bases, its offset in the
  layout, its strand, and the per-column depth. Variable columns stay **multi-allele**
  — never collapsed to a majority base. 100% of the original sequence, depth, and
  sequencing direction are retained.
- **The probe construct is TRANSIENT (line A):** a throwaway search key (it may carry
  one representative base per column) used *only* to locate the region. It is never
  stored, never emitted, never the truth. "No smoothing" governs the *stored sequence*,
  not the disposable search key. (Line B alternative: map the layout's coordinate axis
  instead of a probe — stricter, no ephemeral smoothing; switch to B on request.)

## Components & roles (already built, green)

| Component | Role in the loop |
|---|---|
| `minimizer_cluster` | the cluster: strand-aware (canonical minimizers), self-frequency cap = pre-map repeat proxy, signed union-find (strand tracked, never dropped) |
| `cluster_consensus` (rework) | `members[]/offset/strand/depth` = lossless layout (OUTPUT); `MajorityBase` demoted to probe-helper (TRANSIENT) |
| provenance buckets + `M(pos)` | the **region prior**: bucket says WHERE to map (numt→NUMT loci, rdna→rDNA array, host/M≈1→the sharp locus) → narrows search, kills the whole-genome extension soup |
| `sufficiency` / `consensus_resolution` | the **fit-check**: does the probe collapse to a confident placement (single-dominant / last-K-levels agree)? |
| cascade escalation | misfit → next bucket; exhausted → honest genuine-ambiguous flag |

## Merge-gate (cluster correctness — the worst error is a chimeric/paralog-bridged contig)

- **self-frequency-adaptive min-overlap:** longer overlap required in repeat/low-unique
  regions (pre-map proxy; `M(pos)` refines post-map).
- **read-length-relative + absolute floor:** ≥ X% of the shorter read AND ≥ absolute
  floor (~300 bp long-read; scaled down for short reads where the ambiguity guard
  carries more weight).
- **identity floor in the overlap:** ~90% long-read / ~95% short-read — tolerant of
  *true allele differences* (two reads from different alleles legitimately differ;
  paralog separation is carried by overlap *length* + the ambiguity guard, not identity).
- **WaveCollapse ambiguity guard:** reads overlap but the placement is blurry → do NOT
  merge.
- **revcomp merged with strand tracked:** cDNA/genomic reads come from both strands of
  one locus → detect overlap in both orientations, store strand per read, merge across
  it (track, never drop).
- **transitive consistency:** A–B, B–C ⇒ verify A–C before transitively merging (else
  chimera risk).

## Why it fixes both axes

- **Recall (partial reads):** reads are partial (5'/3' truncation, degradation) — NOT
  whole transcripts (the earlier oversimplification). The construct spans the **union**
  of member coverage → longer than any single read → anchors uniquely → dissolves the
  91 bp-fragment trap.
- **Speed:** map ~M constructs instead of ~N reads (M ≪ N), each to its bucket's
  probable region (not whole-genome) → the 99.8% extension time collapses.

## V1 open items

- bucket → region-prior table (provenance class → candidate loci).
- re-bucket policy: order of buckets to try, exhaustion criterion.
- probe construction: representative-per-column vs longest member.
- propagation: construct placement → per-member offset → per-read CIGAR (each read's
  own bases preserved).
