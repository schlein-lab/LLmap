# Y-Palindrome and Acrocentric Short-Arm T1 Entries

**Status:** initial release — 10 T1 entries added 2026-06-02
**Curator:** schlein-lab
**Schema:** v0.2 (curated.schema.json)

## 1. What's in this batch

Ten new T1 (curated) catalog entries fill two structural-architecture gaps left after the IGHG4 and HLA/FCGR initial batches:

| Locus_id | Architecture | Haplotype class | Anchor reference |
|---|---|---|---|
| `Y_palindrome_P1` | `palindromic_inverted` | `Y_palindrome_P1_canonical` | Skaletsky 2003; Rhie 2023 |
| `Y_palindrome_P2` | `palindromic_inverted` | `Y_palindrome_P2_canonical` | Skaletsky 2003; Rhie 2023 |
| `Y_palindrome_P3` | `palindromic_inverted` | `Y_palindrome_P3_canonical` | Skaletsky 2003; Rhie 2023 |
| `Y_palindrome_P4` | `palindromic_inverted` | `Y_palindrome_P4_canonical` | Skaletsky 2003; Rhie 2023 |
| `Y_palindrome_P5` | `palindromic_inverted` | `Y_palindrome_P5_canonical` | Skaletsky 2003; Repping 2002; Rhie 2023 |
| `chr13_pArm_rDNA_cluster` | `complex_mosaic` | `acrocentric_rDNA_cluster_chr13` | Nurk 2022; Guarracino 2023 |
| `chr14_pArm_rDNA_cluster` | `complex_mosaic` | `acrocentric_rDNA_cluster_chr14` | Nurk 2022; Guarracino 2023 |
| `chr15_pArm_rDNA_cluster` | `complex_mosaic` | `acrocentric_rDNA_cluster_chr15` | Nurk 2022; Guarracino 2023 |
| `chr21_pArm_rDNA_cluster` | `complex_mosaic` | `acrocentric_rDNA_cluster_chr21` | Nurk 2022; Guarracino 2023 |
| `chr22_pArm_rDNA_cluster` | `complex_mosaic` | `acrocentric_rDNA_cluster_chr22` | Nurk 2022; Guarracino 2023 |

All 10 pass `python3 catalog/validate.py`.

## 2. Y palindrome coordinate provenance

**Skaletsky 2003 Table 3** provides arm/spacer/identity for P1–P8 but not chrY coordinates (Skaletsky used contig-level NT_011875 / NT_011903 coordinates that do not lift cleanly into hg38 — the GRCh38-Y assembly contains documented errors in P1–P3).

**Authoritative coordinate source: `chm13v2.0Y_inverted_repeats_v1.bed`** from `s3://human-pangenomics/T2T/CHM13/assemblies/annotation/` (marbl/CHM13). Direct extraction:

| Palindrome | left-arm start | left-arm end | spacer | right-arm start | right-arm end | arm length (bp) | spacer (bp) | total span (bp) | identity % |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| P1 | 24,167,431 | 25,503,792 | 283 kb | 25,787,351 | 27,122,770 | 1,336,361 | 282,559 | 2,955,339 | 99.97 |
| P2 | 23,897,358 | 24,035,399 | 2 kb | 24,037,534 | 24,167,130 | 138,041 | 2,135 | 269,772 | 99.97 |
| P3 | 22,760,852 | 23,044,266 | 169 kb | 23,213,315 | 23,582,346 | 283,414 / 369,031* | 169,049 | 821,494 | 99.94 |
| P4 | 19,356,702 | 19,546,840 | 40 kb | 19,586,496 | 19,776,671 | 190,138 | 39,656 | 419,969 | 99.98 |
| P5 | 18,362,332 | 18,857,802 | 3 kb | 18,861,260 | 19,356,680 | 495,470 | 3,458 | 994,348 | 99.98 |

*P3 right arm is longer than left arm — asymmetric due to amplicon-mosaic insertion.

**GRCh38-Y coordinates: APPROXIMATE.** Rhie 2023 explicitly documents that P1–P3 in GRCh38-Y are discordant with T2T-Y due to assembly errors in GRCh38. P4–P8 are concordant in arm/spacer lengths and identity (99.84–99.96%). For the catalog GRCh38 coords are given as ±200 kb intervals, anchored to the AZFb/AZFc literature (Navarro-Costa 2010 PMC2910558: AZFa ~12.9–13.7 Mb, AZFb ~18.1–24.7 Mb, AZFc ~23–26.8 Mb).

## 3. Acrocentric short-arm coordinate provenance

**Nurk 2022 (T2T-CHM13 v2.0)** provides the first complete acrocentric p-arm sequences. rDNA array coordinates extracted from the paper / browser tracks:

| Chromosome | p-arm length (Mb) | rDNA array start | rDNA array end | rDNA size (Mb) |
|---|--:|--:|--:|--:|
| chr13 | 16.5 | 5,817,416 | 9,348,041 | 3.53 (largest) |
| chr14 | 16.2 | 2,099,537 | 2,817,811 | 0.72 (smallest) |
| chr15 | 17.1 | 2,506,442 | 4,707,485 | 2.20 |
| chr21 | 10.9 | 3,108,298 | 5,612,715 | 2.50 |
| chr22 | 15.1 | 4,793,794 | 5,720,650 | 0.93 |

**Guarracino 2023 (HPRC acrocentric pangenome)** quantifies the pseudo-homologous region (PHR) totals: chr13 4.53 Mb, chr14 6.48 Mb, chr15 0.72 Mb (outlier), chr21 3.79 Mb, chr22 2.81 Mb. SST1 array on 13p11.2/14p11.2/21p11.2 forms a >99% identity SD community (direct orientation on chr13/21, inverted on chr14). A 630-kb SD 1.6 Mb distal of the SST1 array is a candidate Robertsonian-translocation breakpoint for rob(13;21).

**GRCh38 coordinates are explicitly `null`** for all five acrocentric short arms. The hg38 acrocentric p-arms are N-masked (~30 Mb of Ns total across the five chromosomes); UCSC `genomicSuperDups` does not annotate them. This is the canonical example of `pangenome_only` clinical_function: the loci are real, clinically important (Robertsonian translocations cause ~1 in 1300 births and ~5% of Down syndrome cases), and absent from the standard reference.

## 4. Coordinates that were hard to nail down

Of the 10 loci, three required compromise:

1. **GRCh38 coords for Y_palindrome_P1, P2, P3.** Rhie 2023 explicitly says these are discordant in GRCh38-Y because of assembly errors. I wrote them as approximate ±200 kb intervals tied to AZFc literature, but anyone using these coords should prefer the T2T-CHM13 coords as authoritative. Listed as `note: "Approximate"` in the JSON.
2. **GRCh38 coords for chr13/14/15/21/22 p-arm.** Set to `null` — the hg38 p-arms are N-masked and no coordinate is defensible. Documented explicitly in each entry's `architecture_note`.
3. **HPRC_pangenome_v1.1 coords for the Y palindromes.** Left as `null` — per-haplotype Y assemblies vary substantially (the 43-sample Y paper, Hallast/Yang 2023, doi:10.1038/s41586-023-06425-6, shows extensive structural variation that would require per-haplotype entries). Population coords were not extractable.

## 5. Why these 10 and not others

- **Y palindromes P1–P5**: highest clinical relevance (AZFa/AZFb/AZFc spermatogenic-failure deletions). P6–P8 are biologically interesting but lack a recurrent clinical CNV product; deferred for a future batch.
- **Acrocentric chr13/14/15/21/22**: ClinGen lists Robertsonian translocations between all 10 chromosome pairs (13;14, 13;15, 13;21, 13;22, 14;15, 14;21, 14;22, 15;21, 15;22, 21;22) — the five chromosomes are the substrate for one of the most common balanced-translocation classes in clinical cytogenetics (~1 in 1300 live births). All five are absent from GRCh38 and are the canonical `pangenome_only` examples.

## 6. Cross-references

Y-palindromes cross-reference all four siblings via `related_loci`. Acrocentric p-arm entries cross-reference all four siblings. This is the first batch where intra-locus cross-references are used systematically — earlier IGHG4/IGH entries only had pairwise refs.

## 7. Mapping-strategy caveats

The acrocentric entries set `multi_position` as the **default stage-1 fallback** (not stage-3 as in IGHG4). This is honest: rDNA arrays + SST1 are not uniquely mappable at the read level with current short-read or even standard HiFi technology. Unique mapping requires either (a) reads ≥100 kb spanning SD context or (b) per-haplotype assembly. The mapping strategy reflects this — LLmap should report top-K positions rather than force a single answer.

The Y-palindrome entries use spacer/flanking anchors as fallback stage 1, because palindrome arms (99.94–99.98% identity, single discriminating bp every 1–2 kb) are tractable for unique mapping ONLY when reads span the spacer or the palindrome boundary.

## 8. References

- Skaletsky H et al., *Nature* 2003. doi:10.1038/nature01722
- Rozen S et al., *Nature* 2003. doi:10.1038/nature01723
- Repping S et al., *AJHG* 2002. doi:10.1086/342928
- Kuroda-Kawaguchi T et al., *Nat Genet* 2001. doi:10.1038/ng757
- Rhie A et al., *Nature* 2023 (T2T-Y). doi:10.1038/s41586-023-06457-y
- Hallast P, Yang C et al., *Nature* 2023 (43 Y assemblies). doi:10.1038/s41586-023-06425-6
- Nurk S et al., *Science* 2022 (T2T-CHM13 v2.0). doi:10.1126/science.abj6987
- Hoyt SJ et al., *Science* 2022 (T2T repeats). doi:10.1126/science.abk3112
- Guarracino A et al., *Nature* 2023 (acrocentric pangenome). doi:10.1038/s41586-023-05976-y
- Navarro-Costa P et al. 2010 (PMC2910558) — AZF region review with hg18/hg19 coords

## 9. Files added

```
catalog/curated/Y_palindrome_P1.json
catalog/curated/Y_palindrome_P2.json
catalog/curated/Y_palindrome_P3.json
catalog/curated/Y_palindrome_P4.json
catalog/curated/Y_palindrome_P5.json
catalog/curated/chr13_pArm_rDNA_cluster.json
catalog/curated/chr14_pArm_rDNA_cluster.json
catalog/curated/chr15_pArm_rDNA_cluster.json
catalog/curated/chr21_pArm_rDNA_cluster.json
catalog/curated/chr22_pArm_rDNA_cluster.json
docs/design/y_palindrome_acrocentric_t1_entries.md   # this file
```
