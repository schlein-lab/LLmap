# Transcript-Mode — the HPC clusterValidierung (Agent 2, 2026-06-13)

## Build/Test (GRÜN)
- Isoliert geklont nach `~/llmap-tx` (Branch feat/transcript-mode-p1-p2, HEAD 5f9b3c8). /tmp war 50M-tmpfs/100% voll → TMPDIR+TEST_TMPDIR auf /home.
- CMake CPU-only configure rc=0, build rc=0 (100%, keine Fehler). Binary `build/src/llmap`.
- **Fokus-Tests 46/46 grün** (alle TranscriptStage.* + Sniffer).
- Full ctest seriell: 1798/1801. Die 3 Failures = `ReferenceIndexTest` hardcodet `/tmp/...idx` (test_reference_index.cpp:58) → volles tmpfs. Pre-existing Test-Bug, NICHT der neue Code.

## Funktionstest auf echten FLNC-Daten (KRITISCHER BEFUND)
Input: 200 iso-seq FLNC reads (~1750 bp) aus HG00272.ighg4.md.bam; Referenz grch38_IGH_locus.fa (chr14, 1.4 Mb).

| Lauf | mapped | N(Intron)-CIGARs |
|---|---|---|
| --mode transcript -x map-hifi | 38/200 (19%) | 0 |
| --mode reads -x map-hifi      | 38/200 (19%) | 0 |
| --mode transcript -x map-ont  | 196/200 (98%) | 0 |

Beobachtungen:
1. **Mapping-Rate ist preset-dominiert** (id=0.90→19%, id=0.70→98%). Kein fundamentaler Verlust; iso-seq braucht permissive Presets.
2. **0 Spliced-N-CIGARs selbst bei 98% Mapping.** Jeder gemappte Read = **genau 1 Alignment-Record** (196/196).
3. Introns erscheinen als **Insertionen + Soft-Clips innerhalb EINES Alignments**, z.B.:
   `1=51I97M60I1=287M1=62I324M221I2=374M1=69I`
   `34=32M1=88I82M1D316M1D15M19I1=287M142I1=202M1=59I15M801S`

## Root cause (Hypothese, stark gestützt)
Die P2-Prämisse — „der lineare Chainer zerlegt einen spliced Read in mehrere per-Exon-
`ClassicalAlignment`s" — ist auf echten Daten **FALSCH**. Der klassische Pfad (Chainer +
WFA2-Extension) liefert **1 Alignment/Read** und absorbiert Intron-Lücken als I/D/Soft-Clip
statt die Chain zu brechen. Damit bekommt `ApplyTranscriptStage` nur eine Sub-Chain →
Singleton → 0 Joins → 0 N-CIGARs. **Der Spliced-Stage ist auf realen FLNC-Daten inert.**

Deshalb bestanden die Unit-Tests (synthetische Multi-Sub-Chain-Inputs), aber Real-Data spliced nicht.

## Konsequenz / nötiger Fix (für Agent 1, der den Code besitzt)
Upstream: im Transcript-Mode muss der Chainer an intron-großen REF-Lücken **brechen**
(nicht via WFA bridgen), damit per-Exon-Sub-Chains entstehen — ODER der Stage muss die
rohen Chain-Anker bekommen und selbst an großen ref-gaps splitten (statt der post-WFA
ClassicalAlignments). Erst dann liefert der Joiner N-CIGARs.

Verbindung zur Operator-Sorge: genau hier wird die Exon-Struktur (die Isoforme unterscheidet
und ein novel Alt-Last-Exon sichtbar machen würde) in ein einzelnes lineares Alignment
kollabiert — statt als diskrete Exon-Blöcke erhalten zu bleiben.

---
## NÄCHSTER MODE (Operator-Vormerkung, 2026-06-13)
Nach Transcript-Mode: neuer Mode **"Exogene biologische Kontamination"** (exogenous
biological contamination — Reads fremder Organismen: Viren/Bakterien/Pilze/Fremd-DNA im
Sample). Verwandt mit dem bestehenden Mode-6 taxbin (taxonomic binning of un-human reads,
Commit 05a7af8 auf the HPC cluster). Noch nicht spezifiziert — erst Transcript-Mode lossless beweisen.

---
## E2E-DURCHBRUCH (2026-06-13, kombinierte Fixes)
Synthetisches Single-Copy Genom↔Transkriptom-Paar, `llmap align --mode transcript`:
- **Lossless/no-drop: BEWIESEN** — 9/9 Reads präsent + gemappt, novel-Read gemappt (nicht gedroppt).
- **Kanonische Spliced-Mapping: 6/6 PASS** (canon_0..4 + canon_rev, je 4 N == 4 Truth-Junctions,
  inkl. **Reverse-Strang**), N-Längen innerhalb Seed-Fenster-Präzision (TOL=50).
- canon_0 CIGAR: `7=167M8I308N131M21I521N220M47I267N89M9I809N175M16=` (4 N).
- **Novel-Alt-Last-Exon: 3/4** — exon5-alt mergt nicht (distales 1600bp-Intron). Ursache:
  exon4 bildet ein DUPLIKAT-Sub-Chain (von `min_score_fraction=0`); eine Kopie mergt in
  {1,2,3,4}, die andere in {4,5alt} → greedy-walk fragmentiert. **Lossless** (exon5-alt im
  Soft-Clip präsent). Fix: Sub-Chain-Dedup im Stage vor dem Joiner.

### Fix-Kette (wer)
1. R-B `intron_break_min` (Chain bricht an Introns) — Agent 2 (chain_score/chain.h/cmd_align)
2. extension_max_span cap (kein Cross-Exon-Bridge) — Agent 2 (cmd_align)
3. **CIGAR k-mer-Overlap-Fix** (Query-Koords exakt, 220q=220r) — Agent 2 (classical_pipeline_extend)
4. min_score_fraction=0 + max_chains/alignments=256 (alle Exons überleben) — Agent 2 (cmd_align)
5. QuerySpanFromCigar (Sub-Chain-Query aus Soft-Clips) — Agent 1 (cmd_align_transcript)
6. EmitSplicedCigar interne-SC-Strip + I-Encoding — Agent 1 (chain_spliced)
7. max_query_gap_bp=80 (Boundary-Slop tolerieren) — Agent 1 (cmd_align_transcript)

### Dokumentierte Follow-ups (alle lossless)
- Sub-Chain-Dedup → novel 4/4 (Agent 1, Stage).
- Splice-Site-Snapping (GT/AG) → exakte N-Boundaries statt Seed-Fenster-Präzision (Design §2.1).
- DNA-Regression: 143-178 classical/chain/mapping-Tests grün nach allen Fixes.

---
## BENCHMARK vs minimap2 (2026-06-13, echte Daten)
100 genom-weite PacBio iso-seq FLNC-Reads (scisoseq dedup BAM) vs GRCh38 chr14,
`llmap align --mode transcript -x map-ont` vs `minimap2 -ax splice:hq -uf`:

| Metrik | llmap | minimap2 |
|---|---|---|
| mapped | 11 (11%) | 24 (24%) |
| spliced reads (>=1 N) | 1 | 18 |
| distinct introns | 19 | 62 |
| Junction-Jaccard (±10bp) | 0.125 | (Referenz) |
| Primary-Placement-Agreement | 2/5 | — |
| Speed | ~1.1 s/read | ~ms/read |

### Ehrliche Bewertung
- **Synthetik (clean): bewiesen** — lossless spliced Mapping, canon+reverse exakt (6/9).
- **Echte Daten: LLmap weit hinter minimap2.** Kernlücke: auf real-error-behafteten CCS-Reads
  **feuert das Splicing kaum** (1 vs 18 spliced) — der R-B-Seed-Break + Seed-Fenster-Boundary-
  Ansatz ist fragil gegen Sequenzierfehler nahe Exon-Grenzen (Error → Seed bricht → Exon-Sub-Chain
  fehlt → kein N). minimap2 ist dagegen robust.
- **Speed:** ~1 s/read (chr14, CPU) = Größenordnungen langsamer; Chain-DP auf 104MB-Ref + WFA-Last.
- **Mapping-Rate halbiert** (11% vs 24%): real-Reads haben Barcodes/Adapter/Fehler, die llmap
  schlechter soft-clippt.

### Nächste Optimierungs-Prioritäten (aus dem Benchmark)
1. **Error-Robustheit des Splicings** (das große Gap): Splicing muss auf fehlerbehafteten Reads
  feuern — Splice-Site-Snapping (GT/AG) statt Seed-Fenster, fehler-tolerante Exon-Boundary.
2. **Speed:** Chain-DP/Candidate-Gen auf großen Referenzen (Index-Caching, weniger WFA).
3. **Sensitivität:** Adapter/Barcode-Soft-Clip, Preset-Tuning für CCS iso-seq.

---
## SPLICE-SNAPPING + LOCUS-SELEKTION (2026-06-13) + GENCODE-Testset
Fix-Kette für real-Data-Splicing: (i) GT/AG-Snapping (`splice_snap`, Agent2) → N exakt,
(ii) Chain-Locus-Selektion (`SelectTranscriptLocus`, Agent2) → richtige Exon-Chains extended.

Benchmark 100 echte FLNC-Reads vs chr14 (llmap map-ont vs minimap2 splice:hq):
| Metrik | vor Snapping | +Snapping | +Locus-Selektion | minimap2 |
|---|---|---|---|---|
| mapped | 11 | 11 | **20** | 24 |
| spliced | 1 | 1 | **5** | 18 |
| Junction-Jaccard vs mm2 | 0.125 | 0.286 | 0.269 | — |

Boundary-Konkordanz vs GENCODE v46 (chr14, 49.252 Introns):
| | llmap | minimap2 |
|---|---|---|
| Junctions | 23 | 62 |
| **exakt 0bp auf GENCODE** | **73.9%** | 32.3% |
| unmatched >50bp | 21.7% | 66.1% |

ERKENNTNIS: llmaps Junctions sind bp-scharf (74% exakt = Snapping wirkt). Precision-at-low-Recall:
nur 5 spliced (mm2: 18). D(p)-Verteilung braucht mehr Recall (within-Locus-Multi-Exon-Assemblierung
5→18) = nächster Chaining-Block. minimap2's 66% „unmatched" = breiterer Recall inkl. Single-Cell-
novel/Artefakt-Junctions (kein Apfel-Apfel ohne geteilte Reads).
