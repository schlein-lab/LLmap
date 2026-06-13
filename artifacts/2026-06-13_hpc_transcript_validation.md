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
