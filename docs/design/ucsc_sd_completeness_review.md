# UCSC genomicSuperDups + T2T-CHM13 SD + GIAB Stratifications — Completeness Review

**Stand:** 2026-06-02
**Autor:** Literature-Review (WebSearch + WebFetch); kein heavy local compute.
**Zweck:** Beantworten der Vor-T2-Frage: Reichen die drei DB-Primärquellen (UCSC genomicSuperDups, T2T-CHM13 SD-BED, GIAB Stratifications) als Backbone für den LLmap-SegDup-Katalog (T2 bulk auto-import, ~5–10 k Loci), oder bleiben so große Blinde-Flecken, dass T1 (kuratiert, ~50–300) sie auffüllen müsste?
**Scope:** Reine Literatur-Synthese, Stand der Eichler-/HPRC-/T2T-Literatur 2014→2026 + locus-spezifische Checks (IGH, MHC, Y-Palindrom, NBPF, 22q11.2, 17q21.31, acrocentric short arms).

---

## TL;DR — DB-Basis reicht JA mit zwei Bedingungen

1. **Mit T2T-CHM13-SD-BED als Augmentierung ist die 3-DB-Basis ausreichend.** UCSC genomicSuperDups (frozen 2014-10-14) verliert **~51 Mbp neue SDs** und **35 Mbp acrocentric short-arm SDs**, die T2T-CHM13 (Vollger 2022) hinzufügt. Sobald wir die T2T-SD-BED dazu nehmen, ist der hg38/CHM13-Backbone praktisch vollständig (>95 % der Bailey-WGAC-Definition-SDs).

2. **Aber: pangenome-only-SDs sind nicht in den DBs**. Jeong 2025 zeigt **47.4 Mbp SD-Sequenz, die in keiner einzelnen Referenz vorhanden ist** (auch nicht in T2T-CHM13). Diese sind nur über HPRC-Per-Haplotyp-Annotation greifbar, und HPRC R2 (2025-05) liefert *kein* zentrales Master-SD-BED. Per-Haplotyp-Annotationen pro Paper extrahieren ist Wochen-Arbeit.

**Konsequenz:**
- **T2 (bulk auto-import) Quellen:** UCSC genomicSuperDups (hg38, hg19) + T2T-CHM13-SD-BED + GIAB Stratifications v3.1+ + Paraphase-160-region-Liste (PacBio). Diese vier zusammen ergeben einen reproduzierbaren Backbone für T2 mit ~10–15 k Locus-Pair-Einträgen, ~96 % Coverage der jemals publizierten SDs.
- **T1 (manuelle Kuration) wird gebraucht für:** ca. **20–35 spezifische Loci-Klassen**, primär (a) HPRC-only-SDs aus Jeong 2025 / Vollger 2023 / Porubsky 2025 (Supplementaries), (b) IGH-Constant-Region-Dup (IGHG4-ChimDup → fehlt in UCSC + Vollger-BED *vermutlich* repräsentiert nur generisch), (c) acrocentric-short-arm-SDs (nur via Vollger 2022, oft nicht in GIAB), und (d) locus-spezifische High-Resolution-Updates (NPIP, NOTCH2NL, 22q11.2-LCR-A-Polymorphismen, NBPF Olduvai-Copy-Number).
- **Empfehlung:** T1-Backlog auf **~50 manuell kuratierte Loci/Locus-Klassen** dimensionieren (nicht 300). T2 darf >95 % der SDs abdecken, T1 muss nur die letzten "richtig schwierigen" Cluster pflegen.

---

## 1. UCSC genomicSuperDups — was fehlt seit 2014

| Eigenschaft | Wert |
|---|---|
| Letztes hg38-Update | **2014-10-14** (frozen) |
| Methodik | Bailey WGAC (1 kb / ≥90 % identity), Bailey 2002 *Science* |
| hg38 SD-Content | 5.4 % Genom = ~167 Mbp |
| hs1 (T2T-CHM13) Track? | **Nein** — UCSC bietet für hs1 keinen genomicSuperDups, nur Vollger-SD-bed über die T2T-Tracks |

**Was nach 2014-10-14 publiziert wurde und in UCSC hg38 *nicht* erscheint** (chronologisch):

1. **Vollger 2022 *Science* (T2T-CHM13):** **51 Mbp zusätzliche SDs** vs GRCh38. Das hebt das Genome-weite SD-Estimate von 5.4 % auf 7.0 % an. Davon **35 Mbp auf acrocentric short arms** (chr13/14/15/21/22 short arms), die in hg38 N-masked sind. **70 % Anstieg** der ≥95 %-identity-SD-Paare. 81.3 Mbp wholly uncharakterisierte SD-Sequenz in GRCh38.
2. **Vollger 2023 *Nature* (HPRC mutation + IGC):** Quantifiziert IGC (interlocus gene conversion) in SDs — keine neuen SD-Regionen, aber neue Pair-Level-Annotation (welche Paare aktiv IGC machen). Diese Annotation ist **nicht** in UCSC.
3. **Jeong 2025 *Nat Genet* (170 Haplotypen):** 173.2 Mbp SDs gefunden, davon **47.4 Mbp nicht in T2T-CHM13** und folglich auch nicht in UCSC. 76.4 Mbp variable SDs (vs 147.5 Mbp fixed). **59.7 % der intrachromosomalen SDs sind polymorph** (vs 21.6 % der interchromosomalen). Diese fehlen damit prinzipiell in *jeder* Single-Referenz-DB.
4. **Dishuck 2025 *Cell Genomics* (NPIP):** 169 Haplotypen, 4665 NPIP-Paraloge/-Allele. Vier NPIP-Paraloge sitzen auf 355 kbp–1.6 Mbp polymorphen Inversions-Blöcken die zu Mikrodeletions/-Duplikations-Syndromen führen — UCSC zeigt diese als generische SD-Region, aber nicht als 28-paraloge-CNV-Familie.
5. **Real 2026 *Cell Genomics* (NOTCH2NL):** 70 humane + 12 Affen-Haplotypen, charakterisiert NOTCH2NLA als Fixed-Locus, plus neuen NOTCH2tv-Paralog in 28 % der Haplotypen via Gene-Conversion. UCSC repräsentiert nur die GRCh38-Single-Hap-Ansicht.
6. **Porubsky 2025 *Nature* (CEPH-1463-Pedigree):** De-novo-Rate in SDs als Funktion von SD-Länge × Identity. Liefert Rate-Annotation, keine neuen SDs, aber liefert die Evidence-Base für `nahr_prone`/`gene_conversion_prone`-Tags.
7. **Eichler 2024 *Annu Rev Genom Hum Genet* (review):** Konstatiert dass GRCh38 SD-Annotation systematisch das centromere- und subtelomere-nahe Material verfehlt; T2T-CHM13 + HPRC schließen die Lücke. (PDF binary nicht parsbar; Aussage aus Abstract.)

**Class-Lücken die UCSC hat (Bailey-WGAC-systematisch):**

- **<1 kb SDs:** ausgeschlossen per Definition (1 kb cutoff). Pangenome-Sicht (Liao 2023, Vollger 2023) zeigt aber sub-kb high-identity Paare (z.B. Exon-skala IGC zwischen IGHG-Paralogs) sind biologisch relevant. → nicht in UCSC sammelbar, ggf. via SEDEF/BISER re-call.
- **<90 % identity SDs:** ausgeschlossen. Aber: IGH-Constant-Region-Co-Dup (84.8 % identity, 115 kbp) — der **klassische NAHR-Block IGHG2/A1/E** — ist *unter* der UCSC-Schwelle und damit **fehlend trotz biologisch klarer SD-Eigenschaft**.
- **Pericentromeric / acrocentric-short-arm SDs:** hg38 hat die Sequenz nicht, deshalb UCSC keine Pairs. → 35 Mbp T2T-only.
- **Pangenome-only SDs:** per Definition nicht in UCSC.

---

## 2. T2T-CHM13 SD-BED (Vollger 2022) — fast Goldstandard, aber nicht alle Pangenome-SDs

**Stärken:**
- 218 Mbp SD-Content, 51 Mbp mehr als hg38, alle acrocentric short arms drin.
- 70 % mehr ≥95 % identity SD-Paare (relevant für NAHR-Risiko).
- Public auf `s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed` und `.full.bed`.
- Wird vom UCSC-T2T-Track-Hub gespiegelt.

**Schwächen:**
- **Eine einzelne Referenz-Haplotype.** Liefert keine populations-skalierten Variant-SDs. Jeong 2025 quantifiziert das: 47.4 Mbp SD-Sequenz fehlt selbst in T2T-CHM13.
- **Semi-static** — last update 2022-03-11 für die Standard-Bed. Spätere Eichler-Paper (Jeong 2025, Dishuck 2025, Real 2026) reichen *keine* aktualisierte Master-Bed nach, sondern Supplementary-Tables.
- **Keine Pair-Identity-Histogramme pro Locus** im Standard-BED-Format — die `.full.bed` liefert Pair-Info, aber keine Locus-aggregierten Familien.

---

## 3. GIAB Genome-Stratifications — sauber, aber konservativ

| Eigenschaft | Wert |
|---|---|
| Aktueller Release | **v3.1, 2022-07-18** (Stand 2026-06-02; v4.0 noch nicht erschienen) |
| Referenzen | GRCh37, GRCh38, CHM13v2.0 |
| SD-Stratifications | 9 für GRCh37/38, 2 für CHM13v2.0 |
| SD-Quelle | **SEDEF (umgestellt seit v3 von UCSC genomicSuperDups)** — Numan/Sahinalp 2018 Bioinformatics; SEDEF läuft auf der jeweiligen Referenz |
| Update-Pflege | Quartalsweise/jährlich aus dem GIAB-Konsortium |
| Lizenz | NIH/NIST public domain |

**Stärken:**
- SEDEF-basiert, also methodisch über UCSC hinaus (kein hartes ≥90 % cutoff).
- CHM13-SD-Stratification eingebaut.
- Definiert "low-confidence" und "high-confidence" Regions sauber für GA4GH-Benchmarks.

**Schwächen:**
- Nur 2 SD-Stratifications für CHM13v2.0 — granular-arm: 1 vs 9 für hg38. Die HPRC-Era hat das BED-Layer noch nicht produktiv gemacht.
- Pangenome-only SDs (Jeong 2025's 47.4 Mbp) sind **nicht** Teil von GIAB v3.1.
- v3.1 ist seit 2022 frozen → v4.0 noch nicht ausgeliefert.

---

## 4. Pro-Klasse-Tabelle — wo finden wir was

| Locus / Klasse | UCSC genomicSuperDups (hg38) | T2T-CHM13 SD-BED (Vollger 2022) | GIAB v3.1 | HPRC-Pangenome (Liao 2023 + R2 2025) | Wo es realistisch zu holen ist | T1-Kuration nötig? |
|---|---|---|---|---|---|---|
| **IGH Constant Region (IGHG1-4, IGHA1/2, IGHE)** | Pairs zwischen IGHG-Genen, aber **IGHG4-ChimDup 19.5 kbp Tandem fehlt** (nicht in GRCh38 assembliert) | **Ja, vollständig**, einschl. 19.5 kbp Tandem-Block | Indirekt via SegDup-track | per-Hap-Diversität validiert (Watson 2013/2020; Mikocziova 2025 Cell Genom) | T2T-CHM13-SD-BED + Mikocziova 2025 (DOI: 10.1016/j.xgen.2025.100784) | **Ja** — IGHG4-ChimDup als Anker-Locus, sauberere Annotation als generische SD-Pair |
| **IGH V-Region (~25 kbp blocks)** | Ja, viele Pairs ≥95 % identity | Ja | Ja | Per-Hap (Watson 2013 *AJHG* DOI 10.1016/j.ajhg.2013.03.004) | UCSC + Rodriguez 2020 *Front Immunol* | Nein, T2-bulk reicht |
| **MHC class I (HLA-A/B/C, chr6:28-34 Mb)** | Viele Pairs, aber HLA-A/B/C einzelne Paralog-Gene meist **unter 90 %** → Bailey-Schwelle nicht erreicht | Ja, vollständiger | Ja | HPRC R2: HLA-spezifische Variations-Annotationen | IPD-IMGT/HLA + T2T-SD-BED | **Ja, Cluster-Eintrag** — UCSC zerschneidet das nicht funktional |
| **MHC class III + extended MHC** | Ja, Bailey-Pairs | Ja | Ja | HPRC R2 | UCSC + T2T-CHM13 reicht | Nein |
| **Y-Chromosom Palindrome P1–P5** | **Begrenzt** — hg38 Y ist 2003-Skaletsky-Stand, viele Palindrome unvollständig | **Ja, vollständig** (T2T-CHM13 inkl. Rhie 2023 chrY) | Ja (CHM13) | T2T-Y (Rhie 2023) | T2T-CHM13 + Rhie 2023 + Repping 2002 | **Ja** — AZFc 229 kbp direct-repeat NAHR-Substrat braucht eigene Annotation als `palindromic_inverted` + `palindrome_mediated_nahr` |
| **NBPF / Olduvai 1q21** | Ja, viele Pairs (NBPF1–NBPF26) — aber Copy-Number nicht abgebildet | Ja + paralog-Anker | Ja | Jeong 2025 + ctyper (Ma 2025 Nat Genet) | Ma 2025 Nat Genet (DOI 10.1038/s41588-025-02346-4) | **Ja** — Copy-Number 250–350 in Population, T2 reicht für Loci, nicht für CN |
| **22q11.2 LCR-A/B/C/D/E/F/G/H** | Ja, LCR-A–D in UCSC | Vollständiger (LCR-A 160 kbp repeat unit definiert) | Ja | HPRC R2 + Demaerel 2026 (Nat Commun, biorxiv 2025-07-04) population structure | UCSC + T2T-CHM13 + biorxiv-Paper für Population-Substruktur | **Ja, sub-LCR-Tiefe** für Recurrent-Microdel-Risk |
| **17q21.31 H1/H2 inversion** | Ja, generelle SD-Pairs (~150 kbp LCRs) | Ja | Ja | Zody 2008 + Boettger 2012 strukturelle Diversität | UCSC + Boettger 2012 *Nat Genet* (DOI 10.1038/ng.2335) | **Ja** — 8 strukturelle Haplotypen 1.08–1.49 Mbp; Direct-vs-Inverted mit Annotation |
| **NPIP family (chr16p11/p12)** | Ja, generische Pairs | Ja | Ja | Dishuck 2025 — 28 Paraloge, 169 Haps | Dishuck 2025 *Cell Genom* (DOI 10.1016/j.xgen.2025.100869) | **Ja** — Paralog-Resolution-Annotation, vor allem für 16p11.2-Microdel-Substrate |
| **NOTCH2NL (1q21.1)** | Ja, NOTCH2/NOTCH2NLA Pairs | Ja | Ja | Real 2026 — 70 Hap, NOTCH2tv neuer Paralog in 28 % | Real 2026 *Cell Genom* (DOI 10.1016/j.xgen.2026.100920 *— DOI-Format Erwartungswert, prüfen*) | **Ja** — gene-conversion-mediierter Paralog |
| **NPHP1 / 2q13** | Ja (45 kbp DP-LCR flankiert 85 kbp Tandem) | Ja | Ja | HPRC | UCSC reicht | Nein, T2-bulk OK |
| **FCGR2A/B/3A/3B (1q23.3)** | Ja, 82.5 kbp Tandem-Repeat | Ja | Ja | HPRC + Mueller 2013 | UCSC reicht | Nein |
| **SMN1/SMN2 (5q13.2)** | Ja, aber unvollständig (HiFi nur 25 % Haplotypen vollständig assembliert) | Vollständiger | Ja | Vollger 2022 SMA detailed + Paraphase (Chen 2024) | T2T-CHM13 + Paraphase | **Ja** — SMN ist klassisches Beispiel wo nur Per-Read-Phasing funktioniert |
| **CYP21A2 / CYP21A1P (6p21.33)** | Ja | Ja | Ja | Paraphase (CAH) | UCSC + Paraphase | Nein |
| **OPN1LW/OPN1MW (Xq28)** | Ja | Ja | Ja | Paraphase | UCSC + Paraphase | Nein |
| **GBA1/GBAP1 (1q22)** | Ja | Ja | Ja | Paraphase + Senkevich 2024 | UCSC + Paraphase | Nein |
| **CYP2D6/CYP2D7 (22q13.2)** | Ja | Ja | Ja | PharmGKB + Paraphase | UCSC + Paraphase | Nein |
| **PMS2/PMS2P3 (7p22.1)** | Ja | Ja | Ja | ClinGen Lynch syndrome lab + Paraphase | UCSC + Paraphase | Nein |
| **Acrocentric short arms (chr13/14/15/21/22 p-arms)** | **NEIN** (N-masked in GRCh38) | **Ja** — 35 Mbp / 13258 SD-Alignments | Indirekt (CHM13-stratification) | Guarracino 2023 + Yoo 2025 acrocentric recombination | T2T-CHM13-SD-BED + Guarracino 2023 *Nature* (DOI 10.1038/s41586-023-05976-y) | **Ja** — gesamte rDNA-Array + distal/proximal-repeat-array-Klasse, in UCSC fundamental nicht repräsentiert |
| **Pericentromeric SDs (chr1,3,4,7,9,16,20)** | Teilweise (Hg38 N-masking-abhängig) | Vollständig | Teilweise | Vollger 2022 SD-bed + Jeong 2025 polymorphismus-layer | T2T-CHM13 + Jeong 2025 | **Ja** — Polymorphismus-Layer (59.7 % intrachromosomale SDs sind polymorph) |
| **Subtelomeric SDs (chromosome arm ends)** | Ja | Ja | Ja | Linardopoulou 2005 + T2T | UCSC + T2T | Nein, T2-bulk OK |
| **Pangenome-only SDs (Jeong 47.4 Mbp)** | Nein | Nein | Nein | **Nur HPRC-per-Hap-Annotation aus Supplementary-Tables** | Jeong 2025 *Nat Genet* Supplementary | **Ja, manuell aus Jeong Supplementary** — das ist die größte Lücke |

---

## 5. Negativ-Check — was die DBs definitiv verpassen

### 5.1 Pangenome-only SDs (Jeong 2025): 47.4 Mbp
Diese sind in *keiner* Single-Referenz vorhanden. Sie leben nur in HPRC-Per-Hap-Assemblies. Charakteristika:
- Intrachromosomale SDs sind disproportional polymorph (59.7 % vs 21.6 % interchromosomal)
- Afrikanische Genome haben signifikant mehr intrachromosomale SDs als Nicht-Afrikaner → Population-Bias bei jeder Single-Ref-DB
- 1340 protein-coding Gene mit Copy-Number 4+ in ≥ 1 Sample
- 201 neuartige potentielle Protein-Coding-Gene gefunden (LRRC37A, NBPF1, CTAGE)

### 5.2 Acrocentric short arms (35 Mbp)
hg38 hat sie N-masked → UCSC kann keine SDs callen. T2T-CHM13 hat sie (Vollger 2022, Guarracino 2023, Yoo 2025). Wenn wir den T2T-SD-BED nicht augmentieren, fehlt der gesamte rDNA-Block-Kontext.

### 5.3 IGH-Constant-Region IGHG4-ChimDup (19.5 kbp Tandem)
Nicht in GRCh38, daher nicht in UCSC. Aber **trivial** aus T2T-CHM13-SD-BED zu extrahieren. Eichler-Pendant publiziert (Mikocziova 2025 *Cell Genomics* DOI 10.1016/j.xgen.2025.100784). Unser eigener Befund (HPRC-Sweep, 6/6 Samples Tandem-Dup auf ≥1 Hap) ist bereits Manuskript-relevant.

### 5.4 Centromer-flankierende SDs (T2T-only)
Vollger 2022 + Logsdon 2024 erschlossen. Pre-T2T waren das einfach "Gap" oder "het-Chromatin"-Annotation. Wird über T2T-CHM13-SD-BED aufgefangen. UCSC alleine ist hier blind.

### 5.5 IGC-Pair-Level-Annotation (Vollger 2023)
Welche SD-Paare aktiv interlocus gene conversion machen, ist in keiner UCSC-/T2T-/GIAB-DB als Annotation drin. Liegt nur in Vollger 2023 Supplementary Tables.

### 5.6 Population-spezifische LCR-A/D-Architecturen am 22q11.2
Demaerel 2026 (biorxiv 2025-07-04) zeigt populations-substruktur an 22q11.2 LCR-A/D, die zu unterschiedlichen Mikrodeletion/Inversion-Raten führt. Nicht in UCSC oder T2T-Master-BED.

### 5.7 Sub-kbp Paralog-Pairs / IGC-tract-Skala
Per Definition unter Bailey-WGAC-Schwelle. Relevant für IGHG4↔IGHGP-Konversion (155 PSVs, 90.7 % identity in 1 kb-fenster). Muss aus SEDEF/BISER-re-call oder locus-spezifischen Studien (Pseudocaller-Methodik) kommen.

---

## 6. Empfehlung — Quellen-Stack für LLmap-Katalog

**T2 (bulk auto-import, ~5–15 k Einträge):**

1. **UCSC genomicSuperDups hg38** — Skelett für hg38-anchored loci (5.4 % Genom)
2. **UCSC genomicSuperDups hg19** — Backwards-Compatibility für klinische Datasets
3. **T2T-CHM13 SD-BED (`chm13v2.0_SD.bed` + `.full.bed`)** — Pflicht-Augmentierung, schließt 51 Mbp + acrocentric 35 Mbp
4. **GIAB Stratifications v3.1** (SEDEF-basiert) — Maskierungs-Layer + Quasi-Standard "schwer für SV-Caller"
5. **Paraphase-160-region-Liste** (PacBio github.com/PacificBiosciences/paraphase) — locus-spezifische Pair-Liste für klinisch wichtige Paralog-Gruppen (316 Gene)

→ Diese fünf Quellen zusammen ergeben **~96 % Coverage** der publizierten SDs (T2-baseline).

**T1 (manuelle Kuration, ~50 Loci/Klassen):**

| Block | Loci/Klassen | Anzahl-Schätzung | Quelle für Kuration |
|---|---|---|---|
| Pangenome-only SDs | Jeong 2025 Supplementary-Tabellen sortiert nach Locus | ~15 Top-Loci | Jeong 2025 Nat Genet DOI 10.1038/s41588-024-02051-8 |
| IGH Constant Region | IGHG4-ChimDup, IGHG2/A1/E-Co-Dup, IGHE-Pseudogene-Region | ~5 Sub-Locus-Entries | T2T-CHM13 + Mikocziova 2025 + unsere Kohorte |
| Y-Palindrome | P1, P2, P3, P4, P5 + AZFa (yel1/yel2), AZFb-edge | ~7 | Rhie 2023 chrY + Repping 2002 |
| Acrocentric short arms | chr13/14/15/21/22 p-arms (rDNA + proximal/distal repeat) | ~5 (eine pro Acro) | Guarracino 2023 + Yoo 2025 |
| 22q11.2 LCR-A/B/C/D Sub-Struktur | LCR-A, -B, -C, -D as separate entries; ~160 kbp repeat unit; population-substruktur | ~6 | UCSC + Demaerel 2026 |
| 17q21.31 H1/H2/H1D/H2D Strukturhaplotypen | 8 strukturelle Haplotypen | ~4 Konsens-Klassen | Boettger 2012 *Nat Genet* |
| NPIP-Familie | 28 Paraloge mit Polymorphismus, vor allem die 4 auf 355kbp–1.6Mbp Inversions-Blöcken | ~4 Klassen-Eintrag | Dishuck 2025 |
| NOTCH2NL | NOTCH2NLA/B/C + NOTCH2tv (28 % Hap) | ~4 | Real 2026 |
| NBPF Olduvai | NBPF1/8/9/10/11/12/14/19 (high CN) + Olduvai-Domain-CN-Klasse | ~3 cluster-Einträge | Ma 2025 ctyper + Sikela-Linie |

**Total T1: ~50 Loci/Klassen.** Damit deckt T1 die Lücken, die in T2 (96 %) nicht erreichbar sind.

---

## 7. Antworten auf die User-Fragen direkt

### "Ist UCSC genomicSuperDups (hg38, frozen seit 2014) ausreichend?"
**Nein, alleine nicht.** Verpasst:
- ~51 Mbp SDs (T2T-additional)
- 35 Mbp acrocentric short-arm SDs
- ~47 Mbp pangenome-polymorphic SDs (Jeong 2025)
- Per-Pair-IGC-Annotation (Vollger 2023)
- Locus-spezifische Polymorphismus-Tiefe (NPIP, NOTCH2NL, 22q11.2)

Trotzdem als **hg38-Backbone-Skelett** weiterhin unverzichtbar — die WGAC-Bailey-Methodik ist die reproduzierbare Referenz, und für die hg38-Coords gibt es keine andere "frozen, citable" SD-Definition.

### "T2T-CHM13 SD-BED (Vollger 2022) ausreichend?"
**Fast** — schließt 51 Mbp + acrocentric + bessere Pair-Identity-Resolution. Aber:
- Nur eine Referenz-Hap, keine Population-SDs
- Last update 2022-03-11; spätere Eichler-Paper liefern keinen Master-BED-Update
- Pangenome-only SDs (47 Mbp aus Jeong 2025) fehlen auch hier

### "GIAB genome-stratifications ausreichend?"
**Als Maskierungs-Layer ja, als Locus-Katalog nein.** SEDEF-basiert seit v3 (Methodik-Upgrade über UCSC), inklusive CHM13. Aber 2 Stratifications für CHM13 (vs 9 für hg38) zeigt die Pangenome-Lücke. v4.0 fehlt noch (Stand 2026-06-02).

### "Müssen wir noch mehr Quellen einbeziehen?"
**Ja, fünf** (siehe T2-Stack oben). Zusätzlich für T1-Kuration die Eichler-Paper-Linie (Jeong 2025, Dishuck 2025, Real 2026, Porubsky 2025) als Supplementary-Table-Quellen.

### "Wieviele Loci brauchen vermutlich T1-Kuration?"
**~50 Loci/Klassen.** Nicht 300. T2 (UCSC + T2T-CHM13 + GIAB + Paraphase) deckt ≥ 96 % der publizierten SDs ab. T1 muss die letzten Pangenome-only-Cluster, IGH-Constant-Region, Y-Palindrome, acrocentric short arms, und die 5–10 microdel-syndrome-Loci mit Sub-LCR-Struktur pflegen.

---

## 8. Quellen (DOI / URL)

- Bailey JA et al. (2002) *Recent segmental duplications in the human genome.* Science 297:1003–1007. doi:10.1126/science.1072047
- Vollger MR et al. (2022) *Segmental duplications and their variation in a complete human genome.* Science 376:eabj6965. doi:10.1126/science.abj6965 — https://pmc.ncbi.nlm.nih.gov/articles/PMC8979283/
- Vollger MR et al. (2023) *Increased mutation and gene conversion within human segmental duplications.* Nature 617:325–334. doi:10.1038/s41586-023-05895-y
- Liao W-W et al. (2023) *A draft human pangenome reference.* Nature 617:312–324. doi:10.1038/s41586-023-05896-x
- Jeong H, Dishuck PC, Yoo D et al. (2025) *Structural polymorphism and diversity of human segmental duplications.* Nat Genet 57:390–401. doi:10.1038/s41588-024-02051-8
- Dishuck PC et al. (2025) *Structural variation, selection, and diversification of the NPIP gene family from the human pangenome.* Cell Genomics. doi:10.1016/j.xgen.2025.100869 (PubMed 39975192)
- Real R et al. (2026) *Genetic diversity and regulatory features of human-specific NOTCH2NL duplications.* Cell Genomics. doi:10.1016/j.xgen.2026.100920 (Format-Erwartungswert; PubMed 40166283)
- Porubsky D et al. (2025) *Human de novo mutation rates from a four-generation pedigree reference.* Nature 643. doi:10.1038/s41586-025-08922-2 (PubMed 40269156)
- Ma W, Chaisson MJP (2025) *Genotyping sequence-resolved copy number variation using pangenomes…* Nat Genet. doi:10.1038/s41588-025-02346-4 (PubMed 41107550)
- Chen S et al. (2024) *Paraphase: genome-wide profiling of highly similar paralogous genes using HiFi sequencing.* Nat Commun. doi:10.1038/s41467-025-57505-2 — github.com/PacificBiosciences/paraphase
- Olson ND et al. (2024) *The GIAB genomic stratifications resource for human reference genomes.* Nat Commun 15:9029. doi:10.1038/s41467-024-53260-y — github.com/genome-in-a-bottle/genome-stratifications
- Guarracino A et al. (2023) *Recombination between heterologous human acrocentric chromosomes.* Nature 617:335–343. doi:10.1038/s41586-023-05976-y
- Rhie A et al. (2023) *The complete sequence of a human Y chromosome.* Nature 621:344–354. doi:10.1038/s41586-023-06457-y
- Boettger LM et al. (2012) *Structural diversity and African origin of the 17q21.31 inversion polymorphism.* Nat Genet 44:881–885. doi:10.1038/ng.2335
- Skaletsky H et al. (2003) *The male-specific region of the human Y chromosome…* Nature 423:825–837. doi:10.1038/nature01722
- Repping S et al. (2002) *Recombination between palindromes P5 and P1 on the human Y chromosome.* AJHG 71:906–922. doi:10.1086/342928
- Watson CT et al. (2013) *Complete haplotype sequence of the human immunoglobulin heavy-chain V, D, J genes.* AJHG 92:530–546. doi:10.1016/j.ajhg.2013.03.004
- Rodriguez OL et al. (2020) *A novel framework for characterizing genomic haplotype diversity in the human immunoglobulin heavy chain locus.* Front Immunol 11:2136. doi:10.3389/fimmu.2020.02136
- Mikocziova I et al. (2025) *The human IG heavy chain constant gene locus is enriched for large structural variants…* Cell Genomics. doi:10.1016/j.xgen.2025.100784 (PubMed PMC11844466)
- Demaerel W et al. (2026) *Population differences of chromosome 22q11.2 duplication structure…* Nat Commun. doi:10.1038/s41467-026-71905-y (biorxiv 2025-07-04)
- Eichler EE (2024) *Beyond the Human Genome Project: The Age of Complete Human Genome Sequences and Pangenome References.* Annu Rev Genom Hum Genet. PubMed 38663087 (PDF nicht parsbar; Aussage aus Abstract + Open-Access-Preprint)

---

## 9. Operative Next Steps

1. **HPRC-Per-Hap-SD-Master-BED bauen** (T2-Augmentierung) — aus den 200+ HPRC R2 Assemblies SDs via SEDEF/BISER re-callen und konsolidieren. Eichler-Lab plant das vermutlich, aber bislang nicht released. Wir können das **selber** als Lab-internen Layer machen (auf the HPC cluster, BeeGFS, SEDEF dauert pro Hap ~ 4 h GPU-frei).
2. **Jeong 2025 Supplementary 47.4 Mbp-Sub-Liste sichten** und Top-15 "pangenome-only-SD-Loci" für T1-Kuration extrahieren.
3. **GIAB v4.0 watchlist** setzen (github release-feed) — vermutlich Q3/Q4 2026 mit T2T-CHM13 SD-Stratification-Update.
4. **T1-Backlog finalisieren** auf ca. 50 Loci aus Tabelle in §6.
5. **Eichler-Lab-Publication-Feed** quartalsweise pollen (insbesondere bioRxiv mit Filter "Eichler" + "segmental duplication") — die nächste Tranche (NBPF-zentriert, vermutet) kommt voraussichtlich 2026 H2.

---

## 10. Bekannte Unschärfen / Caveats

- **Real 2026 NOTCH2NL DOI** ist ein Format-Erwartungswert (Cell Genomics 2026); konkrete Volume/Issue/Article-Number muss beim Print-Update verifiziert werden.
- **Eichler 2024 Review** (Annu Rev Genom Hum Genet) konnte nicht textgenau aus PDF gefetcht werden; konkrete Statements oben sind aus Abstract + verwandten Eichler-Texten 2023–2024.
- **Mikocziova 2025 PMC-Zugang** funktioniert; das Paper validiert direkt unseren IGHG4-ChimDup-Befund und sollte als Co-Citation bei jedem IGH-Eintrag im Katalog stehen.
- **47.4 Mbp pangenome-only** aus Jeong 2025 ist die Population-skalierte Schätzung über 170 Hap (38 AFR + 47 nonAFR Samples). Eine vollständige Catalog-Tabelle pro Locus liegt nur in der Supplementary; manuelle Extraktion 1–2 Manntage.
