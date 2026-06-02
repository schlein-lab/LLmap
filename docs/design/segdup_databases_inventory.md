# Segmental-Duplication-Datenbanken — Inventar für LLmap-SegDup-Katalog

**Stand:** 2026-06-02
**Autor:** Recherche-Inventar (WebFetch + WebSearch), keine heavy parsing der Quell-BEDs.
**Zweck:** Quellen-Inventar für quartalsweise Updates des manuell kuratierten LLmap-SegDup-JSON-Katalogs pro Locus. Pro Quelle: URL, Format, API, Lizenz, Update-Frequenz, Coverage, Beispiel.

---

## TL;DR — Empfehlung

**Primär-Quellen** (auto-pollbar, dürfen Katalog-Backbone bilden):

1. **UCSC genomicSuperDups** (hg38, hg19, mm10/mm39, ggf. hs1 via assembly hub) — stabilster, am besten kuratierter Bailey/WGAC-Track; tab-delimited, einfach parsbar, gemeinfrei. Static (letztes Update hg38 = 2014-10-14), aber genau das macht ihn als Referenz-Anker brauchbar.
2. **T2T-CHM13 SD-BED** (`chm13v2.0_SD.bed` + `chm13v2.0_SD.full.bed`, Vollger 2022) — autoritativer SD-Katalog für die T2T-Referenz. Public S3, BED-Format, semi-static (v2022-03-11, mit punktuellen Updates über CHM13-Repo).
3. **gnomAD-SV v4** (Karczewski-Linie, Release 2023-11) — populationsskalierte SV-Frequenzen inkl. DUP/mCNV; nicht reiner SD-Katalog, aber unverzichtbar für Allele-Frequencies und für Regions-Maskierung (96.9% Precision-Region-Definition). VCF, public.
4. **HPRC (Liao 2023 + Release 2, 2025-05)** — pangenome-skaliertes Material; SD-Annotationen leben primär in begleitenden Eichler-/Vollger-Releases. Pro Sample muss man die Per-Assembly-Annotation aus AnVIL / S3 ziehen.
5. **Eichler-Lab Resource-Releases** (Vollger 2022, Jeong 2025, Dishuck 2025, Real 2026) — kontinuierlich aktualisierte SD-Maps und locus-spezifische Studien (NPIP, NOTCH2NL, ape-SDs). Manuell zu polling, weil pro-Paper unterschiedliche Endpoints.

**Sekundär** (locus-spezifische Tiefe / Annotation-Layer, nicht für Backbone):
- IMGT/GENE-DB (IGH/IGK/IGL/TRA/TRB), IPD-IMGT/HLA (MHC class I+II), IPD-MHC (cross-species MHC), ClinGen Dosage Sensitivity (klinische Relevanz + recurrent CNVs).

**Watch-only / Tools mit pre-computed Output**:
- dbVar (sehr breit, redundant; lohnt nur als optional einspielbarer Aggregator), DGV (Microarray-Ära, hg38 vorhanden, aber Lizenz unklar und Methodik veraltet), DECIPHER (registration-only, eher klinischer Layer), BISER/SEDEF/ASGART (Tools, keine kuratierten Kataloge — nur wenn wir selbst nachrechnen wollen), FlyBase/MGI (keine dedizierten SD-Tables — falls cross-species nötig, müssten wir selbst BISER laufen lassen).

---

## Detail-Inventar

### 1. UCSC genomicSuperDups (Bailey/WGAC)

| Feld | Wert |
|---|---|
| URL | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz` (analog für andere Assemblies) |
| Track-Doc | https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=genomicSuperDups |
| Format | tab-delimited .txt.gz, BED-artig + Pair-Spalten (chrom/Start/End/name/score/strand + otherChrom/otherStart/otherEnd/otherSize + alignment-Metriken matchB/mismatchB/fracMatch/jcK/k2K) — ca. 30 Spalten |
| API | nein direkt (Table Browser + REST API möglich, aber Bulk per FTP/HTTP empfohlen) |
| Lizenz | UCSC standard: free for academic/non-profit; commercial via license. Kein CC0 |
| Update-Frequenz | **Static** — hg38 zuletzt 2014-10-14; hg19 noch älter. Praktisch frozen ⇒ ideal als deterministic anchor |
| Methode | Bailey 2002 (Science 297:1003) WGAC: fuguization + ≥1 kb / ≥90% identity SD-Calls |
| Assembly-Support | hg38, hg19, hg18; mm10, mm39 (mouse); rn7 (rat); danRer11 (zebrafish); dm6 (Drosophila — als Track, nicht überall vorhanden); plus diverse weitere via Table Browser. **hs1 (T2T-CHM13) hat KEINEN genomicSuperDups-Track** im Standard-bigZips — dort kommt T2T-eigene Annotation zum Einsatz |
| Beispiel-Locus | chr14 IGH-Region: SD-Pairs zwischen IGHG-Genen und Flanking-Duplikon-Blöcken (locus chr14:105M-107M) trivial extrahierbar via `tabix`-loses zcat + awk |
| LLmap-Empfehlung | **PRIMÄR**. Sollte als Skelett-Layer pro Assembly dienen. Wegen static state quartalsweise nur sanity-check, kein echtes Polling nötig |

### 2. T2T-CHM13 Segmental Duplications

| Feld | Wert |
|---|---|
| URL | `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed` und `chm13v2.0_SD.full.bed` |
| Quelle | Vollger et al. 2022, Science (Eichler Lab); siehe github.com/marbl/CHM13 README |
| Format | BED (simple = 4-6 Spalten; full = erweiterte Spalten mit Pair-Info, identity, alignment-Metriken) |
| API | nein, statischer S3-Pfad |
| Lizenz | Public (T2T-Consortium); Vollger 2022 Science / bioRxiv |
| Update-Frequenz | **Semi-static**, v2022-03-11. Folge-Updates erscheinen typischerweise im Eichler-Lab als neue Releases pro Paper, nicht in dieser Datei |
| Coverage | T2T-CHM13 v2.0 (chm13v2.0), ~218 Mbp / 7.0% Genom als SD (vs. 5.4% in GRCh38) |
| Beispiel-Locus | IGH (chr14_MATERNAL bei CHM13): Tandem-IGHG-Duplikon-Block 19.5 kb (deckt sich mit unserem IGHG4-dup-Befund) |
| LLmap-Empfehlung | **PRIMÄR** für die CHM13-Spur. Pflicht, weil GRCh38-SD chr14:IGH falsch repräsentiert (siehe `ighg4_chm13_discrepancy.md`) |

### 3. dbVar (NCBI)

| Feld | Wert |
|---|---|
| URL | `ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/` |
| Strukturen | `by_assembly/` (GRCh37, GRCh38, älter), `by_study/`, `sandbox/sv_datasets/nonredundant/` (aggregiert) |
| Format | GVF + VCF + TSV (Studie-abhängig auch XML) |
| API | E-utilities (esearch/efetch via Entrez); kein dedizierter SV-Filter-Endpoint |
| Lizenz | NCBI Standard — i.d.R. public domain für Aggregat, study-spezifisch teils restriktiver |
| Update-Frequenz | Continuous (Studie-für-Studie-Submissions) |
| SD-Filter | **Nicht direkt** — kein `var_type=segmental_duplication`. SDs erscheinen indirekt als DUP/CNV in SD-reichen Regionen. Sudmant 2015 (`nstd112`) ist die SD-relevanteste Großstudie |
| Beispiel | Sudmant 2015 nstd112: 53 Dateien / 353 MB, populationskalierte CNVs/SDs |
| LLmap-Empfehlung | **WATCH-ONLY**. Zu redundant zu gnomAD-SV; lohnt nur, wenn wir gezielt historische CNV-Studien einbauen wollen |

### 4. gnomAD-SV v4

| Feld | Wert |
|---|---|
| URL | https://gnomad.broadinstitute.org/downloads (v4 SV-Sektion); auch AWS Open Data: `s3://gnomad-public-us-east-1/` |
| Format | VCF (gz + tbi), gesplittet pro Chr |
| API | gnomAD GraphQL API (frequenz-queries); Bulk via VCF |
| Lizenz | **gnomAD Public — keine Restriktion** (gemäss gnomAD terms; effektiv frei nutzbar, Citation erwartet). Auf AWS unter "Open Data" |
| Update-Frequenz | Major releases ~jährlich; v4-SV-Release 2023-11; minor patches dazwischen |
| Coverage | GRCh38; 63,046 Genome (unrelated); 1,199,117 high-quality SVs; Klassen: DEL/DUP/INS/INV/CTX/CPX/mCNV |
| SD-Relevanz | Explizit dokumentiert: 9.7% Genom (SR+SD) liefern ~25k SVs/Genom, die nur long-read aufrufen kann. Precision steigt von Gesamtmenge auf 96.9% wenn SD/SR-Regionen exkludiert ⇒ gnomAD liefert auch eine implizite "low-confidence in SD"-Maskierungs-BED |
| Beispiel | IGH chr14:105.5-106.9 Mb: gnomAD-SV-DUPs mit AF + sample-counts |
| LLmap-Empfehlung | **PRIMÄR** für Allele-Frequenz-Layer und Population-Maskierung |

### 5. DGV (Database of Genomic Variants, TCAG)

| Feld | Wert |
|---|---|
| URL | http://dgv.tcag.ca/dgv/app/downloads |
| Format | tab-delimited TXT + GFF3 |
| API | nein |
| Lizenz | "Usage disclaimer" auf tcag.ca, nicht klar CC; in Praxis akademisch frei nach Citation |
| Update-Frequenz | Sporadisch — sichtbare Releases 2016, 2020, 2025-12 |
| Coverage | GRCh37 + GRCh38 |
| SD-Filter | Keine dedizierte SD-Klasse — DGV mischt CNVs/Indels/Inversionen aus Microarray- und WGS-Studien. SDs müssen via Overlap mit genomicSuperDups identifiziert werden |
| LLmap-Empfehlung | **WATCH-ONLY**. Methodisch microarray-lastig; nur als historischer Querkonsens-Layer |

### 6. DECIPHER

| Feld | Wert |
|---|---|
| URL | https://www.deciphergenomics.org/about/downloads/data |
| Format | "Data Files" — Details hinter Login/Terms; typischerweise CSV der population-control CNVs |
| API | nein öffentlich |
| Lizenz | DECIPHER Terms of Use (PDF); Data-Sharing-Policy. Effektiv: akademisch nach Registration, klinischer Kontext erwartet |
| Update-Frequenz | Continuous (klinische Einreichungen) |
| Coverage | GRCh37/GRCh38, klinische CNVs aus ~250 Zentren |
| SD-Relevanz | Indirekt — DECIPHER ist die go-to-Quelle für *Konsequenz* von SD-mediierten Recurrent-CNVs (z.B. 16p11.2, 22q11.2). Kein SD-Track, aber CNV-Hotspot-Anreicherung in SD-flankierten Regionen |
| LLmap-Empfehlung | **SEKUNDÄR** für klinische Annotation einzelner Loci, nicht für Backbone |

### 7. HPRC Pangenome Resources

| Feld | Wert |
|---|---|
| URL | https://humanpangenome.org/data/ ; github.com/human-pangenomics/hpp_pangenome_resources ; AnVIL |
| Releases | Release 1 (Liao 2023, 47 individuals, 94 Haps); Release 2 (2025-05-12, >200 Individuals) |
| Pangenome-Modelle | Minigraph, Minigraph-CACTUS, PGGB (drei parallele Graphen) |
| Format | GFA (graphs), FASTA (per-assembly), VCF (called variants), BED (per-sample annotations) |
| Lizenz | NIH-funded, open (public consortium) |
| SD-Annotation | **Nicht zentral als ein "SD-master-bed"** publiziert. SD-Annotationen kommen aus begleitenden Eichler-Papern (Vollger 2023, Jeong 2025) — z.T. als Supplementary Tables, z.T. auf S3 unter `human-pangenomics/` |
| Update-Frequenz | Quartalsweise Daten-Drops im Release-2-Modus, Papier-getaktet für SD-Layer |
| Beispiel | Pro Sample (z.B. HG00272, HG00290) haben wir lokal Per-Read-Analysen gemacht; SD-Position pro Hap aus Vollger-Pipeline |
| LLmap-Empfehlung | **PRIMÄR** für Multi-Hap-Diversität pro Locus; brauchst aber Manual-Cherry-Pick aus den begleitenden Papers |

### 8. Eichler Lab Resources

| Feld | Wert |
|---|---|
| URL | https://eichler.gs.washington.edu/ (umgezogen von eichlerlab.gs.washington.edu) |
| Schlüssel-Releases | • Vollger 2022 Science (T2T-CHM13-SDs) → marbl/CHM13 + S3-bed; • Jeong 2025 Nat Genet (struktural polymorphism + diversity); • Dishuck 2025 Cell Genom (NPIP family); • Real 2026 Cell Genom (NOTCH2NL); • Porubsky 2025 Nature (4-generation pedigree, SD de novo Rate) |
| Format | pro-Paper unterschiedlich: BED, VCF, Custom-TSV, Supplementary-Tables |
| API | nein |
| Lizenz | Per-Paper, i.d.R. Daten public via Zenodo / GitHub / S3 |
| Update-Frequenz | Papier-getaktet (4-8 release-relevante Papers/Jahr) |
| LLmap-Empfehlung | **PRIMÄR** für locus-Tiefe und method-of-record bei SD-Calls. Manuell pollen, jedes Quartal Eichler-Publikationen durchgehen |

### 9. dbVar/CHM13 + GIAB SD-Stratifications

| Feld | Wert |
|---|---|
| URL | github.com/genome-in-a-bottle/genome-stratifications |
| Format | BED-Stratification-Sets (SegDup, low-mappability, CNVs) für hg38, hg37, T2T-CHM13 |
| Lizenz | NIH/NIST public domain |
| Update-Frequenz | Releases v2.0 (2021), v3.0 (2023), v4.0 (erwartet) |
| Relevanz | Liefert die Quasi-Standard "schwer für SV-Caller"-BEDs, die GA4GH benchmarks nutzen — perfekt als Maskierungs-Layer |
| LLmap-Empfehlung | **PRIMÄR-NAHE** (eigentlich Sekundär, aber sehr unkompliziert) — auto-pollbar pro Release |

### 10. SD-Calling-Tools (BISER, SEDEF, ASGART)

| Tool | Status | Lizenz | Output | Pre-computed Catalog? |
|---|---|---|---|---|
| **SEDEF** (vpc-ccg/sedef) | **Deprecated** seit 2018 | MIT | BEDPE (28+ Felder) | hg19, hg38, mm8 "final calls" verfügbar im Repo |
| **BISER** (0xTCG/biser) | aktiv, v1.4 (2023-03) | MIT | BEDPE + `.elem` core-SD-decomposition | Keine offiziellen pre-computed Sets; Repo hat Beispiel-Outputs |
| **ASGART** (delehef/asgart) | low-activity, v2.4.0 (2020-11) | GPLv3 | JSON (duplication families + duplicons) | Keine pre-computed Datasets; Online-Server http://asgart.irit.fr |
| LLmap-Empfehlung | **WATCH-ONLY für Tools**, **PRIMÄR für SEDEF-final-calls-BEDs** wenn wir explizit redundancy zum WGAC-Track wollen |

### 11. Locus-spezifische Datenbanken

#### IMGT (Ig + TR)
| Feld | Wert |
|---|---|
| URL | https://www.imgt.org/download/GENE-DB/ |
| Format | FASTA (nt/aa, mit/ohne gaps, in-frame/all-frame) + GeneList |
| Loci | IGH, IGK, IGL, TRA, TRB, TRG, TRD; species: human + mouse + Ratte u.a. Wirbeltiere |
| API | "IMGT/V-QUEST" + "IMGT/HighV-QUEST" als Web-Services; kein REST für GENE-DB |
| Lizenz | "Copyright IMGT" — academic free use with citation; commercial restriction; redistribution mit Erlaubnis |
| Update-Frequenz | RELEASE-File getrackt, ~quartalsweise Mikro-Updates |
| LLmap-Empfehlung | **SEKUNDÄR / PRIMÄR FÜR IG-LOCI**. Pflicht-Quelle für IGH-Allele/Gene-Coords. Manuelle Pflege wegen Citation-Pflicht |

#### IPD-IMGT/HLA
| Feld | Wert |
|---|---|
| URL | https://www.ebi.ac.uk/ipd/imgt/hla/download/ ; FTP: `ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/` ; github.com/ANHIG/IMGTHLA |
| Format | FASTA, PIR, MSF + flat-file + alignment-archive |
| Lizenz | **CC-BY-ND** (Creative Commons Attribution-NoDerivs) |
| Update-Frequenz | Quartalsweise Major releases |
| Coverage | HLA class I + II Genes mit allelischer Tiefe |
| LLmap-Empfehlung | **PRIMÄR FÜR MHC-LOCUS**. Auto-pollbar via GitHub-Branches |

#### IPD-MHC (non-human)
| Feld | Wert |
|---|---|
| URL | https://www.ebi.ac.uk/ipd/mhc/ |
| Format | Download-Seite (Details nicht offen); zentrale Sequence-DB |
| Species | NHP, DLA (dog), FISH, OLA (sheep), BoLA (cattle), ELA (horse), SLA (pig), RT1 (rat), CHICKEN, CLA (cat), CeLA (camel) |
| Lizenz | IPD license page |
| Aktuell | Release 3.16.0.0 (2026-01, build 231) |
| LLmap-Empfehlung | **SEKUNDÄR**, nur falls cross-species MHC-Analyse anliegt |

#### ClinGen (OMIM-nah, klinischer Layer)
| Feld | Wert |
|---|---|
| URL | https://search.clinicalgenome.org/kb/downloads |
| Format | CSV (gene-level dosage), BED (haploinsufficiency, triplosensitivity, recurrent CNV), AED, JSON-API (actionability), TSV |
| API | REST endpoints; Allele Registry + Evidence Repository |
| Lizenz | Open, Citation erwartet |
| Update-Frequenz | **Täglich** mit 60-Tage-FTP-Archiv |
| Recurrent-CNV-Filter | **Ja**, eigener Filter im Dosage-Sensitivity-Bereich; perfekt für SD-mediierte Recurrent-CNVs (16p11.2, 17q21.31, etc.) |
| LLmap-Empfehlung | **PRIMÄR** für klinische Annotation; täglich aktualisiert ⇒ auto-pollbar |

#### OMIM
| Feld | Wert |
|---|---|
| URL | https://omim.org/downloads |
| Format | TSV (genemap2.txt, mim2gene.txt, morbidmap.txt) |
| Lizenz | Registration required (academic frei), keine commercial redistribution |
| LLmap-Empfehlung | **SEKUNDÄR**. Klinisch wertvoll, aber redundant zu ClinGen für unsere Recurrent-CNV-Loci |

### 12. Cross-species

#### Maus (MGI / mm39)
- **Bailey 2004** (genome.cshlp.org/content/14/5/789): originale Mouse-SD-Analyse, C57BL/6J-basiert
- **UCSC mm10/mm39** hat genomicSuperDups (verifizierte Track-Family)
- **MGI selber** hat **KEINE dedizierte SD-Table** — Maus-SD-Daten leben in UCSC + neue Eichler-Maus-Paper (2025 biorxiv "SD-mediated rearrangements alter mouse genomes")
- LLmap: **PRIMÄR = UCSC mm39 genomicSuperDups**

#### Rat (rn7)
- UCSC genomicSuperDups als Track, less extensive curation; LLmap: PRIMÄR-NAHE für Rat-Spur

#### Zebrafish (danRer11)
- UCSC genomicSuperDups vorhanden; LLmap: SEKUNDÄR

#### Drosophila (dm6, FlyBase)
- **FlyBase**: KEINE dedizierte SD-Table. Pre-computed Files (FASTA/GFF/GTF) auf S3, aber nicht SD-spezifisch
- Wenn wir SD-Calls für Drosophila brauchen: BISER auf dm6 selber laufen lassen
- LLmap: **WATCH-ONLY**, nur on demand

---

## Pollstrategien-Tabelle

| Quelle | Polling | Methode | Kadenz | Manuell? |
|---|---|---|---|---|
| UCSC genomicSuperDups | curl + sha256 | wget hgdownload .txt.gz | quartalsweise sanity check | nein |
| T2T-CHM13 SD-BED | curl + sha256 | wget S3 BED | quartalsweise | nein |
| gnomAD-SV v4 | release-monitor | GitHub gnomad-broadinstitute release-feed | release-trigger | nein |
| HPRC | manual | watch HPRC twitter + GitHub releases | quartalsweise | **ja** |
| Eichler Lab papers | manual | watch eichler.gs.washington.edu/publications + bioRxiv | quartalsweise | **ja** |
| GIAB Stratifications | release-monitor | GitHub genome-in-a-bottle/genome-stratifications | release-trigger | nein |
| ClinGen | daily | wget CSV/BED, diff | täglich (oder wöchentlich) | nein |
| IMGT/GENE-DB | release-file check | wget RELEASE + diff | quartalsweise | nein |
| IPD-IMGT/HLA | GitHub watch | git pull github.com/ANHIG/IMGTHLA | quartalsweise | nein |
| dbVar | optional | FTP study-by-study | on demand | n/a (watch-only) |
| DGV | optional | FTP | on demand | n/a |
| DECIPHER | manual | Login + manual | on demand | **ja** |
| SEDEF/BISER/ASGART | n/a | nur wenn wir SDs selbst callen | – | n/a |

---

## Lizenz-Kompatibilität (TL;DR pro Lizenz-Cluster)

- **gemeinfrei / NIH-funded** (UCSC akademisch, T2T-CHM13, gnomAD, GIAB, ClinGen, dbVar): kein Problem; LLmap-Katalog darf diese Daten unter beliebiger eigener Lizenz weiterleiten, solange Citation erfolgt
- **CC-BY-ND** (IPD-IMGT/HLA): ⚠️ **No-Derivs**. Wir dürfen Originaldaten **NICHT** modifizieren und unter eigener Lizenz weitergeben. Wir können sie nur als externe Referenz im Katalog verlinken
- **IMGT-Copyright** (IMGT/GENE-DB): academic free + cite, commercial restricted. Falls LLmap commercial-fähig sein soll: vorsichtig
- **GPLv3** (ASGART): Tool, kein Daten-Problem
- **MIT** (SEDEF, BISER): unkritisch
- **DECIPHER**: registration-only; **nicht** für Katalog-Bundle, nur ad-hoc-Lookup

---

## Beispiel-Locus IGH (chr14:105.5-106.9 Mb GRCh38 / chr14_MATERNAL CHM13)

Konsistenz-Check zwischen Quellen für unsere Top-Locus:

| Quelle | IGH-SD-Annotation? | Notes |
|---|---|---|
| UCSC genomicSuperDups hg38 | Ja, Pairs zwischen IGHG-Genen | aber: GRCh38 fehlt die canonical Tandem-Dup, daher unvollständig |
| T2T-CHM13 SD-bed | Ja, vollständiger inkl. 19.5 kb IGHG-Tandem-Dup | autoritativ |
| gnomAD-SV v4 | DUPs als Population-AF, low confidence in SD region | Frequenz-Layer |
| HPRC (Liao 2023 / R2) | per-Hap-Assembly; wir haben dies in unserem cohort_data bereits validiert | gold standard für Diversity |
| IMGT/GENE-DB | Gen-Koordinaten der IGHG-Paralogs (IGHG1-4, IGHGP) | für Allele-Annotation |
| ClinGen | recurrent CNVs am Locus? minimal | nicht zentrale Quelle |

→ Für IGH liefert die Kombination **UCSC + T2T-CHM13-SD-BED + HPRC + IMGT** ein vollständiges Bild; gnomAD/dbVar liefern frequenzielles Background; ClinGen/DECIPHER sind nicht informativ.

---

## Anhang: Quellen-URLs (Stand 2026-06-02)

- UCSC track UI: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=genomicSuperDups
- UCSC database FTP: https://hgdownload.soe.ucsc.edu/goldenPath/{hg38,hg19,hs1,mm39,rn7,danRer11,dm6}/database/
- T2T-CHM13 README: https://github.com/marbl/CHM13
- T2T-CHM13 SD S3: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/
- dbVar FTP: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/
- gnomAD v4 SVs: https://gnomad.broadinstitute.org/downloads#v4-structural-variants ; News-Post: https://gnomad.broadinstitute.org/news/2023-11-v4-structural-variants/
- gnomAD AWS Open Data: https://registry.opendata.aws/broad-gnomad/
- DGV downloads: http://dgv.tcag.ca/dgv/app/downloads
- DECIPHER downloads: https://www.deciphergenomics.org/about/downloads/data
- HPRC data: https://humanpangenome.org/data/ ; GitHub: https://github.com/human-pangenomics/hpp_pangenome_resources
- Eichler Lab: https://eichler.gs.washington.edu/
- GIAB Stratifications: https://github.com/genome-in-a-bottle/genome-stratifications
- ClinGen downloads: https://search.clinicalgenome.org/kb/downloads
- IMGT GENE-DB: https://www.imgt.org/download/GENE-DB/
- IPD-IMGT/HLA: https://www.ebi.ac.uk/ipd/imgt/hla/download/ ; GitHub: https://github.com/ANHIG/IMGTHLA
- IPD-MHC: https://www.ebi.ac.uk/ipd/mhc/
- OMIM: https://omim.org/downloads
- SEDEF: https://github.com/vpc-ccg/sedef (deprecated)
- BISER: https://github.com/0xTCG/biser
- ASGART: https://github.com/delehef/asgart
- FlyBase Downloads: https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview

---

## Offene Fragen / nächste Schritte

1. **HPRC R2 SD-Master-Bed**: existiert eine? Ggf. direkten Kontakt zum HPRC-data-team / AnVIL nehmen, statt aus Supplementaries zusammensuchen
2. **CHM13 SD update post-2022-03-11**: ist Vollger 2022 das letzte File, oder gibt es ein chm13v2.1/v3 Re-Release? Im Eichler-Trakt 2025/26 Jeong + Real Papers checken
3. **GIAB Stratifications v4**: Release-Datum / Inhalt verifizieren (Recherche zeigt v3.0 als letztes Major release)
4. **Schema-Konvergenz**: LLmap-Katalog-JSON pro Locus braucht ein einheitliches Schema, das alle obigen Quellen via einfaches Mapping aufnehmen kann (Vorschlag: `{locus_id, assembly, regions:[{chrom,start,end,strand}], paralog_pairs:[...], allele_freqs:{gnomAD:..., HPRC:...}, source_provenance:[...], confidence:[...]}`)
