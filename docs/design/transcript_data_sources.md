# Transcript Data Sources for IGHG4-ChimDup Silence Test

**Status:** Recherche-Snapshot 2026-06-02
**Kontext:** HPRC iso-seq LCLs (49 Samples) sind klonal-isotype-selektiert (siehe `hprc_lcl_clonality_blocks_ighg4_test`-Memory). Für **polyklonale IGHG4-Expressions-Tests** und insbesondere für die Frage "ist ChimDup-Locus transkriptionell stumm vs. canonical?" brauchen wir andere Datenquellen — PBMC, Tonsillen, sortierte B-Zellen, idealerweise mit BCR-Pairing.

---

## TL;DR

**(a) Lokal auf the HPC cluster ready-to-use:**
- 1× PBMC PacBio Iso-Seq (BioIVT anonymer Spender, 13 GB BAM, 44.3M dedup reads, 18 IGHG4-Reads in 12 cells — **methodisch validiert, statistisch grenzwertig**) — `/beegfs/u/<user>/ighg4_dv/pbmc_pacbio/`
- 49× HPRC LCL iso-seq (363 GB, **klonal-biased, für IGHG4 unbrauchbar** — siehe Memory `hprc_lcl_clonality_blocks_ighg4_test`)
- Mausembryos GSE119945 + E-MTAB-6967 (Cao/Pijuan-Sala, scrnaseq, **nicht relevant** — Maus, kein IGH)
- NRW-Drop M47835/M47836/M47837 ONT (5.3 TB Roh): **DNA, kein RNA-seq, irrelevant**

**Lokales Fazit:** Effektiv 1 polyklonale PBMC-Probe vorhanden. Für seriöse ChimDup-Silence-Statistik brauchen wir mindestens 5-10 zusätzliche polyklonale Samples.

**(b) Top-3 Empfehlung (downloaden):**

| Rang | Accession | Typ | Reason |
|------|-----------|-----|--------|
| 1 | **PRJNA922682** (FLAIRR-seq) | PacBio Iso-Seq targeted IGH, 10 healthy PBMC + B-cell donors | **Goldstandard** für IGHG4 — voller IGH-Konstantbereich (CH1-3 + Hinge), allele-resolved Isotype, Q60 accuracy. Kleinste Download-Last (~15 GB). |
| 2 | **EGAS00001006779** (Pathogen-PBMC) | PacBio Iso-Seq bulk, 1 donor, 5 conditions | Voll-Length PBMC-Iso-Seq, deep coverage, *aber* nur 1 Donor → keine zwischen-Donor-Statistik. Ergänzt FLAIRR. EGA controlled access (DAC) — Antrag nötig. |
| 3 | **PacBio Kinnex PBMC public** (`DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p`) | PacBio Revio Kinnex 10x 5'/3', 20k cells | Single-cell, frei downloadbar von pacbcloud.com, ~155 GB per rep. Validiert IGHG4-cell-counts auf vergleichbarer Tech wie unsere BioIVT PBMC. |

**(c) Wieviel muss downloaded werden:**
- FLAIRR-seq PRJNA922682: ~15 GB (30 Gbases via sra-toolkit, 46 SRA experiments) — **~6-12h** auf the HPC cluster SLURM
- Kinnex PBMC public 1 rep: ~155 GB (BAM + FASTA + index) — **24h** wget/aria2c
- Pathogen-PBMC EGAS00001006779: zuerst DAC-Antrag (Wochen), dann ~30-50 GB

**(d) IgG4-RD Datasets:** Ja, mehrere relevant.
- **HRA003750** (Lu et al. 2023 *JCI Insight*): BD Rhapsody scRNA + AbSeq, 9 IgG4-RD Patienten + 7 Controls, 61k PBMCs — **controlled access** (West China Rheumatology DAC). Strategischer Wert für ChimDup-Selektions-Analyse: IgG4-RD-Patienten zeigen IGHG4-Expansion → wenn ChimDup-Carrier in der IgG4-RD-Kohorte überrepräsentiert sind, wäre das ein Selektions-Signal.
- Mehrere weitere IgG4-RD scRNA-seq Studies (Submandibulardrüse PMC10877852, peripheral T-helper PMC10530310) — alle 10x oder BD-Rhapsody-Plattform, short-read, nicht ideal für IGHG4-CH1-SNP-Diskriminierung aber für **B-cell-Counting + isotype label** brauchbar.

---

## 1. the HPC clusterInventar

### 1.1 PBMC / B-Zell-RNA

| Pfad | Inhalt | Größe | Tech | Status |
|------|--------|-------|------|--------|
| `/beegfs/u/<user>/ighg4_dv/pbmc_pacbio/` | BioIVT PBMC, scisoseq.5p--3p dedup BAM, 44.3M reads, 39.5k IGH-primary, 18 IGHG4-reads in 12 cells | 13 GB BAM + 234 MB index | PacBio Sequel single-cell iso-seq (pb_sc_isoseq, isoseq groupdedup v4.1.1) | Paralog-mapping + kmer-paralog TSV bereits berechnet (siehe `PBMC.kmer_paralog.tsv`) |
| `/beegfs/u/<user>/<group>/shared/references/transcriptome_longread/hprc_isoseq/` | 20 HPRC LCL-FLNC BAMs (+a13_work) | 363 GB | PacBio HiFi iso-seq | **Klonale Verzerrung** → für IGHG4 nicht brauchbar |
| `/beegfs/u/<user>/<group>/<lab>/ighg_rna_analysis/` | 49 HPRC samples mit per-sample IGHG-Mapping + DV-Calls | inferred ~100 GB | minimap2 → IGHG-ref + DeepVariant | Selbe LCL-Klonalitäts-Limitation |
| `/beegfs/u/<user>/ighg4_dv/isoseq_*` | Mapping/Analysis-Outputs (paralog_aln, chr14_aln, padded, personalized, NanoCount, per_paralog) | varies | nach-HPRC-Pipeline | Derivative von HPRC LCLs — nicht eigenständige Daten |

### 1.2 Single-Cell RNA-Seq (nicht IGH-fokussiert)

| Pfad | Inhalt | Größe | Relevanz |
|------|--------|-------|----------|
| `/beegfs/u/<user>/<group>/shared/datasets/scrnaseq/cao_2019/` (GSE119945) | Mouse organogenesis sci-RNA-seq3 | 0 (sra_cache leer, fastq dir leer — Download nie abgeschlossen) | **Irrelevant** (Maus) |
| `/beegfs/u/<user>/<group>/shared/datasets/scrnaseq/pijuan_sala_2019/` (E-MTAB-6967) | Mouse gastrulation 10x | 0 (Download nie abgeschlossen) | **Irrelevant** (Maus) |
| `/beegfs/u/<user>/<group>/shared/datasets/geo/GSE132042/` | Tabula Muris Senis multi-organ 10x | 20 GB | **Irrelevant** (Maus, sowieso keine Plasma/Mature-B-Zellen) |

### 1.3 NRW-Drop (Kubisch)

| Pfad | Inhalt | Größe | Relevanz |
|------|--------|-------|----------|
| `/beegfs/g/kubisch/<user>/incoming/datastorage_nrw_bc813032/raw/M47835/` | ONT pod5 (PromethION) für Fam47835 | 5.3 TB raw | **Irrelevant für IGHG4-RNA** (DNA-Sequenzierung) |
| `/beegfs/g/kubisch/<user>/incoming/datastorage_nrw_bc813032/Analysis/M478*LR/` | wf-human-variation Outputs | ~280 GB | DNA-SV/SNV-Calls, nicht Transcript |

### 1.4 sra-toolkit Status

- `sra_toolkit` Conda-Env: **bestätigt vorhanden** (in `_scripts/download_scrnaseq.slurm` aktiviert)
- Download-Workflow: `prefetch <accession> -O <OUTPUT_DIR> --max-size 1000000000 -a ascp` → `fasterq-dump -e 8 ... --split-files`
- BeeGFS quota: 2.4 PB frei, kein limit-issue
- Allokationen: `your_slurm_account` (8 CPUs / 8h job) etabliert

---

## 2. GEO/SRA-Kandidaten — publizierte Datasets

### 2.1 Top-Empfehlungen mit Detailmetadaten

#### **(A) FLAIRR-seq — PRJNA922682** (Ford et al. 2023, *J Immunol* 210:1607)

| Attribut | Wert |
|----------|------|
| Accession | **PRJNA922682** (12 BioSamples, 46 SRA experiments) |
| Sample types | PBMC (healthy), purified B cells (healthy), whole blood (1× COVID-19 time course) |
| Sample count | 10 healthy donor PBMCs + B-Zell-Replikate + 1 COVID-Verlauf |
| Tech | PacBio HiFi (CCS), 5'-RACE targeted IGH amplification, full-length IgG + IgM transcripts (~1.5 kb) |
| Bases | 30 Gbases total |
| Compressed size | **14.96 GB** |
| Coverage | IGHV+IGHD+IGHJ+**IGHC** (volles CH1-CH3 + Hinge mit Isotype-Auflösung); Q60 mean accuracy → SNP-level für CH1-Diskriminierung möglich |
| Access | **Public NCBI SRA** (prefetch + fasterq-dump direkt) |
| License | NCBI open |
| Citation | DOI: 10.4049/jimmunol.2200825 |

**Warum #1:** Targeted IGH-Iso-Seq mit voller Konstant-Region-Coverage, expliziter IGHG-Subisotyp-Auflösung im Paper, PBMC-polyklonal — exakt unsere Frage. Klein genug für Download in einer Nacht.

#### **(B) Pathogen-PBMC Iso-Seq — EGAS00001006779** (Boahen et al. 2024, *iScience*)

| Attribut | Wert |
|----------|------|
| Accession | **EGAS00001006779** (PacBio long-read) + EGAS50000000007 (QuantSeq short-read), PXD045237 (Proteomik) |
| Sample types | PBMC, healthy donors, 5 conditions (LPS, Poly(I:C), heat-inactivated *S. aureus*, *C. albicans*, RPMI control) |
| Sample count | **1 donor** für long-read, 5 donors für QuantSeq + Proteomik |
| Tech | PacBio Iso-Seq (full-length) — Sequel/Revio |
| Coverage | Full transcriptome, kein Targeted — IGHG4-Coverage abhängig von B-Zell-Anteil in PBMC und Stimulation |
| Access | **EGA controlled** — Data Access Committee (DAC) Antrag notwendig |
| Size | Inferred 30-50 GB |
| Citation | DOI: 10.1016/j.isci.2024.110525 |

**Warum #2:** Voll-Length-Iso-Seq mit deep PBMC-Coverage und Stimulations-Bedingungen (Pathogen-Stimulus → erwartete Class-Switch-Aktivierung → ggf. erhöhte IGHG4-Reads). EGA-DAC-Antrag dauert 2-6 Wochen.

#### **(C) PacBio Kinnex PBMC Tutorial — `DATA-Revio-Kinnex-PBMC-20kcells`** (PacBio Public)

| Attribut | Wert |
|----------|------|
| URL | `https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p-rep1/` |
| Sample types | PBMC (Standard 10x Chromium GEM-X 3' library, Kinnex-amplified) |
| Sample count | 2 replicates 20k cells + 1× 10k cells (3p), zusätzlich 5p-Variante und Parse-Variante |
| Tech | PacBio Revio + Kinnex single-cell RNA kit |
| Coverage per rep | Dedup BAM 61 GB, mapped BAM 34 GB, fasta 45 GB, total ~155 GB per rep |
| Access | **Public via HTTPS** (downloads.pacbcloud.com) |
| License | PacBio public data — Citation-recommended |
| Cells | 20 000 (Standard PBMC pre-mix) |

**Warum #3:** Single-cell + voll-Length-Transcripts, vergleichbar mit unserer BioIVT PBMC aber 2 Replikate + frische Tech (Revio + Kinnex statt Sequel + MAS-Seq). Free download, kein Antrag.

### 2.2 Ergänzende Kandidaten (Tier 2)

#### **(D) Tonsil Atlas — E-MTAB-13687** (Massoni-Badosa et al. 2024, *Immunity*)

| Attribut | Wert |
|----------|------|
| Accession | **E-MTAB-13687** (ArrayExpress) |
| Sample types | Whole tonsil tissue (Mandeln), 17 donors (10 discovery + 7 validation), Kinder/Adulte mix |
| Sample count | 17 donors, 357 206 cells |
| Tech | 10x Chromium scRNA-seq + scATAC-seq + Multiome + CITE-seq + **scBCR-seq** + spatial |
| Coverage | scBCR liefert per-cell Isotype-Label inklusive IGHG4 — aber 10x short-read, CH1-SNPs nicht direkt diskriminierbar |
| Access | ArrayExpress raw + Zenodo (processed Seurat) + Bioconductor `HCATonsilData` |
| License | Open with attribution |
| Citation | DOI: 10.1016/j.immuni.2024.01.006 |

**Strategischer Wert:** Großer polyklonaler **Tonsillen**-Datensatz mit pairing zwischen BCR-Isotype und Cell-State. Limitation: short-read, kann zwischen canonical-IGHG4 und ChimDup-IGHG4 (CH1-/Hinge-SNPs) wahrscheinlich nicht diskriminieren ohne Custom-Probe-Design. Aber: **Cohort-Background** für Allel-Frequenz-Vergleich.

#### **(E) IgG4-RD scRNA-seq — HRA003750** (Lu et al. 2023, *JCI Insight*)

| Attribut | Wert |
|----------|------|
| Accession | **HRA003750** (GSA-Human, NGDC China) |
| Sample types | PBMC, 9 treatment-naive IgG4-RD + 7 HC |
| Sample count | 16 donors, 61 379 PBMCs |
| Tech | BD Rhapsody + AbSeq (CITE-protein) |
| Coverage | scRNA-seq, 3'-tag-basiert — IGHG4-CH1-SNP-Diskriminierung schwierig |
| Access | **GSA-Human controlled access** — HDAC001850 (West China Rheumatology) |
| License | Restricted, DAR-Antrag notwendig |
| Citation | DOI: 10.1172/jci.insight.167602 |

**Strategischer Wert:** **IgG4-RD-Phänotyp-Kohorte**. Wenn ChimDup-Haplotyp ein Selektions-Treiber für IgG4-RD ist, müssten ChimDup-Carrier-Frequenz in dieser Kohorte ggü. HPRC angereichert sein. Aktuell nicht für direkten Expression-Test geeignet (3'-tag short-read).

#### **(F) IGHG4 DRG B-cells — medRxiv 2024.11.12.24317004** (Mancilla Moreno et al.)

- Sample: 8 Patienten chronischer Halsschmerz, DRG-Tissue (nicht PBMC/Tonsille)
- 1 Patient mit hoher IGHG4-Expression im C2 DRG (B-Zell-Infiltration)
- Tech: Visium spatial + scRNA-seq
- Accession: nicht öffentlich gefunden — Data-Availability-Section in PDF auslesen
- **Sehr klein**, eher Bestätigung-im-Tissue als statistische Power

### 2.3 Nicht empfohlen (gecheckt, ungeeignet)

- **GSE119507** Single-cell tonsil **CD4 T-Zellen** — keine B-Zellen-Coverage
- Mehrere IgG4-RD-Studies auf Submandibulardrüse — short-read tissue, kein CH1-Diskriminierungs-Potenzial
- GSE174188 SLE PBMCs — short-read scRNA-seq, nicht IGH-fokussiert

---

## 3. Download-Workflows

### 3.1 SRA (PRJNA922682)

```bash
# On the HPC cluster via SLURM (your_slurm_account, ~6-12h, 50 GB Cache)
sbatch --account=your_slurm_account --partition=std --cpus-per-task=8 --time=12:00:00 \
       --wrap="
source ~/miniforge3/etc/profile.d/conda.sh
conda activate sra_toolkit
OUT=/beegfs/u/<user>/<group>/shared/datasets/sra/PRJNA922682
mkdir -p \$OUT
# 1. SRR-Liste via eutils
esearch -db sra -query 'PRJNA922682[BioProject]' | efetch -format runinfo | csvcut -c 1 | tail -n+2 > \$OUT/srr_list.txt
# 2. prefetch alle Runs
prefetch --option-file \$OUT/srr_list.txt -O \$OUT/cache --max-size 50g
# 3. BAM aus PacBio HiFi extrahieren (NICHT fasterq-dump — wir wollen unaligned BAM)
for srr in \$(cat \$OUT/srr_list.txt); do
  sam-dump --aligned-region none --no-aligned \$OUT/cache/\$srr/\$srr.sra | samtools view -bS - > \$OUT/\$srr.bam
done
"
```

### 3.2 PacBio Public Kinnex (HTTPS)

```bash
# Direct wget, ca. 155 GB, kein Auth, ~24h auf the HPC cluster
sbatch --account=your_slurm_account --time=24:00:00 --wrap="
OUT=/beegfs/u/<user>/<group>/shared/datasets/pacbio_kinnex_pbmc
mkdir -p \$OUT && cd \$OUT
BASE=https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p-rep1
for f in 2-DeduplicatedReads/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam{,.pbi} 2-DeduplicatedReads/scisoseq.mapped.bam{,.bai}; do
  wget -c \$BASE/\$f
done
"
```

### 3.3 ArrayExpress Tonsil (E-MTAB-13687)

```bash
# Aria2c parallel, ~300-500 GB raw fastq, browse first to pick BCR-modality only
# https://www.ebi.ac.uk/biostudies/files/E-MTAB-13687/
# Idealerweise nur scBCR-fastqs (BCR Library) — meist ~10-50 GB statt full multimodal
```

### 3.4 EGA (EGAS00001006779) — Controlled

```bash
# 1. DAC-Antrag via https://ega-archive.org/datasets/EGAD00001009998
#    Begründung: IGHG4-ChimDup Stummheits-Test, citation Boahen+Belios+pseudocaller
# 2. nach Approval (Wochen): pyega3 download mit credentials
pyega3 -c 4 fetch EGAD00001009998 --output-dir /beegfs/u/<user>/datasets/ega_pathogen_pbmc
```

### 3.5 GSA-Human (HRA003750) — Controlled

- Antrag bei West China Rheumatology DAC (HDAC001850) via https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003750
- Nach Approval: HTTPS-Links für direct download
- Strategie: Antrag erst stellen wenn Lokal-Pipeline für BD-Rhapsody-3'-tag validiert ist

---

## 4. Lokal-vs-Remote-Inventar-Übersicht

| Quelle | Status auf the HPC cluster | Reproduzier-Aufwand |
|--------|-------------------|---------------------|
| PBMC BioIVT PacBio (intern) | **READY** 13 GB | analysis_complete (12 G4-cells) |
| HPRC LCL Iso-Seq | **READY** 363 GB | analysis_complete (klonal-Bias-Verständnis) |
| **PRJNA922682 FLAIRR** | **MUSS HERUNTERGELADEN** | 1× SLURM job, ~12h, 15 GB output |
| **PacBio Kinnex PBMC** | **MUSS HERUNTERGELADEN** | 1× SLURM wget, ~24h, 155 GB per rep |
| **E-MTAB-13687 Tonsil BCR** | **MUSS HERUNTERGELADEN** (BCR-only subset) | ~30 GB, 24h |
| **EGAS00001006779 Pathogen** | **DAC-ANTRAG NÖTIG** | 2-6 Wochen + 50 GB |
| **HRA003750 IgG4-RD** | **DAC-ANTRAG NÖTIG** | Wochen + ~100 GB |
| DRG-IGHG4 medRxiv 2024 | Data-Availability nicht öffentlich extrahiert | Korrespondenz-Autor anschreiben |
| Mausembryos (lokal) | Download-Stubs vorhanden, **irrelevant** | — |
| NRW Kubisch-Drop | DNA, nicht RNA | irrelevant |

---

## 5. Empfehlungs-Reihenfolge für ChimDup-Silence-Test

### Phase 1 — sofort (kein DAC nötig)
1. **PRJNA922682 (FLAIRR)** herunterladen → 10 healthy PBMC donors mit **CH1-3-Sequenz** voller Länge bei Q60. Mapping auf chimDup-vs-canonical-Referenz mit LLmap, per-read genotyping.
2. **PacBio Kinnex Public PBMC** (1-2 reps) herunterladen → single-cell mit Kinnex (ähnlich BioIVT aber Revio). Zellzahl-Statistik für 12-Cell-Limit bestätigen.

### Phase 2 — DAC-Anträge parallel laufen lassen
3. **EGAS00001006779** Antrag bei EGA (Boahen DAC) — Begründung Cross-Donor-Validation, sobald genehmigt analysieren.
4. **HRA003750** Antrag bei West China — strategisch für ChimDup-Selektions-Test in IgG4-RD-Kohorte (nicht für direkten Expression-Test).

### Phase 3 — Eigene Datengewinnung (langfristig)
5. **Tonsillen-Iso-Seq** kollaborativ — keine bestehende public Source hat das: Targeted-IGH-Iso-Seq auf Mandeln. Mit Mount Sinai (FLAIRR-Gruppe) oder PacBio Kollaboration sondieren. Belios-Paper-Strategie.

---

## 6. Open Questions / Risiken

- **FLAIRR donor-IDs:** Sind die 10 healthy PBMCs polyklonal genug? Paper liest sich so, aber finale Validation erst nach Download.
- **Kinnex public IGHG4-coverage:** 20k cells, ~1-2% B-cells × 10% class-switched × <1% IGHG4 → ~20-40 cells expected, vergleichbar zu unserem BioIVT (12 cells). Reicht für Replikation, nicht für robuste ChimDup-Stratifikation.
- **EGA DAC-Wartezeit:** Erfahrungswerte 2-6 Wochen, manchmal Monate. Phase 1 nicht blockieren lassen.
- **CH1-SNP-Diskriminierung in short-read Kohorten** (Tonsil-Atlas, IgG4-RD): braucht custom probe design oder targeted re-sequencing — nicht naiv aus 10x-Daten extrahierbar.
- **Reproduzierbarkeit:** Belios-Paper sollte den Lokal-PBMC-Befund (12 cells) plus FLAIRR-Validation enthalten, nicht stark auf single-donor anlehnen.

---

## Referenzen

- Ford et al. 2023, *J Immunol* — FLAIRR-seq method, BioProject PRJNA922682, DOI: 10.4049/jimmunol.2200825
- Boahen et al. 2024, *iScience* — Pathogen-stimulated PBMC long-read, EGA: EGAS00001006779, DOI: 10.1016/j.isci.2024.110525
- Massoni-Badosa et al. 2024, *Immunity* — Tonsil Atlas E-MTAB-13687, DOI: 10.1016/j.immuni.2024.01.006
- Lu et al. 2023, *JCI Insight* — IgG4-RD scRNA, GSA-Human HRA003750, DOI: 10.1172/jci.insight.167602
- Mancilla Moreno et al. 2024 medRxiv — IGHG4 DRG, DOI: 10.1101/2024.11.12.24317004
- Memory cross-refs: `pbmc_paralog_finding`, `hprc_lcl_clonality_blocks_ighg4_test`, `nrw_s3_dataset`, `hprc_isoseq_data`, `ighg4_dose_effect_transcriptome`, `mikocziova_2025_ighg4`
