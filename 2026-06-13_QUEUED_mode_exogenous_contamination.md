# QUEUED — Mode "Exogene Kontamination & Artefakt-Provenienz" (Operator, 2026-06-13)

Start NACH Transcript-Mode-e2e-Grün. Operator-Frage: "Das sollten wir gezielt
abbilden und mitbedenken. Wie ginge das am besten?" → Architektur-Empfehlung unten.

## Operator-Taxonomie (vollständig, mit Referenzen — NICHT verlieren)

Gemeinsamer Nenner ALLER Punkte: sie **täuschen niederfrequente „Varianten"/Signale
vor**, die in Wahrheit durch einen **alternativen Read-Ursprung** oder einen **bekannten
Chemie-/Biologie-Prozess** erklärt sind.

**A. Exogene Biologie**
- EBV/LCL-Artefakte: Großteil der Referenz-Samples (1000G, HapMap, teils GIAB) aus
  EBV-immortalisierten LCLs → EBV-Episomen/-Integrationen + LCL-erworbene somatische
  Mutationen & CNVs als niederfrequente „Varianten".
- Mikrobiom & Reagenzien: Hautkeime (Cutibacterium, Staphylococcus); Mycoplasma
  (systematisch in Zelllinien-RNA-seq, Olarerin-George/Hogenesch GEO-Screen); „Kitome"
  (Salter et al.: Bradyrhizobium/Ralstonia, v.a. low-biomass); Stuhl/Saliva → Nahrungs-/
  Pflanzen-DNA.
- Index-Hopping/Barcode-Swapping: Patterned-Flowcells (ExAmp, NovaSeq) ~0,1–2 % Reads
  im falschen Sample → echte Varianten aus Nachbarproben „bluten ein".
- Spike-ins: PhiX (oft nicht rausgerechnet), ONT DCS/Lambda, ERCC (RNA).

**B. Damage-/Chemie-Artefakte (täuschen Low-VAF-SNVs vor; strangverzerrt)**
- 8-oxoG: G→T durch akustisches Shearing (Costello/Chen); Strandbias, niedrige VAF;
  GATK Orientation-Bias-/FoxoG-Filter.
- FFPE-Deaminierung: C→T / G→A, niederfrequent, strangverzerrt.
- Cytosin-Deaminierung an Read-Enden (ancient/degraded DNA), C→T-Profil.
- RNA-Editing: A→I (als A→G, v.a. Alus durch ADAR), C→U (APOBEC) — sieht aus wie SNV.

**C. Mapping-/Referenz-Artefakte (Segdup-Terrain)**
- NUMTs: nukleäre mtDNA-Segmente → falsche „Heteroplasmie" niederfrequent.
- Paralog-/Segdup-Misalignment: hochidentische Paraloge falsch platziert → falsche
  Low-VAF-SNVs (IGHG4, Retina-Segdups).
- Prozessierte Pseudogene: GBAP1 vs GBA1 (Thesis-relevant!), PMS2/PMS2CL, SMN1/2 — im
  Transkriptom intronlose „Expression".
- Decoy/unplaced contigs, rDNA, Satelliten: hs37d5-Decoy, kollabierte Dups in der Ref.

**D. Chimäre & Concatemere**
- Ligations-/PCR-Chimären → falsche SVs/Fusionen.
- ONT: ungetrimmte Adapter → chimäres Mapping; Concatemere/Fused-Reads (zwei Moleküle/Pore).
- RT-Template-Switching → artifizielle trans-Spleiß-/Fusion-cDNA (Single-Cell TSO).

**E. Transkriptom-Spezialitäten**
- Internal Priming auf genomischen Poly-A (falsche 3'-Enden); gDNA-Kontamination
  (intronische/intergene Reads); unvollständige rRNA-Depletion; Globin-Reads (Vollblut);
  circRNA-Backsplice (sieht aus wie Chimäre); Readthrough/conjoined genes (SIDT2-TAGLN);
  Ambient-RNA-„Soup" (scRNA).

**F. Klinisch heikel**
- Klonale Hämatopoese (CHIP), somatischer Mosaizismus im Blut → täuschen germline bei
  niedriger VAF vor. Tumor-in-Normal-Kontamination.

## Architektur-Empfehlung (Agent 1)

**KEIN monolithischer „Kontaminations-Mode".** Stattdessen eine **Provenienz-/
Erklärungs-Schicht**, orthogonal über alle Modi — weil all diese Artefakte aus LLmap-Sicht
EINE Struktur teilen: *eine konkurrierende Erklärung für die Read-Herkunft/-Chemie, die
besser ist als „Host-Variante".* Genau dafür ist lossless WaveCollapse gebaut: statt
argmax-Collapse auf die Host-Referenz + Variant-Call, hält LLmap **Wahrscheinlichkeitsmasse
über konkurrierende Hypothesen** und labelt die dominierende — das Artefakt ist dann keine
Variante, sondern Masse, die einem anderen Ursprung/Prozess gehört.

**Drei Mechanismus-Klassen (je auf vorhandene LLmap-Bausteine gemappt):**

1. **Konkurrierende-Referenz-Hypothesen** (A exo, C mapping):
   - Exo-Katalog als zusätzliche Kandidaten-Quelle in der WaveCollapse-Likelihood
     (EBV, Mycoplasma, Kitome-Taxa, PhiX/ERCC/Lambda, Decoy). Kataloge existieren teils
     (SILVA/UNITE/BOLD/PR2, Viren, Bakterien).
   - NUMT/Paralog/Pseudogen/Segdup = **LLmaps Kernkompetenz schon heute** (PSV, segdup,
     igh_locus). NUMT = mtDNA↔Kern-Hypothese; Pseudogen (GBAP1/GBA1) = intronlos vs
     intronhaltig (Transcript-Mode kennt Exon-Struktur!).

2. **Kontext-konditionierte Fehler-/Editier-Modelle** (B damage/editing):
   - Ein Mismatch in einem Damage-/Editier-KONTEXT wird als *Prozess* abgewertet, nicht
     als Variante: 8-oxoG (G→T, Strangbias/FoxoG), FFPE (C→T strang), ancient (C→T
     Read-Ende), A→I (A→G in Alu/ADAR), C→U (APOBEC). LLmap hat dafür schon `rnamod`
     (m6A/A-to-I, AID-Footprint) — das ist der Keim. Erweitern zu einem allgemeinen
     **kontext+strang-konditionierten Substitutions-Prior**.

3. **Struktur-/Split-Hypothesen** (D chimera, E circRNA/readthrough):
   - Chimera-Modul (Block 7) + VDJ-Mask existieren; circRNA-Backsplice kennt
     `schema_transcript` schon (XK=MAPPED_CIRCULAR). Adapter/Concatemer = Split-Read +
     Adapter-Trim. Index-Hopping = Barcode-/Cross-Sample-Signal (single-cell-Tags da).

**Output-Doktrin (= Lossless auf Artefakte angewandt):** kein Read/Call wird still
gedroppt; jeder bekommt einen **Provenienz-Tag mit Konfidenz**:
`host-variant | paralog-collapse | numt | exo:EBV | exo:mycoplasma | spikein:phix |
damage:8oxoG | edit:A2I | chimera | index-hop | pseudogene:GBAP1 | chip-somatic`.
Ein Downstream-Variant-Caller filtert/flaggt dann nach Provenienz statt blind Low-VAF-SNVs
zu callen. Das ist exakt die Transcript-Mode-Doktrin (labeln statt verwerfen), generalisiert.

**Verifikations-Achsen** (analog mode_architecture.md): Host-Mismatch + Exo-Match +
Coverage-Inseln + Strang-Bias + Kontext-Signatur (Damage/Editing) + GC/Codon-Fremdheit +
Barcode-Kontext + Split-Geometrie.

**Warum LLmap das besser kann als kraken2/GATK-Filter:** die machen je EINEN Artefakt-Typ
isoliert + argmax. LLmap vereint sie in EINER lossless Likelihood mit Konfidenz — und löst
die gefährlichste Klasse (Paralog/NUMT/Pseudogen-Collapse) ohnehin schon als Kern.

→ Bei Start: Design-Doc `docs/design/llmap_provenance_axes.md` (analog mode_architecture),
  Operator-Review, dann inkrementell — zuerst die Klassen, die schon Bausteine haben
  (Paralog/NUMT/Pseudogen via PSV+Transcript; Damage/Editing via rnamod-Erweiterung).

## Output-Design (Operator-Präferenz, 2026-06-13): separater quantifizierbarer Bucket

Operator: "separat in einem Bucket ausgeben, mit den inhaltlichen Flags versehen,
dadurch quantifizierbar."

Konkrete Form:
- **Separater Output-Kanal** (nicht in den Haupt-BAM gemischt): ein „Provenienz-Bucket".
  - `<out>.provenance.bam` (oder ein eigenes Bucket im BucketPyramid) — alle als
    Artefakt/Kontamination/Mismapping geflaggten Reads landen hier, NICHT still gedroppt.
  - `<out>.provenance.parquet` — eine Zeile pro geflaggtem Read:
    `read_id, provenance_class, sub_source, confidence, evidence_axes[], host_alt_pos,
     exo_ref, strand_bias, context_motif`.
- **Inhaltliche Flags** = die Provenienz-Klassen (host-variant/paralog/numt/exo:*/
  spikein:*/damage:*/edit:*/chimera/index-hop/pseudogene:*/chip-somatic) als BAM-Tag
  (z.B. `PV:Z:exo:mycoplasma` + `PQ:f:<conf>`).
- **Quantifizierung** = Summary-Tabelle `<out>.provenance_summary.tsv`:
  pro Klasse `count`, `fraction_of_total`, `mean_confidence`, `bases`. So wird
  „0,8 % PhiX / 1,2 % Index-Hop / 2,1 % NUMT-Collapse / 0,3 % Mycoplasma" direkt ablesbar.
- **Generalisiert `lossless_aggregator`**: das macht heute schon Invarianten-Accounting
  (kein Read verloren) — die Provenienz-Quantifizierung ist die inhaltliche Erweiterung
  desselben Accountants (jeder Read MUSS in genau einem Bucket landen: host-call ODER
  provenance-class; Summe == Input = beweisbar lossless).

Vorteil für die Klinik/Thesis: der Variant-Caller sieht nur den sauberen Host-Bucket;
der Provenienz-Bucket ist separat auswertbar (Kontaminations-QC + biologisch echte
Sonderfälle wie A→I-Editing oder CHIP getrennt quantifiziert, nie als germline-SNV vermengt).

## Realdaten-Validierung (Operator, 2026-06-13) — NACH Mode-Fertigstellung
Operator: „1-2 Longread-Genome + Short-read-Genome auf diese Special Buckets testen."
- **Input (auf the HPC cluster/BeeGFS):** Longread = HG002 HiFi + HG002/ONT; Short-read = HG002 /
  NA12878 Illumina WGS. Ein **LCL-stämmiges** Sample (1000G/HapMap) ist der beste Stresstest
  (sollte messbares EBV-`exo` + NUMT-`para` + Damage-`dmg`-Spektrum zeigen).
- **Lauf:** `llmap align` mit always-on Provenienz-Layer → pro Sample `contamination_spectrum.parquet`.
- **Erwartete Readouts pro Bucket-Familie:** `host`-Fraktion (Großteil), `exo:*` (EBV bei LCL,
  PhiX-Spike-in-Reste), `para`/`numt`/`pseudo` (Segdup/mt-Confounds), `dmg:*`/`edit:*` (8-oxoG-/
  Deaminierungs-/ADAR-Fraktion — Long- vs Short-read unterscheiden sich hier stark), `chim`/`xsample`.
- **Plattform-Vergleich:** Longread vs Short-read zeigt unterschiedliche Artefakt-Profile
  (z.B. ONT-Concatemere/`chim` vs Illumina-8-oxoG/Index-Hop) → validiert, dass die Detektor-Klassen
  plattform-spezifisch korrekt greifen.
- **Σ-Invariante als Lossless-Beweis:** Σ(alle PV-Klassen) == N_input pro Sample.
- Ausführung = the HPC cluster (Agent 2's Daten-Lane); Detektor-Klassen = Agent 1.
