# Provenienz-Klassen — was wir NOCH NICHT bedacht haben (Brainstorm, 2026-06-13)

Bisher: A exogen · B Damage/Editing · C Mapping-Confusion · D Chimär/Concatemer ·
E Transkriptom · F klinisch (CHIP/Mosaik/Tumor-in-Normal). Hier die fehlenden Familien.
Jede: **Signatur → wie sie eine „Variante" vortäuscht → LLmap-Detektions-Hook.**

## G. Organismale Chimärität / Cross-Individuum (NEU — klinisch hochrelevant)
- **Maternal-fetale Mikrochimärie** (cffDNA in Mutter / maternale Zellen im Kind): Reads aus
  einem ZWEITEN Individuum-Genom → fake low-VAF-Varianten (echte Allele eines anderen Menschen).
- **Transplant-/Transfusions-Chimärität:** Spender-DNA im Empfänger (Knochenmark-Tx → Blut ist
  Spender-Genotyp!). Riesige Falle für „germline"-Calls aus Blut nach KMT.
- **Twin-/Blut-Chimärismus**, vanishing twin.
  → Hook: Reads tragen ein **konsistentes alternatives Genotyp-Set** (anderes Individuum), nicht
    zufällige Fehler — erkennbar an gehäuften, ko-segregierenden Nicht-Host-Allelen über Loci.
    `xindiv:maternal` / `xindiv:donor`. (Unterscheidet sich von Index-Hopping = technisch.)

## H. Immunlocus-somatische Diversität (NEU — direkt IGH/TCR, dein Kernthema)
- **V(D)J-Rekombination:** IGH/IGK/IGL/TCR-Reads sind somatisch umgelagert → dichte „Varianten"/
  Breakpoints, die KEINE germline sind, sondern echte Immun-Diversität. LLmap hat schon die
  VDJ-Mask (Block 7) + junction_hunter.
- **Somatische Hypermutation (SHM):** AID-getrieben, hohe Punktmutationsdichte in V-Regionen →
  sieht aus wie Mutations-Hotspot. LLmap hat `aid_footprint` (rnamod)!
- **Klassenwechsel-Rekombination (CSR):** S-Region-Junctions (das IGHG4-Thema).
  → Hook: `vdj:recomb` / `vdj:shm` / `vdj:csr` — endogene Biologie, flaggen statt als germline callen.

## I. Referenz-seitige Confounder (NEU — die Referenz SELBST als Fehlerquelle)
- **Referenz-Fehler/Assembly-Lücken:** falsche Base in GRCh38 → jeder Read „variant" dort.
- **Kollabierte Dups in der Referenz** (die Referenz hat 1 Kopie wo 2 existieren) → alle Paralog-
  Reads stapeln sich → fake Heterozygotie. (hs37d5-Decoy-Thema, aber breiter.)
- **Referenz-Allel = Minor-Allel:** an manchen Loci ist die Referenz das seltene Allel → „Variante"
  ist eigentlich das häufige Allel. Kontext für Interpretation.
- **Fehlende Insertion in der Referenz:** Reads über eine nicht-referenzierte Insertion werden
  soft-geclippt/falsch aligned → fake SNVs an den Clip-Rändern.
- **Liftover-Artefakte** (GRCh37↔38↔T2T): Koordinaten-/Allel-Flips.
  → Hook: `ref:error` / `ref:collapsed` / `ref:minor` — braucht eine Referenz-Qualitäts-Annotation
    (T2T/Pangenome als Wahrheit; LLmap hat Pangenome-GAF-Bridge).

## J. Multiplizität / Amplifikation (NEU — VAF-Inflation)
- **PCR-Duplikate / optische Duplikate / PCR-Jackpotting:** dieselbe Molekül-Kopie n-fach →
  künstlich hohe/niedrige VAF. (Normalerweise vor-gefiltert, aber oft unvollständig.)
- **WGA-Artefakte (MDA, Single-Cell):** Strang-Displacement-Chimären + Allelic Dropout +
  Amplifikations-Bias → fake Allel-Imbalance, Chimären.
- **UMI-/Barcode-Kollisionen** (Single-Cell): zwei Moleküle teilen UMI → fälschlich gemerged.
  → Hook: `dup:pcr` / `dup:optical` / `wga:chimera` / `umi:collision`. LLmap hat single-cell-Tags.

## K. Mobiles/strukturelles Genom (teils via junction_hunter abgedeckt)
- **Somatische Mobile-Element-Insertionen (L1/Alu-Retrotransposition):** novel Insertion → Reads
  sehen chimär/geclippt aus → fake SV. Aktiv in Tumor + Gehirn.
- **ecDNA / eccDNA (extrachromosomale zirkuläre DNA):** Onkogen-Amplikons auf Zirkeln →
  zirkuläre Junctions, hohe Kopie → fake CNV/Fusion.
- **Genkonversion zwischen Paralogen:** nicht-reziproker Transfer → Paralog-Allel erscheint als
  „Variante" im Ziel-Paralog. Subtiler als Misalignment (C) — echte Sequenz, falscher Ursprung.
  → Hook: `mei:L1`/`mei:alu` / `ecdna` / `geneconv`. NAHR ist schon in junction_hunter (Mode-5).

## L. Reagenz-/Synthetik-Sequenzen (erweitert A über PhiX hinaus)
- **Klonierungs-Vektoren** (pUC, lambda, BAC-Backbone), **E. coli** aus rekombinanten Enzymen/
  Plasmid-Prep, **Adapter-/Primer-Dimere**, **Oligo-/gBlock-Kontamination**, synthetische Spike-in-
  Controls (Sequins, SIRVs in RNA).
  → Hook: `synth:vector` / `synth:ecoli` / `adapter-dimer`. Panel-basiert (wie exo).

## M. Echte Low-VAF-Biologie (KEIN Artefakt — aber muss von Artefakt getrennt sein)
- **Echte mtDNA-Heteroplasmie** (vs NUMT-Artefakt C!) — der biologische Zwilling der NUMT-Falle;
  beide low-VAF auf mt — MÜSSEN getrennt werden. `mt:heteroplasmy` (real) vs `numt` (Artefakt).
- **Entwicklungs-Mosaizismus** (post-zygotisch, gewebespezifisch) — echt, aber low-VAF.
- **Telomer-/Zentromer-/Satelliten-Array-Reads:** nicht eindeutig platzierbar (HOR-Arrays,
  α-Satellit) → MAPQ0-Stapel. `repeat:satellite` (teils C-rDNA).

## Querschnitt-Hook: Strangbias + Ko-Segregation als universelle Diskriminatoren
Viele dieser Klassen teilen Signaturen, die LLmap schon hat/leicht aggregiert:
- **Strangbias** (Damage, ref-error) — site-level.
- **Ko-Segregation** (xindiv, chimerism, gene-conversion): mehrere „Varianten" auf DENSELBEN Reads
  → ein zweites konsistentes Haplotyp-Set, nicht zufällige Fehler. Das ist der stärkste
  Diskriminator für „anderes Individuum/Paralog/Kopie" — und passt exakt zu WaveCollapse
  (Masse kollabiert auf ein kohärentes alternatives Set).

## Konsolidierung mit Agent 2 (unabhängig gebrainstormt → Superset)
Agent 2's ⭐-Ergänzungen, die ich NICHT hatte (allesamt must-add):
- ⭐ **Xenograft/PDX-Maus-Kontamination** — 10-50% Maus-Reads in humanen Tumor-Xenografts.
  Riesig in der Krebsgenomik. `exo:mouse`. (Mein größter Miss.)
- ⭐ **Zelllinien-Kreuzkontamination / HeLa / fehl-ID'd Linien** (ATCC-Skandal) — ein zweites
  humanes Genom gibt sich als das Sample aus. `xindiv:cellline`.
- ⭐ **Bisulfit / enzymatische Methylierungs-Konversion** (C→T überall) — ohne Flag ist in
  Methylierungs-Seq JEDE Variante falsch. `B`-Evidenz, mode-level.
- ⭐ **Sample-Swap / Fehletikettierung** — das GANZE Sample ist die falsche Person
  (Genotyp-Konkordanz/Sex/Ancestry). `xsample:swap`.
- ⭐ **Ancestry-Referenz-Mismatch** — nicht-europäisches Genom vs GRCh38 → systematischer
  Ref-Allel-Bias an ancestry-divergenten Stellen. `ref:ancestry`.
- **HERV-Integrationen, AAV/Lenti/CRISPR-Transgen-Vektor** (Gene-Therapy), **PCR-Jackpot**
  (HIGH-VAF, nicht low!), **Plattform-Kontext-Fehler** (ONT/PacBio-Homopolymer, Illumina-GGC).

Meine Uniques, die Agent 2 nicht betonte:
- **Ko-Segregation als universeller Diskriminator** (s.u. Querschnitt) — architektonisch das Wichtigste.
- **Liftover-Artefakte** (GRCh37↔38↔T2T Koordinaten-/Allel-Flips). `ref:liftover`.
- **UMI-/Barcode-Kollisionen** (Single-Cell). `umi:collision`.
- **Telomer/Zentromer/α-Satelliten-HOR-Array-Reads** (MAPQ0-Stapel). `repeat:satellite`.
- **Referenz-Allel = Minor-Allel** (Interpretations-Kontext). `ref:minor`.

## Enum-Design-Empfehlung (Antwort auf „Enum erweitern?")
**NICHT alle ~40 Spezifika ins Enum hardcoden.** Stattdessen offen-endig wie `ModificationKind`
(Block 2.6: uint16 + `NovelUnclassified` + free-form `CustomKindTag`):
- **Geschlossenes Enum nur für die ~9 Mechanismus-FAMILIEN:** host · exo · paralog_confound ·
  base_artifact · structural · cross_individual · reference_artifact · immune_somatic · amplification.
- **Free-form Sub-Detail-String** für die Spezifika (`exo:mouse`, `exo:ebv`, `ref:ancestry`,
  `dmg:8oxoG`, `vdj:shm`, …) — wächst ohne Enum-Churn/ABI-Bruch.
So bleibt die Σ-Invariante sauber (jeder Read → genau eine Familie) UND erweiterbar.

## Priorisierung (höchster Wert für LLmaps Kontext)
1. **H VDJ/SHM/CSR** — direkt dein IGH-Kern, Bausteine da (VDJ-Mask, aid_footprint, junction_hunter).
2. **C↔M NUMT-Artefakt vs echte mt-Heteroplasmie trennen** — Brot-und-Butter-Falle, biologisch wichtig.
3. **G Cross-Individuum (Transplant/maternal)** — klinisch dramatisch (KMT-Blut = Spender-Genotyp).
4. **I Referenz-Confounder** — die Referenz selbst; T2T/Pangenome als Wahrheit (Bridge da).
5. **J Duplikate/WGA/UMI** — VAF-Inflation, single-cell-Tags da.
6. **L Vektor/E.coli/Dimer** — billig panel-basiert mit A zusammen.
