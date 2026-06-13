# Mode "Exogene biologische Kontamination" — Design-Input (Operator, 2026-06-13)

Roter Faden über ALLE Punkte: **Quellen, die niederfrequente „Varianten" vortäuschen**
(fake low-VAF SNV/SV/Heteroplasmie). Ziel des Modes: solche Reads/Sites als das markieren,
was sie sind — **flaggen, nie still filtern** (gleiche Lossless-Philosophie wie Transcript-Mode).

## Operator-Taxonomie (vollständig)

### A. Exogene Sequenz (mappt auf NICHT-Host)
- **EBV / LCL-Artefakte:** 1000G/HapMap/teils GIAB aus EBV-immortalisierten LCLs →
  EBV-Episomen/-Integrationen + LCL-erworbene somatische Mutationen & CNVs als low-freq „Varianten".
- **Mikrobiom & Reagenzien:** Hautkeime (Cutibacterium, Staphylococcus), Mycoplasma (systematisch
  in Zelllinien-RNA-seq, Olarerin-George/Hogenesch GEO-Screen), „Kitome" (Salter et al.:
  Bradyrhizobium/Ralstonia, v.a. low-biomass), Nahrungs-/Pflanzen-DNA (Stuhl/Saliva).
- **Index Hopping / Barcode Swapping:** Patterned-Flowcells (ExAmp/NovaSeq) ~0.1–2% Reads ins
  falsche Sample → echte Nachbarproben-Varianten „bluten ein". Paradigmatische Fake-low-VAF-Quelle.
- **Spike-ins:** PhiX (oft nicht rausgerechnet), ONT DCS/Lambda, ERCC (RNA).

### B. Damage- & Chemie-Artefakte (täuschen Low-VAF-SNVs vor; strangverzerrt)
- **8-oxoG:** G→T durch akustisches Shearing (Costello/Chen), Strandbias → GATK Orientation-Bias/FoxoG.
- **FFPE-Deaminierung:** C→T / G→A, low-freq, strangverzerrt.
- **Cytosin-Deaminierung an Read-Enden:** ancient/degraded DNA, C→T-Profil.
- **RNA-Editing:** A→I (liest A→G, ADAR in Alus), C→U (APOBEC). Sieht aus wie SNV, ist Edierung.

### C. Mapping-/Referenz-Artefakte (Segdup-Terrain — LLmap-Kern)
- **NUMTs:** nukleäre mtDNA-Segmente → falsche „Heteroplasmie" bei niedriger Frequenz.
- **Paralog-/Segdup-Misalignment:** hochidentische Paraloge falsch platziert (IGHG4, Retina-Segdups).
- **Prozessierte Pseudogene:** GBAP1/GBA1 (Thesis-relevant), PMS2/PMS2CL, SMN1/2. Im Transkriptom
  intronlose „Expression".
- **Decoy/unplaced contigs, rDNA, Satelliten:** hs37d5-Decoy, kollabierte Dups in der Referenz selbst.

### D. Chimäre & Concatemere
- Ligations-/PCR-Chimären → falsche SVs/Fusionen.
- ONT: nicht getrimmte Adapter → chimäres Mapping; Concatemere/Fused-Reads (zwei Moleküle/Pore).
- RT-Template-Switching → artifizielle trans-Spleiß-/Fusion-cDNA (TSO, v.a. Single-Cell).

### E. Transkriptom-Spezialitäten
- Internal Priming auf genomischen Poly-A (falsche 3'-Enden), gDNA-Kontamination (intronische/
  intergene Reads), unvollständige rRNA-Depletion, Globin-Reads (Vollblut), circRNA-Backsplice
  (sieht aus wie Chimäre), Readthrough/conjoined genes (SIDT2-TAGLN), Ambient-RNA-„Soup" (scRNA).

### F. Klinisch heikel (Cross-Sample / somatisch)
- Klonale Hämatopoese (CHIP) & somatischer Mosaizismus im Blut → täuschen germline bei low VAF vor.
- Tumor-in-Normal-Kontamination.

## Architektur-Meinung (Agent 2) — siehe Room-Diskussion
Kein Monolith. 5 mechanistisch verschiedene Detektor-Klassen, alle als **Evidenz-Kanäle in die
bestehende Multi-Signal-Fusion-Engine (Block 4.5)**; Output = per-Read + per-Site Artefakt-
Posterioren als Tags (flag, nie filter). Wiederverwenden: PSV/WaveCollapse (C), Mode-6 taxbin (A),
Chimera-Detektor/Block 7 (D). Neu & leichtgewichtig: Strandbias-/Damage-Profil-Detektor (B).
LLmap liefert EVIDENZ, nicht den klinischen germline/somatik-Call.
