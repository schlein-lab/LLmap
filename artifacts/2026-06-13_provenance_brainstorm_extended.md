# Provenance classes — extended brainstorm (Agent 2, 2026-06-13)

Beyond the v0.1 taxonomy (A exogen / B Damage+Editing / C structural / D cross-sample).
⭐ = aufnahmewürdig ins Enum jetzt. Mechanismus-Klasse in [Klammern].

## A — Konkurrierende Referenz (Origin-Bucket)
- ⭐ Xenograft/PDX-Maus-Kontamination (Tumor-Stroma, 10-50% Maus-Reads) [A]
- ⭐ Zelllinien-Kreuzkontamination / HeLa / fehl-ID'te Linien [A]
- ⭐ Chimärismus: Transfusion / Transplantation / materno-fetal + Mikrochimärismus (zweites humanes Genom) [A]
- Vektor/Transgen/Plasmid/Selektionsmarker (Gene-Therapy AAV/Lenti, CRISPR-Backbone) [A]
- Aktive Mobile Elemente (LINE-1/Alu/SVA) + HERV-Integrationen [A/C]
- ⭐ Echte mtDNA-Heteroplasmie (von NUMT TRENNEN, nicht beide `numt`) [A]
- Prophagen / integrierte Viren / ecDNA / eccDNA (zirkulär) [A/C]
- Wolbachia / Endosymbionten / Apicoplast / Chloroplast (Nicht-Mammalia-Samples) [A]
- Carryover/Bleed vom VORIGEN Run (Flowcell/Instrument-Memory) [A/D]
- Genomische QC-Spike-ins (NA12878-Mix), Lambda/T7/DCS-Controls [A]

## B — Base-Level-Artefakte (per-Site-Evidenz)
- ⭐ Bisulfit / enzymatische Methyl-Konversion (C→T überall) — Methyl-Seq [B]
- 5mC→T-Deaminierung an CpG (Signatur + Artefakt) [B]
- PCR-Jackpot / Early-Cycle-Polymerase-Fehler (HOHE scheinbare VAF) [B]
- Plattform-Kontext: ONT/PacBio-Homopolymer-Indels, Illumina-GGC/Dark-G, Dephasing [B]
- Mapping-/Ref-Allel-Bias an Indels → Alt-Allel-Dropout [B/E]
- GC-/Amplifikations-Bias (Coverage-Nonuniformität, MDA) [B]

## C — Strukturell (Split / Junction)
- ⭐ Referenz-Assembly-FEHLER (kollabierte Dups, Gaps, GRCh38-Issues) — eigene Klasse E [E]
- Rolling-Circle/MDA-Concatemere (Single-Cell-WGA) [C]
- Adapter-Dimer / No-Insert / Adapter-Read-Through [C]
- Repeat-Expansion / STR-Stutter (PCR-Slippage) [C/B]
- Soft-Clip-Artefakte an SV/Repeat-Grenzen [C/E]

## D — Cross-Sample / Population / Individuum
- ⭐ Sample-Swap / Fehletikettierung (ganzes Sample falsche Person; Genotyp/Sex/Ancestry-Check) [D]
- ⭐ Referenz-Ancestry-Mismatch (nicht-EUR vs GRCh38 → Ref-Allel-Bias) [D/E]
- PCR-/Optical-Duplikate (scheinbare VAF-Inflation) [D]
- Within-Run-Barcode-Kollision (≠ Hopping) [D]
- Reference-IS-minor-Allele-Sites (Referenz trägt seltenes Allel → alle „Varianten") [E]

## E — NEU: Referenz-Artefakte (die Referenz ist falsch, nicht der Read)
- Assembly-Kollaps/Gaps, Liftover-Fehler (GRCh37↔38), Alt-Contig-Mishandling, Decoy-Leakage,
  N-Basen/soft-masked-Regionen, chimäre Referenz (GRCh38 = Mosaik mehrerer Individuen),
  SD-CNV-Polymorphismus (Copy-Number variiert → Tiefe-Artefakt), PAR X/Y-Ambiguität.

## F — NEU: „Echte Biologie, die täuscht" (FLAGGEN, nicht buckten/filtern)
- ⭐ V(D)J-Rekombination + somatische Hypermutation + Class-Switch (IGH/TCR) — LLmap hat VDJ-Maske [F]
- Genkonversion (Paralog→Paralog echter Austausch) [F]
- Entwicklungs-/Gewebe-Mosaizismus, mosaik-CNV, LOY (Y-Verlust im Altersblut), Aneuploidie-Mosaik [F]
- Allel-spezifische Expression / Imprinting / X-Inaktivierung (RNA) [F]

## Schlüssel-Erkenntnis fürs Schema
Zwei neue Mechanismus-Klassen über A/B/C/D hinaus:
- **E (Referenz-Artefakt):** der Read ist korrekt, die Referenz falsch → braucht eigene Flag-Familie
  (`refartefact`), denn weder Origin-Bucket noch Base-Error noch Cross-Sample.
- **F (echte-Biologie-Confounder):** muss GEFLAGGT, aber NIE als Kontamination gebuckt/gefiltert werden
  (sonst strippt man Immunbiologie/Mosaizismus = echte Signale). Posterior + explizites `ambiguous`-Flag.
