# Per-Position Splice-Determinism + Boundary-Konkordanz-Testset (Operator-Idee, 2026-06-13)

Inspiriert von AlphaGenome-Talk (ESHG): Transkript-Predictions trafen Exon-Grenzen gut.
Operator-Frage: wie "blurry" sind diese Grenzen, wenn man ALLE Read-Mappings auf Coverage-
Ebene überlagert (verschiedene Isoformen → verschiedene Grenzen)? → Wunsch nach einer
**per-Position-Metrik**, die angibt, wie eindeutig eine Stelle gemappt/gespliced ist.

## Die Metrik: D(pos) ∈ [0,1] — Splice-Determinismus / Isoform-Konsens

Pro genomischer Position `pos`: überlagere alle spliced Reads, die `pos` überdecken,
und klassifiziere ihren lokalen Splice-Zustand an `pos`:
- **included**  — Read aligned `pos` als Exon-Körper (M an der Stelle)
- **excluded**  — Read überspringt `pos` innerhalb eines Introns (N überdeckt `pos`)
- **boundary**  — `pos` ist exakt eine Donor/Acceptor-Kante; bp-genaue Übereinstimmung zählt

`D(pos)` = Anteil der informativen Reads, der dem **Modus** (häufigster Zustand) entspricht:
- konstitutives Exon, alle inkludiert → D = 1.00
- in 9/10 excluded → Konsens "excluded", D = 0.90
- bp-genaue Grenzen teils verschoben → D fällt weiter (z.B. 0.76)
- 5/5 inkludiert/excluded → maximal ambig, D = 0.50

Äquivalent: `D(pos) = 1 − H_norm(pos)`, mit H = normalisierte Entropie der lokalen
Splice-Zustände. = per-Base verallgemeinertes PSI, das Inklusion/Exklusion UND
Grenz-Präzision in einen Determinismus-Score zusammenfasst.

**Warum das exakt zu LLmap passt:** WaveCollapse hält pro Read Wahrscheinlichkeitsmasse
über Platzierungen/Splice-Zustände statt argmax. `D(pos)` ist genau die **per-Position
aggregierte Collapse-Konfidenz**. Lossless heißt: wir collapsen die Isoform-Diversität NICHT
weg, sondern quantifizieren sie. Das ist ein natürlicher Output des Transcript-Mode (wir
haben jetzt spliced Alignments + Junctions).

## Output
- Per-Base-Track `<out>.splice_determinism.bedgraph` (D je Position).
- Per-Junction/Exon-Tabelle: PSI, Donor/Acceptor-bp-Streuung, n_reads, modal_state.
- Genom-weite Verteilung + Auflösungs-Report (wie scharf sind die Grenzen?).

## Variant-Relevanz (Operator-Frage: sinnvoll?) — JA
Die Interpretation einer Variante ist **konditioniert auf D(pos)**:
- Variante in konstitutivem Exon (D≈1) = anderer Impact-Prior als in nativ alternativ
  gesplictem Bereich (D≈0.76).
- Eine Splice-affecting/sQTL-Variante, die D *senkt* (Exon-Skipping induziert), ist direkt
  als Effekt messbar (ΔD Sample vs. Baseline).
- noncoding eQTL/sQTL: D + lokale Junction-Usage sind die natürlichen Readouts.

## Testset (Operator-Wunsch, auf the HPC cluster nach Snapping-Validierung)
1. Echte iso-seq-FLNC gegen GRCh38 + **GENCODE**-Annotation mappen (`--mode transcript`).
2. **Boundary-Konkordanz:** gemappte Donor/Acceptor vs. GENCODE-Exon-Grenzen (bp-genau +
   ±k-Toleranz-Jaccard) — wie gut treffen die Mappings die DB-Grenzen?
3. **Transkript-Verteilung:** welche Isoformen/Transkripte, in welchem Verhältnis.
4. **Auflösung:** `D(pos)`-Track + Grenz-Streuung-Histogramm.
Vergleichbar zur AlphaGenome-Beobachtung (treffen unsere Mappings die Grenzen so gut?).

## GTEx / eQTL / sQTL — zu früh für LLmap?
**Gestaffelt — die Metrik JA jetzt, die Prädiktion mittelfristig:**
- **Jetzt (natürliche Transcript-Mode-Erweiterung):** `D(pos)` + PSI + Boundary-Konkordanz
  pro Sample aus den eigenen Reads. Kein neues großes Modul, nutzt vorhandene spliced
  Alignments. Das ist der richtige nächste Schritt, sobald Transcript-Mode solide ist.
- **Mittelfristig (eigenes Phasen-Projekt):** GTEx-abgeleitete **Population-Baseline** pro
  Gen/Variante (erwartete D/PSI-Verteilung) selbst aus bestehenden Daten rechnen →
  Sample-vs-Baseline-ΔD für Impact-Schätzung; dann eQTL/sQTL-Prädiktion. Das braucht
  GTEx-Ingestion + statistisches Modell + Validierung → eigener Block, NICHT jetzt.
- **Empfehlung:** das Per-Sample-`D(pos)`-Schema jetzt so designen, dass es später direkt
  gegen eine Population-Baseline vergleichbar ist (gleiche Koordinaten/Einheiten). Dann ist
  der Sprung zu GTEx/QTL inkrementell, nicht ein Rewrite.

**Kurzfazit:** Die Metrik ist sinnvoll, neuartig, LLmap-nativ und JETZT machbar. Die
QTL-Prädiktion ist die richtige mittelfristige Ambition, aber ein eigener Block.
