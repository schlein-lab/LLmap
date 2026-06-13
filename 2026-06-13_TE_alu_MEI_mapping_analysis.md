# Transposable Elements / Alu / MEI im Mapping — ist das intellektuell gewürdigt? (2026-06-13)

Operator-Frage: sind Alu-/TE-Insertionen im Mapping intellektuell ausreichend gewürdigt?

## Ehrliches Verdikt: NEIN, noch nicht — nur auf Flag-Ebene erkannt, nicht im Mapping-Kern gelöst.
- In der Provenienz-Taxonomie ist es benannt (Familie K: „somatische L1/Alu/SVA + HERV").
- Der `chimera_detector` (Block 7) fängt MEI-Breakpoints als **generische** Chimäre — ohne
  TE-Spezifika (keine TSD/poly-A/TE-Familie-Erkennung).
- Die Paralog/PSV/WaveCollapse-Maschinerie ist *konzeptionell* das richtige Werkzeug für die
  Referenz-TE-Mehrdeutigkeit, aber **nicht** speziell für die TE-Skala (~1M Alus) engineered/validiert.
→ Ein Flag ohne echte Detektions-/Mapping-Logik ist hohl. TE verdient mehr.

## Warum es zählt (keine Corner-Case)
- ~50% des humanen Genoms sind TE; ~1.1M Alu-Kopien (300bp, 80-95% intra-Familie-Identität),
  ~500k L1 (6kb), SVA, HERV.
- MEI = große SV-Klasse: polymorphe Alus/L1 in der Population + **somatische L1-Retrotransposition**
  (Tumor, Gehirn) — fakt direkt „Insertionen/Chimären/falsche Translokationen" bei niedriger VAF.

## Zwei VERSCHIEDENE Sub-Probleme (oft vermischt)
**(a) Referenz-TE-Mehrdeutigkeit** (das TE ist in der Referenz, fixiert): ein Read aus einer Alu
„sieht" ~1M nahezu identische Referenz-Alus → Multi-Mapping / MAPQ0. Standard-Mapper picken eine
(falsch) oder droppen (MAPQ0).
**(b) Polymorphe/novel MEI** (TE im Sample, NICHT in der Referenz, oder umgekehrt): Reads über den
Breakpoint = Flanke + TE-Konsensus + TSD + poly-A. Die TE-Hälfte mappt mehrdeutig → fake Chimäre/
Soft-Clip-Cluster/falsche Variante.

## Wie LLmaps Architektur das INTELLEKTUELL RICHTIG lösen würde
**(a) → WaveCollapse-Spread-Mass = ehrliche Positions-Unsicherheit (Lossless-Doktrin):**
Genau hier ist LLmap besser als alle Standard-Mapper. Ein TE-interner Read soll **keine gefakte
eindeutige Platzierung** kriegen, sondern **Wahrscheinlichkeitsmasse über die TE-Kopien behalten**;
ein Read, der die Alu + eine **eindeutige Flanke** überspannt, kollabiert via Flanke auf den echten
Locus. MAPQ spiegelt die **echte** Mehrdeutigkeit (nicht ein gefaktes Unique-Hit oder ein 0-Drop).
Das ist die Stage-1-Self-Interference (Read-vs-Read-Cluster = TE-Reads clustern) + Referenz-
WaveCollapse-Wertversprechen — nur noch nicht TE-skaliert engineered/validiert.

**(b) → dedizierter MEI-Detektor** (Geschwister von junction_hunter/chimera_detector):
Breakpoint-Read = unique-Flanke + TE-Konsensus-Match + **TSD** (Target-Site-Duplication, 5-20bp
flankierende Wiederholung) + **poly-A-Schwanz** (L1/Alu) + 5'-Trunkierung (L1). Diese Signaturen
machen aus „generische Chimäre" → „polymorphe AluYa5-Insertion, 15bp TSD, 30bp poly-A". Provenienz-
Klasse `mei:alu`/`mei:l1`/`mei:sva`/`mei:herv` (strukturell, Layer-1 — eigene Familie, nicht `chim`).

**TE-Familien-Anchor-Katalog:** Alu-Subfamilien (AluY/AluSx/AluYa5…), L1HS, SVA, HERV-K als
Referenz-Anker (RepeatMasker/Dfam) → die TE-Hälfte wird **klassifiziert** (Familie/Subfamilie),
nicht nur „Repeat". Erlaubt subfamily-aware Mapping (junge AluY = polymorph-verdächtig).

## Empfehlung
TE/MEI verdient einen **eigenen Design-Block** (analog zur Mode-Architektur), NICHT nur ein
Provenienz-Flag. Konkret, inkrementell:
1. **Sofort billig:** `mei:*` als eigene Provenienz-Familie (statt MEI in `chim` zu verstecken) +
   TSD/poly-A-Signatur im chimera/MEI-Wrapper. (Kann ich bauen — Geschwister zu meiner Klasse D.)
2. **Mittel:** TE-Familien-Anchor-Katalog (Dfam/RepeatMasker) als Anchor-Quelle.
3. **Kern (das intellektuell Befriedigende):** WaveCollapse explizit für TE-Skala —
   TE-interne Reads behalten ehrliche Spread-Mass, Flanken-Kollaps, MAPQ = echte Unsicherheit.
   Das ist der Punkt, wo LLmaps Lossless-Wave-Particle-Philosophie TE *fundamental* besser macht
   als minimap2/BWA (die faken Placement oder droppen).

**Kurz:** aktuell auf Flag-Ebene erkannt, im Mapping-Kern NICHT gewürdigt. Der richtige Weg ist
nicht „noch ein Flag", sondern die WaveCollapse-Spread-Mass für TE-interne Reads (ehrliche
Unsicherheit statt gefaktes Mapping) + ein TE-bewusster MEI-Detektor mit TSD/poly-A/Subfamilie.

## KRITISCHE Verfeinerung (Agent 2) — MEI ist NICHT ein Topf: Layer-1 vs Layer-3
Wie bei numt-vs-mthet und VDJ darf MEI nicht konflatiert werden:
- **Referenz-TE-Confusion** (Read könnte aus N Referenz-Kopien stammen, falsch platziert) →
  **Layer-1 `mei` ARTEFAKT** (Mapping-Confound, richtig zu buckten).
- **Echte non-Referenz-Insertion** (die Person hat hier WIRKLICH eine neue Alu/L1 — polymorph
  oder somatisch) → **Host-Biologie, die wie ein SV aussieht** → **Layer-3 `bio:mei`** (flaggen,
  NIE als Kontamination buckten) UND ein **Befund** (klinisch relevant: somatische L1-Retro-
  transposition in Tumor/Hirn — der low-VAF-Fall, der den Operator umtreibt), kein Wegwerf-Flag.
Spezialtools dafür: MELT, xTea, TLDR (Longread). LLmap würde es nativ + lossless integrieren.

## MEI-Vertiefungs-Block (priorisierter Follow-up, datengetrieben nach Genom-Test)
(a) TE-Familien-Consensus-Katalog (Dfam/RepeatMasker) + **TE-PSV-Auflösung** (Layer-1 `mei`-
    Confidence statt generisches Repeat — Subfamilien-divergente Positionen = PSV für TE-Familien).
(b) **MEI-Insertions-Detektor:** split-read (Flank + TE-Consensus) + **TSD** + **poly-A** →
    non-Referenz-Insertion → setzt **Layer-3 `bio:mei`** + meldet Insertions-Locus als Befund.
(c) Layer-1/Layer-3-Trennung in der Provenienz-Auflösung (Referenz-Confusion vs echte Insertion).

## Vereinheitlichender Insight: TE-Mapping == Provenienz == derselbe EM-Block
TE-Mapping und der Provenienz-Mode sind **nicht zwei Probleme, sondern dieselbe Kern-Engineering**:
*Wahrscheinlichkeitsmasse über kohärente Alternativen halten, statt einen Unique-Call zu faken.*
- Ein **TE-interner Read** = Masse über ~1M Alu-Kopien, kollabiert via eindeutiger Flanke.
- Ein **Provenienz-Read** = Posterior über {host, paralog, exo, numt, …}.
- Strukturell identisch: beides ist **Spread-Mass + Ko-Segregations-/Flanken-Kollaps im Referenz-
  WaveCollapse**. Die Kontaminations-Buckets-im-EM (gerade gebaut) und die TE-Skala-Auflösung sind
  derselbe Mechanismus.
→ Der TE-Block (iii) ist KEIN Detour, sondern der nächste kohärente EM-Schritt auf der Maschinerie,
  die der Provenienz-Mode schon etabliert hat. Das ist der eigentliche LLmap-Mehrwert ggü.
  minimap2/BWA: ehrliche Unsicherheit (Spread-Mass) statt gefaktes Unique-Mapping — bei TE *und*
  Provenienz mit demselben lossless Wave-Particle-Kern.

## Sequenz (Operator bestätigt „1,2,3 OK")
1. Evidenz-Wiring (Agent 2) → Genom-Test (long+short) — datengetriebene Frequenz von mei/chim/para/numt.
2. TE-Vertiefung als EM-Block: (i) `mei:*`-Familie + TSD/poly-A/Subfamilie-Signatur (Agent 1, billig),
   (ii) TE-Familien-Katalog (Dfam/RepeatMasker), (iii) WaveCollapse TE-skaliert (Spread-Mass).
3. Datengetrieben dort vertiefen, wo der Genom-Test mei/chim als häufig zeigt.
