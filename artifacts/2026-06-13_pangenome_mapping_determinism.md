# Pangenom-Mapping-Reproduzierbarkeit als Konfidenz-Achse (Operator-Idee, 2026-06-13)

Operator: ein Mapping ist **wahrscheinlicher**, wenn dasselbe Mapping über viele Pangenom-
Samples ebenso sinnvoll möglich ist — **unwahrscheinlicher/ehrlich-unsicher**, wenn es überall
ähnlich blurry ist. *Read vs interne Reads desselben Genoms (Stage 1) vs Referenz (Stage 2)
**vs hunderte Genome verschiedener Populationen** = die bessere Metrik.*

## Die Achse: Cross-Population-Mapping-Reproduzierbarkeit
Drei Fälle pro Kandidaten-Platzierung:
- **Reproduzierbar-eindeutig:** der Read kollabiert in ~allen Pangenom-Samples sauber auf den
  homologen Locus → **hohe Konfidenz** (konservierter, eindeutig mappbarer Locus).
- **Reproduzierbar-blurry:** der Read bleibt in ~allen Samples delokalisiert (Spread-Mass) →
  **ehrlich niedrige Konfidenz, ABER intrinsisch** (der Locus ist überall repetitiv) → NICHT
  bestrafen, sondern ehrliche Unsicherheit + Flag „intrinsisch ambig übers Pangenom".
- **Sample-divergent:** hier sauber, übers Pangenom blurry/anders/abwesend → **informativ**:
  echte Variante/SV an dem Locus ODER verdächtige (population-private) Platzierung.

→ Das Pangenom fügt eine **Reproduzierbarkeits-Dimension** zur Collapse-Konfidenz hinzu. Es
trennt „ambig in diesem Genom, aber Locus ist sauber" von „ambig weil der Locus bei ALLEN
Menschen repetitiv ist".

## Das ist der MAPPING-Zwilling von D(pos)
- `D(pos)` (Transcript-Mode): per-Position **Splice**-Determinismus (1−Entropie der Splice-Zustände).
- **`M(pos)` (neu): per-Position MAPPING-Determinismus übers Pangenom** = Cross-Sample-Collapse-
  Reproduzierbarkeit (Anteil der Pangenom-Samples, in denen der Locus eindeutig kollabiert vs
  delokalisiert bleibt). Dieselbe Lossless-Entropie-Philosophie, jetzt über Populationen.

## Implementierung — wie es sich in WaveCollapse einfügt
1. **Precompute `M(pos)` ∈ [0,1]** durch LLmap-Collapse übers Pangenom (HPRC, hunderte
   Haplotypen versch. Populationen) → per genomischer Position die Reproduzierbarkeit des
   eindeutigen Collapse. **Nutzt exakt die Pangenom-Baseline-Compute, die wir gerade aufsetzen**
   (Agent 2's `provenance-baseline`-Pangenom-Lauf).
2. **Als Prior-Faktor in die EM-Update-Formel** (SPEC §2.1): neben π_AI, π_bio, K kommt ein
   **π_pangenome(b)** = M(pos_b) — die Bucket-Likelihood wird gewichtet mit der Pangenom-
   Reproduzierbarkeit dieses Buckets.
3. **Output:** per-Read ein `XM`/`MAPQ`-Beitrag aus `M(pos)` + ein per-Locus-Track
   `pangenome_mappability.bedgraph` (analog zu `splice-determinism`).

## Was es vereint (eine Engine, viele Achsen)
- **Stage 1** (Read vs interne Reads desselben Genoms — Self-Interference): „sehen sich die Reads
  selbst?"
- **Stage 2** (Read vs Referenz): „passt's auf die eine Referenz?"
- **NEU Stage 3** (Read-Collapse vs Pangenom): „reproduziert sich das Mapping über Populationen?"
- + **Provenienz-Baseline** (Pangenom = Erwartungsverteilung) + **D(pos)** (Splice-Determinismus):
  ALLES sind ehrliche Reproduzierbarkeits-/Entropie-Metriken auf derselben WaveCollapse-Maschinerie.
  Spread-Mass (Block 2 iii) ist der lokale Fall; `M(pos)` ist der Cross-Population-Fall.

## Der lossless-kritische Punkt (nicht „bestrafen", sondern EHRLICH)
„Blurry überall → unwahrscheinlicher" heißt **nicht** Penalty/Drop (der minimap2/BWA-Fehler:
gefaktes Unique-Hit oder MAPQ0-Wegwerf). Es heißt: **ehrliche niedrige Konfidenz + Flag
„intrinsisch ambig"** — die Masse bleibt verteilt, MAPQ spiegelt die echte Cross-Population-
Unsicherheit. Sample-divergent wird als Befund (echte Variante vs verdächtig) gemeldet, nicht
verworfen. Das ist die Lossless-Doktrin, auf die Pangenom-Achse angewandt.

## Beziehung zu vorhandener Infra
- CHANGELOG: „Pangenome index format" (geplant) + „Pangenome GAF bridges" (Block 8, gebaut) →
  die Pangenom-als-Referenz-Infra ist teilweise da.
- Die `provenance-baseline`-Pangenom-Compute liefert die Sample-Menge; `M(pos)` ist ein zweiter
  Aggregat-Output desselben Laufs (Mapping-Determinismus statt Provenienz-Spektrum).

## Empfehlung
Eigener Block nach Block 2/3, aber **architektonisch der natürliche nächste EM-Schritt**: derselbe
Spread-Mass-Mechanismus (Block 2 iii) übers Pangenom aggregiert = `M(pos)`-Prior. Reihenfolge:
(1) Pangenom-Lauf läuft eh (Provenienz-Baseline) → `M(pos)` als zweiter Aggregat-Output gratis,
(2) π_pangenome in die EM-Kernel-Update einhängen, (3) MAPQ/Track-Output + Sample-divergent-Befund.
