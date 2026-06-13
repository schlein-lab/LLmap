# Multi-Level-Konsens + Assemble-then-Map (Operator, 2026-06-13)

Zwei Verfeinerungen der Kaskade. Beide stark; (2) ist der größte Speed+Recall-Hebel bisher.

## (1) Hinlänglich = wenn 2–3 aufeinanderfolgende Ebenen dasselbe Ergebnis liefern
Statt eines Einzel-Stufen-Kriteriums: der Read ist gelöst, wenn **N aufeinanderfolgende
Kaskaden-Ebenen auf dieselbe Platzierung konvergieren**. Konsens = Konfidenz = Stop.
- **Robuster als mein Einzel-Stufen-`IsProvablyResolved`** — es ist eine Kreuz-Validierung:
  zwei unabhängige Ebenen, die übereinstimmen, sind verlässlicher als eine. **Adressiert genau das
  Falsch-zu-früh-Exit-Risiko**, das Agent 2 benannt hat.
- Die Ebenen existieren schon als **Bucket-Pyramide (coarse L2 → fine L0)**: stimmen coarse-EM und
  fine-EM überein → Stop. Oder seed-chain ↔ EM. Oder (mit (2)) contig-map ↔ read-map.
- Implementierung: das Resolution-Kriterium der Kaskade wird „letzte K Ebenen-Ergebnisse identisch"
  statt eines Single-Stage-Schwellwerts. Mein `sufficiency`-Modul trackt dann Agreement über Stufen.

## (2) Assemble-then-Map: Stage-1-Reads-gegeneinander → Mini-Contigs → die mappen (nicht die Reads)
Operator: Runde 1 mappt Reads gegeneinander (Self-Interference) → das erzeugt Mini-Contig-
Assemblies → die Read-Anzahl drastisch reduzieren, Zeit sparen, und durch **längere Einheiten
eindeutiger mappen**.

**LLmap HAT die Infra schon halb (Phase 2–3):**
- **Stage-1 Self-Interference** (`allpair`): Read-vs-Read-Clustering (FAISS-kNN + Leiden) → Cluster
  ähnlicher Reads. ✅ gebaut.
- **Cluster-Representative** (Medoid) + **Member-Propagation** (Phase 3.5: „cluster-rep mapping →
  member alignment via cheap intra-cluster banded-WFA delta"). ✅ gebaut.
- **Operator-Verfeinerung:** statt nur den *Medoid-Read* zu mappen, **assembliere den Cluster zu
  einem Consensus-Mini-Contig** (länger + fehler-gemittelt als jeder Einzel-Read) → mappe den
  Contig → Reads erben die Platzierung via Member-Propagation.

**Warum das BEIDE Probleme löst:**
- **Speed (der Extension-Engpass löst sich auf):** statt ~10M Reads × 17 Extensions → ~100k Contigs
  × wenige Extensions. Die 99,8%-Extension-Zeit kollabiert, weil man **Größenordnungen weniger
  Einheiten** extendet. Das ist der eigentliche 1000×-Hebel — nicht GPU, nicht nur Early-Exit,
  sondern **weniger, längere Einheiten**.
- **Recall (das `molecule/0`-91bp-Problem verschwindet):** ein Contig, assembliert aus ALLEN Reads
  eines Transkripts, ist eine **lange Einheit, die die volle Multi-Exon-Struktur trägt** → mappt als
  EIN spliced Unit eindeutig, spannt die Introns. Das 91bp-Fragment-Problem (Locus-Selektion pickt
  ein kurzes Stück) löst sich, weil der Contig genug Länge hat, um den echten Locus eindeutig zu
  ankern. Längere Einheiten = höhere `M(pos)`-Mappbarkeit per Konstruktion.

**Die zwei Ideen greifen ineinander:** die Kaskade läuft auf **Contigs statt Reads** (2), und das
Resolution-Kriterium ist **Multi-Level-Konsens** (1). Weniger Einheiten × Early-Exit-bei-Konsens =
multiplikativer Speedup, bei *besserer* Eindeutigkeit (längere Einheiten + Kreuz-Validierung).

## Ehrliche Caveats
- **Cluster-Korrektheit:** mis-geclusterte Reads → falscher Contig → falsche Platzierung für alle
  Member. Die Stage-1-Konsens-Qualität wird kritisch (aber: WaveCollapse-Self-Interference ist genau
  dafür gebaut, und Member-Propagation verifiziert per banded-WFA-Delta gegen den Contig).
- **Coverage-abhängig:** high-coverage-Loci → große Contigs → riesige Ersparnis; Singletons/low-cov
  → kein Contig, individuell mappen (Fallback). Also: Hybrid — Contigs wo möglich, Reads wo nötig.
- **Assembly-Kosten:** der Cluster→Consensus-Schritt kostet, aber **einmal pro Cluster**, amortisiert
  über alle Member. Netto massiv positiv bei realer Coverage.

## Empfehlung
Das ist die **architektonisch sauberste Lösung für Speed UND Recall zugleich** — und sie reitet auf
Stage-1/Member-Propagation, die schon stehen. Reihenfolge:
1. **Cluster→Consensus-Mini-Contig** (neuer Schritt nach `allpair`-Clustering; Consensus aus den
   Cluster-Reads).
2. **Map Contigs statt Reads**, Reads erben via bestehender Member-Propagation.
3. **Resolution = Multi-Level-Konsens** (mein `sufficiency` → Agreement über Stufen/Ebenen).
4. Das ersetzt sowohl den GPU-Umweg als auch das Single-Stage-Early-Exit durch eine fundamentalere
   Reduktion (weniger, längere, eindeutigere Einheiten + Kreuz-Validierung).
Komplementär: der Recall-Fix (Agent 2) ist dann „mappe den Contig richtig" statt „rette 17 Fragment-
Chains" — das `molecule/0`-Problem wird durch die Contig-Länge gegenstandslos.
