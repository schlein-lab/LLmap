# Adaptive Feature-Cascade — welche Features wann, bucket-/wahrscheinlichkeits-getrieben (Operator, 2026-06-13)

Operator: ein Algorithmus, der entscheidet **wann man welche Features anwendet**, basierend auf
Buckets + Wahrscheinlichkeiten, mit dem Ziel, den Read **bei der ersten (billigsten) Anwendung zu
lösen**. → Eine per-Read **adaptive Kaskade** mit Early-Exit.

## Das Prinzip
Nicht alle teure Maschinerie auf jeden Read. Stattdessen: Features in **Kosten-Reihenfolge**, jeder
Read durchläuft sie, **steigt aus, sobald seine Wahrscheinlichkeitsmasse kollabiert** (= gelöst),
und das **nächste Feature wird vom aktuellen Bucket-/Wahrscheinlichkeits-Zustand gewählt**.

**Das ist NICHT die binäre Triage** (die der Operator zu Recht abgelehnt hat — sie würde leichten
Reads die WaveCollapse-Ehrlichkeit verweigern). Hier gibt es **keine Vorab-Klassifikation**: jeder
Read tritt in dieselbe Pipeline ein, eskaliert aber **dynamisch nur so weit wie nötig**.

## Die zwei Entscheidungen
1. **Resolution-Kriterium (wann aussteigen?)** = das **WaveCollapse-Collapse-Kriterium selbst**
   (SPEC §2.2: Entropie/Dropout). Masse kollabiert auf einen dominanten Bucket (niedrige Entropie,
   großer Score-Gap) → gelöst → Exit. Spread → weiter.
2. **Feature-Selektion (welches Feature als nächstes?)** = vom **Bucket-Muster** diktiert:
   - Masse über **Paralog-Buckets** verteilt → **PSV** anwenden.
   - Masse über **Repeat/TE-Buckets** → **M(pos)/Spread-Mass** (Pangenom-Prior).
   - **kein guter Bucket** → **AI/Foundation-Embedding**.
   - **exogen-aussehend** (kein Host-Bucket, aber Exo-Panel-Hit) → **Provenienz-Detektoren**.
   - **intron-artige Lücke** → **Spliced-Chain-Joiner** (Transcript).
   Der Bucket-Zustand **wählt das auflösende Feature** statt blind alle zu fahren.

## Feature-Kosten-Reihenfolge (billig → teuer)
seed-chain (µs) → coarse-EM (Bucket-Pyramide L2) → fine-EM (L0) → PSV/Paralog → Spliced-Joiner →
AI/Foundation-Embedding → Provenienz-Detektoren. Jede Stufe nur, wenn die vorige nicht gelöst hat.

## Das existiert schon HALB — der Checkpoint-Dispatcher (Phase 7)
LLmap hat bereits den `checkpoint_dispatcher`, der **entscheidet WANN zur LLM eskaliert wird**, auf
Basis des Bucket-/Chain-Zustands (AmbiguousChain, unknown region, paralog disambiguation, SD
expansion, novel insertion). **Das IST eine Feature-Selektions-Entscheidung — nur für EIN Feature
(LLM).** Die Operator-Idee = **diesen Mechanismus auf ALLE Features verallgemeinern** (EM-Level,
PSV, AI, Provenienz, M(pos)), nicht nur die LLM. Die Infra (Dispatcher, der auf Bucket-State
reagiert) ist da; sie muss von „1 Feature" auf „Feature-Kaskade" erweitert werden.

## Warum das Speed löst (ohne Triage, ohne Ehrlichkeits-Verlust)
- ~95% der Reads lösen auf **Feature 1 (seed-chain)** → minimap2-Speed.
- Nur die Residual-Ambigen eskalieren — und genau die KRIEGEN die teuren Features (kein Bypass).
- Net: µs/Read im Schnitt, **kein Read wird entwertet** (jeder kriegt genau die Features, die seine
  Ambiguität auflösen). Das ist der „ohne Triage"-Speedup, den der Operator wollte — adaptiv, nicht
  binär. Komplementär zum GPU-Batch-EM (das macht die EM-Stufe selbst schnell; die Kaskade ruft sie
  seltener auf).

## Policy: erst Hand, dann gelernt
- **V1 (Hand):** Kosten-geordnete Kaskade + Collapse-Kriterium als Exit + Bucket-Muster→Feature-
  Regeln (oben). Deterministisch, sofort baubar.
- **V2 (gelernt):** welches Feature reduziert die Entropie am meisten für einen gegebenen Bucket-
  State — eine kleine gelernte (oder LLM-)Policy. Das ist die „LL"-Augmentierung im Namen.

## Beziehung zu allem Gebauten
Vereint: WaveCollapse (Resolution-Kriterium) · Bucket-Pyramide (coarse→fine = billig→teuer) ·
alle Feature-Module (chain/PSV/AI/Provenienz/M(pos)/Joiner) als Kaskaden-Stufen ·
checkpoint_dispatcher (die Eskalations-Entscheidung, verallgemeinert). Es ist die **Meta-Schicht
über der ganzen Architektur** — und löst Speed + bewahrt Lossless gleichzeitig.

## Empfehlung
Eigener Block, aber **architektonisch der Schlussstein**: der `checkpoint_dispatcher` wird zum
generischen `feature_cascade_dispatcher` (Bucket-State → nächstes Feature | resolved). Reihenfolge:
(1) Resolution-Kriterium aus dem bestehenden Collapse-Check exponieren, (2) Feature-Registry mit
Kosten + Bucket-Muster-Trigger, (3) Dispatcher-Loop (apply → check-resolved → select-next), (4)
später gelernte Policy. Komplementär zu GPU-Batch-EM (Speed der Stufe) + Recall-Fix (Qualität der
chain-Stufe).
