# LLmap 1000× Speedup — Vorschläge (Operator-Frage, 2026-06-13)

Ausgangslage: ~1 s/Read (Benchmark). minimap2: µs/Read. 1000× nötig zum Konkurrieren.
**Kern-Einsicht: 1000× kommt NICHT aus Mikro-Opts, sondern daraus, die teure Maschinerie
NICHT auf die ~95% einfachen Reads zu werfen.** Triage/Tiered-Architektur ist der Hebel.

## SCHRITT 0 (zuerst, kein Raten): Profilen
Die Phase-8-`Profiler`-Infra existiert. Wo gehen die 1 s/Read hin? Kandidaten:
- WaveCollapse-EM-Iterationen? · WFA2-Extension? · der 459-Chain-Soup (Real-Data-Befund)? ·
  per-Read-LLM/Checkpoint-Dispatcher? · ONNX-Foundation-Embedding? · Stage-1-Self-Interference?
**Erst messen, dann den dominanten Kostenpunkt angreifen.** Wahrscheinlichste Schuldige: LLM/ONNX
per-Read (Sekunden!) + EM auf einfachen Reads + die Chain-Soup.

## Die Hebel (multiplikativ), nach Wirkung sortiert

**(A) Read-Triage / Fast-Path — DER größte Hebel (~20–50×)**
~95% der Reads sind eindeutig (Top-Chain dominiert, Score-Gap groß, unique Locus). Die brauchen
KEINE WaveCollapse, KEIN LLM — nur classical seed-chain-extend = **minimap2-Speed**. Nur die
genuinely-ambigen Reads (Paralog/Repeat/Kontamination-Kandidaten) gehen in die teure WaveCollapse/
AI-Schicht. Das ist exakt LLmaps Philosophie: der Mehrwert ist auf den HARTEN Reads, nicht den
leichten. → ein `triage`-Gate nach dem Chaining: dominanter Hit → Fast-Path-Output; sonst → EM.

**(B) Kein per-Read-LLM/Foundation im Batch (~10–100×, falls drin)**
Wenn die 1 s/Read einen LLM-Consult oder ONNX-Embedding pro Read enthält, ist DAS der Killer.
Batch-Default: `--classical-only` + WaveCollapse-on-demand; LLM-Checkpoint NUR auf genuinely-
stuck Reads (selten), nie per-Read. Foundation-Embedding nur im AI-Pfad (s.u. GPU-Batch).

**(C) Frühes Kandidaten-Pruning — kill den 459-Chain-Soup VOR der Extension (~5–20× auf Repeats)**
Real-Data: ~459 Chains/Read auf repeat-reichen Regionen, meist spurious. Locus-Selektion (Block 2)
prunt schon, aber wenn sie NACH der Extension läuft, bleibt die Extension-Kosten. → Pruning auf
Minimizer-/Chain-Score-Ebene VOR der teuren WFA2-Extension (nur Top-k Loci extenden).

**(D) GPU-Batched AI-Pfad (~100–1000× für den AI-Teil)**
Wenn der AI-Pfad GEBRAUCHT wird (harte Reads), batchen: Foundation-Embedder + FAISS-ANN bei
batch=10k → SPEC-Ziel **10 µs/Read** (statt per-Read-CPU-Fallback = ms–s). Die Architektur KANN
das (Phase 1/2); 1 s/Read ist vermutlich der nicht-gebatchte/CPU-Debug-Pfad.

**(E) M(pos) als Triage-Signal (gratis aus der Operator-Idee)**
Der precomputed `M(pos)`-Mappbarkeits-Track sagt VORAB, welche Loci ambig sind. Read landet in
high-M(pos) (unique) → Fast-Path; low-M(pos) (Repeat) → WaveCollapse. So spart man die per-Read-
Ambiguitäts-Berechnung selbst — die Triage-Entscheidung ist ein O(1)-Lookup. Schöne Vereinigung:
M(pos) ist Konfidenz- UND Speed-Signal.

**(F) Systems (großteils schon da, Phase 8): konstante Faktoren**
SIMD-Minimizer, Arena-Allocator, mmap-Index, Work-Stealing-Threads — drin. Plus: 47 Samples ×
Cores (Batch-Parallelität). Konstant-Faktor, nicht 1000×, aber multipliziert auf den Rest.

## Net-Rechnung
(A) Fast-Path ~30× × (B) kein-per-Read-LLM ~30× = **~900× allein für den Common Case**, + (C)
Pruning auf Repeats + (D) GPU-Batch für den AI-Rest → **>1000× erreichbar**. Der Schlüssel ist
(A)+(B): die 95% einfachen Reads laufen mit minimap2-Speed, WaveCollapse/AI NUR auf den harten.

## OPERATOR-KORREKTUR: 1000× geht OHNE Triage (der bessere Weg)
Triage (A) ist ein Hebel, aber NICHT nötig — und sie hat einen Preis: einfache Reads bekämen
die WaveCollapse-Ehrlichkeit nicht (nur klassische Platzierung). **Der No-Triage-Pfad behält die
volle lossless WaveCollapse auf ALLEN Reads und holt den Speed aus dem Batching:**
- **PRIMÄR: GPU-gebatchte fused WaveCollapse-EM.** SPEC §2/§3 zielt schon darauf: fused CUDA-
  Kernels (`em_iteration`/`collapse_check`/`refinement`), **10 µs/Read bei batch=10k**. Die ~1 s/Read
  sind der CPU-Fallback/nicht-gebatchte Pfad — der GPU-Batch-Pfad IST die Architektur. Jeder Read
  läuft durch die echte WaveCollapse, amortisiert über den Batch → µs/Read, ohne Bypass/Entwertung.
- Das ersetzt (A) Triage als primären Hebel; (A) bleibt optional. (B) kein-per-Read-LLM/ONNX +
  (C) frühes Pruning + (D) GPU-Batch-Embedding/ANN + (F) Systems bleiben gültig + multiplizieren.
- Architektur (WaveCollapse-Core) ist korrekt (Operator bestätigt) — der Fix ist die EM-
  IMPLEMENTIERUNG (GPU/SIMD-Batch), nicht das Architektur-Design.

## Reihenfolge
1. **Profilen** (Phase-8-Profiler) → dominanten Kostenpunkt bestätigen (Vermutung: LLM/ONNX
   per-Read + EM auf easy Reads).
2. **Triage-Gate** (A) + Batch-Default ohne per-Read-LLM (B) — der ~900×-Block.
3. **Pre-Extension-Pruning** (C) + **M(pos)-Triage-Lookup** (E).
4. **GPU-Batch-AI-Pfad** (D) für die harten Reads.
Dann ist der pangenom-weite minimap2-Benchmark compute-feasibel.
