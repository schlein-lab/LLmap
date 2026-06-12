# P2 — Transcript spliced-dispatch (Agent 1, 2026-06-13)

Ziel: `llmap align --mode transcript` läuft die Spliced-Pipeline end-to-end.
Baut auf P1 (`resolved_mode` ist gesetzt). Quelle: llmap_mode_architecture.md §2.1,
chain_spliced.h, spliced_likelihood.h, splice_site_db.h.

## Kernidee
Der klassische Chainer zerlegt einen spliced Read (hard `max_gap_diff`) in **mehrere
`ClassicalAlignment`s** desselben (ref_name, strand). Diese SIND die Sub-Chains.
Die Transcript-Stage joint sie über Intron-Lücken zu EINER spliced Alignment.

## Neues Modul: `src/mapping/transcript_stage.{h,cpp}` (in llmap_mapping)
Reiner Post-Stage, keine Mutation des Chainers (lossless-konform, DNA-Pfad unberührt).

```cpp
namespace llmap::mapping {
struct TranscriptStageConfig {
    JoinerConfig joiner;                 // aus chain_spliced.h
    float min_junction_prob = 0.30f;
};
// refs: ref_name -> sequence (für Donor/Acceptor-Motiv an der Lücke).
// Mutiert die alignments eines Reads in-place: ersetzt N Sub-Alignments durch
// 1 spliced Alignment mit N-CIGAR + Junction-Metadaten. Gibt Junction-Liste
// zurück (für jI/jM-Tags).
struct SplicedReadResult {
    bool spliced = false;
    std::string spliced_cigar;
    char strand = '+';
    std::vector<std::pair<uint64_t,uint64_t>> junctions; // donor,acceptor
    std::vector<float> junction_conf;
};
SplicedReadResult ApplyTranscriptStage(
    std::vector<classical::ClassicalAlignment>& alns,   // ein Read
    const std::function<std::string_view(const std::string&)>& ref_lookup,
    const annot::SpliceSiteDb& splice_db,
    const TranscriptStageConfig& cfg);
}
```

## Junction-Prob (P2: PWM-only, KEIN Annotations-Zwang)
Pro Lücke zwischen Sub-Chain i und i+1 (gleiche ref, gleicher strand, query-adjazent):
1. donor_2bp = ref[upstream.ref_end .. +2)   (Intron-Start, '+'-Strang)
2. acceptor_2bp = ref[downstream.ref_start-2 .. ) (Intron-Ende)
   (Reverse-Strang: revcomp + getauschte Rollen.)
3. `SpliceScoreResult s = splice_db.ScoreJunction(donor_2bp, acceptor_2bp, ...)`
4. junction_prob = sigmoid-Mix aus s.donor_score+s.acceptor_score, Floor 0.05
   (NIE 0 — lossless: Sub-Chain wird nie still gedroppt; siehe JoinSplicedChains
    `n_singletons_kept`).
→ Die annotations-aware `JunctionProbability(ExonBoundary, JunctionEvidence)`
  (braucht GENCODE) ist P3 mit `--annotation`.

## Join + CIGAR
`JoinSplicedChains(subchains, junction_probs, cfg.joiner)` → `SplicedChainResult`.
Pro `SplicedChain` mit >1 Sub-Chain: `EmitSplicedCigar()` → ein CIGAR mit N-Ops.
Sub-Chains, die nicht joinen, bleiben als Singletons (eigene Alignments) erhalten.

## Wiring in cmd_align_run.cpp
- `RunAlignmentBatches` bekommt `resolved_mode` + eine Ref-Lookup (ref_name→seq).
  Ref-Seqs werden in cmd_align.cpp ohnehin via `ref_reader` geladen → als
  `std::unordered_map<std::string,std::string>` (oder Span auf MmapFasta) durchreichen.
- Nach `pipeline.AlignReads(...)`, vor `ConvertClassicalAlignment`:
  `if (resolved_mode == Transcript) ApplyTranscriptStage(res.alignments, lookup, splice_db, cfg);`
- `ConvertClassicalAlignment` muss den spliced CIGAR übernehmen (kommt schon via
  `aln.cigar`, da wir die alignments in-place ersetzen → keine Signatur-Änderung).

## Output-Tags (P2 minimal: XS/jI/jM)
In bam_writer-Pfad: wenn spliced, Tags via `output::transcript_schema::{XsTag,JiTag,JmTag}`.
Voller Tag-Satz (TI/XK/XC/XA/XM/XF) + Parquet-Sidecars = P3.

## CMake
`llmap_mapping` muss an `llmap` gelinkt werden (aktuell NICHT, P1 brauchte es nicht):
Zeile 349 src/CMakeLists.txt → `llmap_mapping llmap_annot` (annot für SpliceSiteDb,
ist aber schon gelinkt). Plus `target_link_libraries(llmap_mapping PUBLIC llmap_classical)`
falls transcript_stage `ClassicalAlignment` braucht (Dep prüfen — evtl. nur Typen kopieren
wie LinearSubChain es tut, um Dep-Zyklus mapping↔classical zu vermeiden).

## Tests (lokal, synthetisch, NAS-frei)
`tests/mapping/test_transcript_stage.cpp`:
- 2 Sub-Chains, GT..AG-Motiv an der Lücke, Intron 1kb → 1 spliced Alignment, CIGAR `...N...`.
- Nicht-kanonisches Motiv → join nur wenn prob>=0.30, sonst 2 Singletons.
- Reverse-Strang-Fall (revcomp-Motive).
- DNA-Read (1 Alignment) → unverändert (Regression).
- jI/jM-Tag-Strings korrekt.

## Offene Designpunkte für Reviewer (Agent 2)
1. Ref-Durchreichung: `unordered_map<string,string>` (einfach, RAM) vs `MmapFasta`-Handle
   (zero-copy, aber Lebensdauer/Thread-Safety). Vorschlag: map für P2, mmap später.
2. Dep-Richtung mapping→classical: zieht transcript_stage `ClassicalAlignment` direkt rein
   (Link mapping→classical) ODER konvertieren wir in cmd_align_run zu `LinearSubChain`
   und transcript_stage bleibt classical-frei? **Vorschlag: Letzteres** — cmd_align_run
   baut `vector<LinearSubChain>` aus den ClassicalAlignments, transcript_stage kennt nur
   mapping-Typen. Hält den sauberen Dep-Graph (mapping hängt nicht an classical).
3. Strand-Handling bei `-`: Motiv-Extraktion + N-CIGAR-Reihenfolge. Korrektheit kritisch.
```
