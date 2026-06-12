# input_sniffer — Design-Sketch (Agent 1, 2026-06-12)

Zweck: `TranscriptMode::Auto` → konkreter Modus auflösen. Quelle: `docs/design/llmap_mode_architecture.md` §3. Dependency-light (stdlib + llmap_core).

## API (`src/io/input_sniffer.h`, namespace `llmap::io`)

```cpp
enum class FileFormat : uint8_t { Unknown, Fasta, Fastq, Sam, Bam };

struct FastaStats {            // nur bei FASTA befüllt
    uint64_t n_seqs   = 0;     // gezählt im Sample
    uint64_t median_len = 0;
    uint64_t n50      = 0;
    bool     sampled_truncated = false;  // >sample_limit Seqs vorhanden
};

struct SniffResult {
    FileFormat            format = FileFormat::Unknown;
    core::TranscriptMode  mode   = core::TranscriptMode::GenomeReads;
    std::string           reason;   // für [mode-detect] log line
    std::optional<FastaStats> fasta_stats;
};

// Format aus magic bytes / erster Zeile. Liest nur die ersten ~1 KiB.
FileFormat SniffFormat(const std::string& path);

// FASTA-Stats über die ersten <=sample_limit Sequenzen (default 1000).
FastaStats ComputeFastaStats(const std::string& path, uint64_t sample_limit = 1000);

// Voll-Auflösung. mode_override != Auto  ->  direkt zurück (Override gewinnt, geloggt).
// reads_path + assembly_path beide gesetzt -> ReadsVsAssembly.
SniffResult ResolveMode(const std::string& primary_path,
                        core::TranscriptMode mode_override,
                        bool has_reads, bool has_assembly);
```

## Heuristik (Doc §3, deterministisch, KEIN Pilot-Pass v1)

Reihenfolge in `ResolveMode`:
1. `mode_override != Auto` → return override (reason="override").
2. `has_reads && has_assembly` → ReadsVsAssembly.
3. `format = SniffFormat(primary)`:
   - **BAM/SAM**: `@PG`/Header-Token in {isoseq3,lima,splice,STAR,cDNA,FLNC} → Transcript; sonst GenomeReads.
     (Pilot-Pass „≥10% long-N CIGAR" = **out of scope v1**, als TODO markiert — braucht Alignment-Pass.)
   - **FASTQ**: Basename matcht `(?i)flnc|isoseq|cdna|rna` → Transcript; sonst GenomeReads.
   - **FASTA**: stats; `median>50k && n50>100k && n_seqs<5000` → Assembly;
     sonst `median in [300,15000] && n_seqs>50000` → Transcript (FLNC-as-FASTA);
     sonst `median<300` → GenomeReads (short-read FASTA); sonst GenomeReads (low-conf, geloggt).
   - **Unknown** → GenomeReads (reason="format unknown, default").

## Magic bytes (SniffFormat)
- BAM: gzip `0x1f 0x8b` UND (sicherste, dep-freie Heuristik) Pfad endet `.bam` ODER inflated beginnt `BAM\1`. v1: gzip-magic + `.bam`-Endung → Bam; gzip ohne `.bam` → behandelt als evtl. gz-FASTA/FASTQ (dekomprimiert erstes Byte? **v1: nicht** — gz-Input ist separates Thema, reason="gzip non-bam: unsupported sniff, default GenomeReads"). FRAGE an Reviewer: reicht das, oder gz-FASTQ häufig genug für v1?
- SAM: erste Zeile beginnt `@HD`\t / `@SQ`\t / `@PG`.
- FASTQ: erstes Zeichen `@` UND (3. Zeile beginnt `+`).
- FASTA: erstes Zeichen `>`.

## N50 / median (ComputeFastaStats)
- Sammle Längen der ersten `sample_limit` Records (memory: nur Längen, nicht Seq).
- median = sortierte Mitte; n50 = kleinste Länge L mit kumsum(>=L sortiert desc) >= total/2.
- `sampled_truncated=true` wenn Datei mehr Records hat → Heuristik bleibt valide (Stichprobe).

## Tests (`tests/io/test_input_sniffer.cpp`)
Pro Branch ein synthetischer Mini-Input (in-memory tmpfile via `std::tmpfile`/`testing::TempDir`):
- FASTA assembly (3 Contigs à 200 kb) → Assembly
- FASTA FLNC-signature (60k Seqs à ~1.5 kb) → Transcript  *(skaliert: nur Längenliste, keine echten 60k Seqs nötig — n_seqs via Header-Count)*
- FASTA short-read (median 150) → GenomeReads
- FASTQ basename `sample.flnc.fastq` → Transcript
- FASTQ basename `reads.fastq` → GenomeReads
- SAM `@PG ... minimap2 -ax splice` → Transcript
- SAM `@PG ... bwa` → GenomeReads
- Override transcript bei FASTA-assembly-Input → Transcript (override gewinnt)
- has_reads+has_assembly → ReadsVsAssembly
- median/n50-Berechnung Unit (bekannte Längenliste → erwarteter N50)

## Offene Fragen an Reviewer (Agent 2)
1. gz-FASTQ-Sniff v1 nötig, oder reicht plain + `.bam`-Endung?
2. Pilot-Pass (long-N CIGAR) wirklich out-of-scope v1? (braucht Mini-Alignment, zieht Scope auf)
3. FLNC-as-FASTA-Schwelle `n_seqs>50000` mit Stichprobe `sample_limit=1000` nicht direkt messbar → wir bräuchten echten Record-Count (billig: nur `>`-Zeilen zählen, ganze Datei einmal scannen). OK so?
