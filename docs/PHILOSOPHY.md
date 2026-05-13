# LLmap Philosophy — the double meaning of LL

## LosslessLmap

Every input read produces an `AlignmentRecord`. There is no "read got dropped because MAPQ=0" failure mode. The record's `status` is one of:

- `Mapped` — collapsed to a single bucket with confidence ≥ τ
- `Tentative` — best candidate below threshold, but with seeds and a probability distribution preserved
- `Unmapped` — no seeds anywhere; record carries a `RejectionReason`

This is enforced **at construction time** in C++ via factory functions `make_mapped`, `make_tentative`, `make_unmapped`. It is checked at runtime via `is_lossless_consistent()`. It is verified at integration-test time via the bulk-invariant test in `tests/unit/test_alignment_record.cpp`.

The lossless property is the foundation on which everything else stands. If `count(input_reads) != count(AlignmentRecord)`, LLmap has a bug — not a feature.

## LLM-Map

Claude (Anthropic) is a first-class architectural component, not an afterthought. Claude is invoked as a **tool-using agent** with `Bash`, `Read/Write`, `WebFetch`, and `CUDA Codegen` capabilities — not as a per-read voter (which would be prohibitively expensive).

Four Claude sessions per analysis (see [CLAUDE_AGENT.md](CLAUDE_AGENT.md)):

- **Session A** at index-build: Claude reasons about the reference, fetches UCSC tracks, writes custom preprocessing scripts, generates `biology_prior.json`
- **Session B** at sample-init: Claude reads FASTQ metadata, picks preset + parameters
- **Session C** on EM stall: Claude investigates, can write a custom CUDA kernel, sandboxed compile + hot-load
- **Session D** post-run: Claude generates the diagnostic report

LLmap is the first production bioinformatics tool that lets Claude generate CUDA kernels in the algorithmic hot path. The kernels are sandboxed (bubblewrap, no network, symbol allow-list, resource limits, signed audit log) but real.

## The photon analogy is real math, not branding

Reads are not particles to be located. They are wavefunctions over a hierarchical bucket space, evolving under a Hamiltonian composed of sequence likelihood, coverage coupling, AI prior, and biology prior. The wavefunction collapses to a determined position only when probability mass concentrates above threshold.

See [PHYSICS.md](PHYSICS.md) for the formal mapping: path integrals, symmetry breaking, RG flow, decoherence-T2.

## Design ethos

- **Honesty about limits**: LLmap's `What it does NOT do` section in the README is non-trivial. We will explicitly return `uninformative` for reads that cannot be disambiguated.
- **Drop-in compatibility**: BAM-Compat output by default; `samtools`/`bcftools`/IGV unchanged. Researchers can adopt LLmap without rewriting pipelines.
- **Speed is non-negotiable**: AI features must never slow the GPU pipeline. Claude runs async, additive only.
- **Failure modes that page**: any non-trivial dropout, stall, or anomaly produces a Zyrkel notification. Silent failure is the only unacceptable failure.
