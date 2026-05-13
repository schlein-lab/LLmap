# Claude as a Tool-Using Agent

LLmap uses Anthropic's Claude (via the Anthropic API or the `claude` CLI subprocess) as a tool-using agent — not as a per-read inference model. Claude has access to `Bash`, `Read/Write`, `WebFetch`, `WebSearch`, and a sandboxed `CUDA Codegen` capability. It is invoked in four deep sessions per analysis.

## Session A — Index-Build Agent (1× per reference, amortized)

**Trigger**: `llmap index --llm sample-aware` or any higher mode.

**Goal**: produce `biology_prior.json` with annotated bucket-level hints.

**Tool-use loop** (typical ~50 internal tool calls):

1. `Bash: bedtools intersect` reference vs RepeatMasker → identify repeat-dense regions
2. `Bash: samtools faidx` extract candidate paralog regions for inspection
3. `WebFetch` UCSC SD-track for current GRCh38/CHM13
4. `WebSearch` recent papers on paralog catalogs for known difficult loci
5. `Write` custom Python preprocessing script for region-specific bucket boundaries
6. `Bash: python custom_bucketize.py` → emit regional overrides
7. `Write biology_prior.json`

**Output schema**:

```json
{
  "version": "1.0",
  "reference_sha256": "...",
  "buckets": {
    "12345": {
      "level": "L2",
      "annotation": "IGH-Constants C-region",
      "prior_weight": 1.2,
      "paralog_partner_bucket": 12346,
      "expected_coverage_multiplier": 2.0,
      "claude_rationale": "Known duplication +650kb in HPRC samples..."
    }
  },
  "regional_overrides": {
    "chr14:105M-107M": {
      "sub_bucket_granularity_kb": 10,
      "max_iter": 25,
      "convergence_threshold": 0.95
    }
  }
}
```

**Cost**: ~$5 amortized over many runs.

## Session B — Sample-Init Agent

**Trigger**: before each `llmap align` run, when `--llm` ≥ `sample-aware`.

**Goal**: choose preset + tune parameters based on dataset properties.

**Tool-use loop** (~20 tool calls):

1. `Bash: seqkit stats` read length distribution
2. `Bash: fastqc` quality profile
3. `Read` sample-sheet metadata if present
4. Reason: "this is iso-seq from PBMC B-cells → expect 5'-truncation, class-switched IGH dominance"
5. `Write sample_params.json` with preset, thresholds, expected coverage profile

**Cost**: ~$1.

## Session C — Diagnostic Agent (with CUDA codegen)

**Trigger**: EM stalls — convergence rate < 10% per iteration for 3 consecutive iterations, OR > 5% of reads non-converged on Level 2.

**Goal**: investigate and resolve.

**Tool-use loop** (~100 tool calls):

1. `Bash: llmap dump-wave-state --batch <id> -o stalled.json`
2. `Read stalled.json` → analyze bucket distribution
3. `Write investigate.py` for custom diagnostic plotting
4. `Bash: python investigate.py` → identify pattern (e.g. bimodal P over two paralogs with insufficient PSV evidence)
5. `WebFetch` literature on the genomic region
6. **(V1.0 capability)** Decide a custom CUDA kernel is needed. `Write custom_kernel.cu` with paralog-aware specialization
7. `Bash: scripts/sandbox_compile.sh custom_kernel.cu -o custom_kernel.so` — invokes the CUDA sandbox
8. Pipeline hot-loads the sandboxed `.so` for the stalled batch only
9. Resumed EM iteration uses the custom kernel

**Cost**: $5-15 per stall, only when triggered.

## Session D — Reporter Agent

**Trigger**: post-run, when `--llm` ≥ `sample-aware`.

**Goal**: generate diagnostic markdown report + update memoization cache.

**Tool-use loop** (~30 tool calls):

1. `Bash: samtools flagstat`, `mosdepth`, etc.
2. `Read` final wave-state summary
3. Reason over patterns: "IGHG4 region shows 4.2× coverage with bimodal cluster — consistent with mosaic duplication"
4. `Write` markdown report
5. `Bash: rocksdb_put` cache generalizable findings

**Cost**: ~$2.

## CUDA Sandbox

Claude-generated CUDA code is sandboxed. V1.0 protections:

1. **Static AST analysis** rejects any `.cu` file containing: system calls, file IO, network access, exec/fork, raw pointer arithmetic outside designated regions
2. **Compile-in-container**: `nvcc` invoked inside rootless `bubblewrap`, no network, read-only filesystem except `/tmp/sandbox-$pid/`
3. **Symbol allow-list**: generated kernels can only call into LLmap's `llmap_kernel_helpers` namespace (a pre-vetted CUDA helper library)
4. **Resource limits**: every generated kernel accepts a `kernel_budget_ns` argument; exceeded → abort + fall back to deterministic path
5. **Audit log**: every generated kernel signed with Claude session ID, sha256, generation timestamp; retained for post-hoc review in `generated_kernels/audit.log`

## Cost summary

| Mode | Sessions | Cost / sample | Speed impact |
|---|---|---|---|
| `off` | — | $0 | baseline |
| `index-only` | A | $0.05 (amortized) | 0 |
| **`sample-aware`** (default) | A + B + D | ~$3 | <1 min (async) |
| `self-healing` | A + B + C-on-demand + D | $5-15 | 0 (async) |
| `research` | A + B + multiple C + D + literature-deep-dive | $20-50 | 0 (async) |

**Critical architectural invariant**: Claude **never** runs synchronously in the GPU hot path. All Claude output is **additive bias** injected into the next EM iteration *if and only if* it is ready in time. The pipeline never blocks waiting for Claude.
