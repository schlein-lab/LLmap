You are continuing the LLmap autonomous build, iteration N. LLmap is a C++/CUDA sequence mapper with the WaveCollapse algorithm — see LLmap_SPEC.md for the full specification and STATE.md for current progress.

## Your job each iteration

1. **Read STATE.md** in the current working directory. Find the `current_task` block — that is what to do this iteration.
2. **Do one concrete substep** of that task. Pick a small unit of work that completes in ≤ 10 minutes wallclock. Examples:
   - implement one header + impl file + matching gtest
   - write one CUDA kernel + a unit test
   - extend an existing module with a missing feature
   - fix a failing test
3. **Build and test locally**:
   ```bash
   cd build || (mkdir -p build && cd build && cmake -DLLMAP_ENABLE_CUDA=OFF -DLLMAP_USE_NATIVE_ARCH=OFF ..)
   cmake --build . -j$(nproc)
   ctest --output-on-failure
   ```
   If the build or tests fail, FIX them before committing. Do not commit broken code.
4. **Update STATE.md**:
   - Move the completed substep to "Phase progress"
   - Set the next substep as `current_task`
   - Append to "Driver invocations log" with iteration number, timestamp, action, outcome
   - Set `last_successful_iteration` to N
5. **Commit** with a message like `phase X.Y: <what>` and let the driver push.

## Heavy compute → Hummel-2

If a task needs GPU, real reference data, or > 8 GB RAM:
- SSH to `hummel-login`
- Workspace: `~/llmap` (the cloned repo) and `/beegfs/u/bbg6775/llmap/` for data/builds
- Submit via SLURM `sbatch`
- DO NOT block the driver waiting — submit job, log job ID, mark task as `awaiting_slurm_job_<id>` in STATE.md, return. Next iteration will poll the job.

If Hummel is down (check `.hummel_status` file), do only local-host work this iteration: docs, host-side C++ code, tests that don't need GPU.

## Constraints — do not break

- Build must always pass `cmake --build .` cleanly on the local host (CUDA OFF mode).
- `ctest --output-on-failure` must pass for every commit (zero failing tests).
- The lossless invariant test in `tests/unit/test_alignment_record.cpp` MUST stay green — it encodes the V1.0 promise.
- Never commit foundation-model weights, BAM files, or anything > 10 MB to git. Use `.gitignore` for those.
- Never push to `main` if tests fail.
- If you hit an algorithmic dead-end you can't resolve, write the blocker in STATE.md → `## Blockers / open questions`, send a Zyrkel notification via the `mcp__zyrkel__notify` tool, and exit cleanly. The next iteration will see the blocker.

## HARD RULE — before any new code each iteration

Run this command FIRST every iteration:

```bash
find src -name '*.cpp' -exec wc -l {} \; | awk '$1 > 400' | sort -rn
```

If the list is non-empty (which it currently is — 10 files), this iteration MUST split at least one of them per the modular rule below, in addition to (or instead of) advancing the current phase task. **A Phase 3.x commit alone is not sufficient if any src/*.cpp > 400 LOC exists.** Bundle the refactor into the same commit or do it as a preceding commit.

Refactor target naming: split file `foo.cpp` into `foo_<aspect>.cpp` (e.g. `faiss_wrapper.cpp` → `faiss_wrapper_index.cpp` + `faiss_wrapper_search.cpp` + `faiss_wrapper_serialize.cpp`). All sub-files compiled together via the same CMake target. Tests still mirror the original logical unit.

Acceptance for every iteration: at least one of {(a) the monolith list shrank, (b) the list was empty already}. If neither, the iteration is rejected — undo your changes and split one file instead.

## Modular architecture — NO MONOLITHS

This is an explicit user-set constraint. Build small, focused modules — never god-classes or 1000-line files.

- **One concept per header**: if a header has more than ~3 unrelated classes/structs, split it.
- **Soft cap 400 LOC per .cpp file**: above that, refactor into helpers in a sibling file. Hard limit 700 LOC.
- **Free functions over methods** when the function doesn't own state.
- **Clear module boundaries**: anything in `src/core/` cannot depend on `src/ai/` or `src/claude_agent/`. The dependency graph flows core → classical → ai → reference_collapse → output → cli.
- **Headers ≤ 200 LOC**: if a header grows past 200, the public surface is too wide — split it.
- **PIMPL when a header would otherwise need a heavy third-party include** (e.g. ONNX Runtime, FAISS). The pattern is already used in `foundation_embedder.cpp` — follow it.
- **One CMake target per logical unit**: don't pile everything into one static lib. Sub-libraries help compile-time and enforce boundaries.
- **Avoid utility/helper grab-bag files**: `utils.h` is a code smell. If you find yourself adding a third unrelated function to one file, split it.
- **Test files mirror source files 1:1**: `src/foo/bar.cpp` ↔ `tests/unit/test_bar.cpp`. No "test_misc.cpp".

When in doubt: prefer two 200-line files over one 400-line file.

## Style

- C++23 modern style; `std::optional`, no raw new/delete, RAII.
- No comments explaining WHAT — only WHY for non-obvious decisions.
- gtest for unit tests. One assertion per behavioral aspect.
- Match existing file layout: `src/<module>/<name>.{h,cpp}` + `tests/unit/test_<name>.cpp`.

## Output to driver

When done, just exit normally. The driver reads `STATE.md` and the git log to verify progress. If you cannot make progress, write to STATE.md `## Blockers` and exit.

Now: read STATE.md and start.
