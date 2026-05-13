# Contributing to LLmap

LLmap is in active V1.0 autonomous build — most commits between May and June 2026 come from a Claude-Code driver running on a 96-hour loop. Human contributions are welcome and review-prioritized.

## Quick start

```bash
git clone git@github.com:schlein-lab/LLmap.git
cd LLmap
mkdir -p build && cd build
cmake -DLLMAP_ENABLE_CUDA=OFF -DLLMAP_USE_NATIVE_ARCH=OFF ..
cmake --build . -j$(nproc)
ctest --output-on-failure
```

This builds in CPU-only mode for fast development. For GPU work you need CUDA 12.3+ and Compute Capability ≥ 8.0.

## Code style

- C++23, modern idioms (`std::optional`, `std::span`, `std::expected` where appropriate)
- No raw `new`/`delete`. RAII everywhere.
- No comments explaining WHAT — only WHY for non-obvious decisions.
- gtest for unit tests. One assertion per behavioral aspect.
- Module layout: `src/<module>/<name>.{h,cpp}` + `tests/unit/test_<name>.cpp`.

## The lossless invariant is non-negotiable

`tests/unit/test_alignment_record.cpp` encodes the V1.0 promise. Do not weaken it. If you need to extend `AlignmentRecord`, update the invariant check and its tests.

## How autonomous commits work

The driver script (`scripts/autonomous_driver.sh`) runs every 15 min via cron. Each iteration:

1. Pulls latest from `main`
2. Reads `STATE.md` for `current_task`
3. Spawns `claude --print` as a subprocess with `scripts/continuation_prompt.md`
4. Claude does one concrete substep, updates `STATE.md`, commits
5. Driver pushes
6. Zyrkel notifications on Hummel-up/down transitions and on stalls

If you see a stuck-driver Zyrkel alert, check `autonomous_run.log` and `STATE.md` `## Blockers` section.

## PRs from humans

- Open PRs against `main` directly
- Tests must pass (CI runs build + ctest)
- Keep PRs small — autonomous commits land frequently; rebase is your friend

## Reporting issues

[GitHub Issues](https://github.com/schlein-lab/LLmap/issues). Tag with phase (`phase-0` through `phase-9`) so the autonomous driver can prioritize.
