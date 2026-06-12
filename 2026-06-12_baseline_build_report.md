# LLmap Baseline-Build-Report (Agent 2, 2026-06-12)

## Toolchain
- `cmake` fehlte komplett → via `pip3 install --user --break-system-packages cmake` (4.3.2) nach `~/.local/bin`.
- `make` /usr/bin, gcc 13.3.0, kein ninja → Makefile-Generator.
- Configure: `cmake -S . -B build -DLLMAP_ENABLE_{CUDA,FOUNDATION,FAISS,CLAUDE,BENCH}=OFF` (CPU-only). Sauber.

## Build
- `cmake --build build -j$(nproc)` → exit 0. Binary unter `build/src/llmap` (nicht `build/llmap`); `./src/llmap --version` → `llmap 1.0.0`, commit 7f04e76d.

## Tests
- Erstlauf: **95% (1690/1779)**, 89 Fehler — alle CLI-Integrationstests (`VersionCliTest`, `LlmapCliTest`, `ScQcReport`, `Check*`, 2x `RealReferenceTest`).
- **Root cause (Harness-Bug, kein Code-Regress):** `gtest_discover_tests` setzt WORKING_DIRECTORY = `build/tests`. `GetLlmapBinary()` in `tests/unit/test_version.cpp` + `test_llmap_cli.cpp` listet nur `src/llmap`, `../build/src/llmap`, hardcoded `/home/<user>/llmap-local/build/src/llmap` — keiner löst von `build/tests` auf → popen-Fallback `llmap` → exit 127.
- **Fix (nur Test-Harness, kein Production-Code):** Kandidat `../src/llmap` ergänzt (build/tests → build/src/llmap). 2 Dateien.
- Nachlauf: **100% (1779/1779), 0 Fehler.**

## Status
Grüne Baseline steht. Nichts committed (Freigabe ausstehend). Bereit für input_sniffer-Implementierung.
