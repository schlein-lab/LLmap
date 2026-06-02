#!/usr/bin/env python3
"""Validate every T1 entry against catalog/schema/curated.schema.json.

Usage:
    python3 catalog/validate.py              # validate all
    python3 catalog/validate.py <file.json>  # validate one
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

try:
    import jsonschema
except ImportError:
    sys.exit("ERROR: pip install jsonschema (or 'conda install -c conda-forge jsonschema')")


CATALOG = Path(__file__).resolve().parent
SCHEMA_PATH = CATALOG / "schema" / "curated.schema.json"
CURATED_DIR = CATALOG / "curated"


def load_schema():
    with SCHEMA_PATH.open() as f:
        return json.load(f)


def validate_file(path: Path, schema: dict) -> tuple[bool, str]:
    try:
        with path.open() as f:
            entry = json.load(f)
    except json.JSONDecodeError as e:
        return False, f"JSON parse error: {e}"
    try:
        jsonschema.validate(entry, schema)
    except jsonschema.ValidationError as e:
        return False, f"Schema validation failed at {'/'.join(map(str, e.absolute_path))}: {e.message}"

    # Cross-field invariants beyond JSON-Schema reach
    expected_id = path.stem
    if entry.get("locus_id") != expected_id:
        return False, f"locus_id '{entry.get('locus_id')}' does not match filename '{expected_id}'"
    return True, "OK"


def main(argv: list[str]) -> int:
    schema = load_schema()

    if len(argv) > 1:
        files = [Path(p).resolve() for p in argv[1:]]
    else:
        files = sorted(CURATED_DIR.glob("*.json"))

    if not files:
        print(f"WARN: no T1 files in {CURATED_DIR}")
        return 0

    failures = 0
    for f in files:
        ok, msg = validate_file(f, schema)
        sigil = "✓" if ok else "✗"
        print(f"{sigil} {f.relative_to(CATALOG.parent)}: {msg}")
        if not ok:
            failures += 1

    print(f"\n{len(files) - failures}/{len(files)} valid; {failures} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
