# test_fictive — minimal made-up reference module

- **organism**: `test_fictive` (deliberately fictional — no real biology)
- **reference build**: none — there is no real genome behind this module
- **last_updated**: 2026-05-15
- **maintainer**: LLmap project (knowledge layer / docs)

## Purpose

This module exists as the copy-paste reference example referenced by
`knowledge/EXTENDING.md`. It demonstrates the smallest valid LLmap
organism module: one `MODULE.md`, three `regions/*.json`, one
`classifier_rules.json`. Nothing more is required to make a module
"well-formed" against the schema.

The three regions are invented from scratch and intentionally bear no
relationship to any real organism, repeat family, or gene. They map to
three of the universal taxonomy IDs (`tandem_repeat`, `coding`,
`low_complexity`) so anyone reading the schema can see how to slot a new
region into each common category.

## Made-up biology

We invent a thumb-nail "organism" called `test_fictive` that has:

- a **planet**: Klorpheus IV
- a **chromosome count**: 3 (all linear, named `klorf-1`, `klorf-2`, `klorf-3`)
- a **genome size**: 4.2 megabasepairs
- a **lineage**: synthetic, no biological basis

### Region inventory

| Region | Universal ID | What it is (fictive) |
|--------|--------------|------------------------|
| `zorg_repeat` | `tandem_repeat` | 89 bp pericentromeric tandem array, GC 0.42, Shannon ~1.8 |
| `klorf_coding` | `coding` | protein-coding regions using an "alien" codon usage that prefers CGG / GGC / TAT |
| `null_spacer` | `low_complexity` | homopolymer-ish AT-rich spacers between coding cassettes |

### Files

- `MODULE.md` (this file)
- `regions/zorg_repeat.json`
- `regions/klorf_coding.json`
- `regions/null_spacer.json`
- `classifier_rules.json`

## How to read this module

Treat each file as a copy-paste template. The structure is exactly the
same as a real organism module — only the contents differ:

1. **Region entries** in `regions/*.json` follow `region_entry_schema`
   from `knowledge/schema.json`. Required keys: `name`, `taxonomy_id`,
   `description`, `feature_signature`, `mapping_hints`, `sources`.
2. **Classifier rules** in `classifier_rules.json` are evaluated in
   priority order; the first matching predicate wins. A `default_region`
   is named at the top in case nothing matches.
3. **Sources** — for invented data we use a single placeholder string
   `"synthetic, no biological basis — see EXTENDING.md example."` so the
   field is non-empty and the schema is satisfied.

If you want to author your own organism module, copy this folder, rename
it, and replace the contents. Then follow `EXTENDING.md` step-by-step.
