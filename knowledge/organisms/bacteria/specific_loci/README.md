# Bacteria specific_loci — Layer 2 region database

Per-species, instance-level locus annotations for ten reference bacterial genomes
used by LLmap. Each subdirectory holds JSON cards for one species; each card
describes one named, coordinate-bound locus on the indicated GenBank
reference accession.

This is Layer 2 of the LLmap three-layer region-knowledge model — see
`../../../human/specific_loci/README.md` for the layer description.

## Bacterial vs. mammalian conventions

Bacteria use circular chromosomes (and accessory replicons) referenced by
GenBank/RefSeq accession rather than assembly+chromosome. The `coordinates`
block therefore uses an `accession_<NC_...>` key instead of `grch38`/`chm13`:

```json
"coordinates": {
  "GenBank_NC_000913.3": {
    "chr": "NC_000913.3",
    "start": 4035531,
    "end":   4040906,
    "strand": "+"
  }
}
```

Many bacterial mobile elements (SCCmec types, integron classes, Tn3-family
transposons, SXT ICE) are *element families* whose individual occurrences
are isolate-specific. Those are stored as **reference cards** with the
sentinel `coordinates.<id>.chr == "<element>_consensus"` and `start == 0,
end == 1`, plus `extras.reference_card_only == true`. They contribute
mapping-hint context without claiming a coordinate.

## Mapping-hint conventions

- `rrn` operons — `lambda_scale: 1.5`, `max_occ: 20000`,
  `require_psv_disambig: true`. Near-identical paralogs; ITS spacer
  provides PSV signal.
- IS elements — `lambda_scale: 1.4`, `max_occ: 30000`,
  `is_element: true`, `require_psv_disambig: true`.
- Pathogenicity islands / prophages / ICEs — `lambda_scale: 1.3-1.4`,
  `max_occ: 15000`. Lineage-restricted but with paralog signal across
  isolates.
- PE/PPE/PGRS (M. tb) — `lambda_scale: 1.6`, `max_occ: 50000`,
  `anchor_weight_scale: 0.3`. The PGRS tandem extensions are essentially
  impossible to map uniquely with short reads.
- CRISPR arrays — `lambda_scale: 1.7`, `max_occ: 10000`,
  `tandem_array_handling` implicit via the repeat-spacer structure.
- V. cholerae chromosomal superintegron — `lambda_scale: 1.8`,
  `max_occ: 50000`, `tandem_array_handling: true`.
- DNA uptake sequences (Neisseria DUS) — `lambda_scale: 1.4`,
  `max_occ: 80000`, treated as ubiquitous short repeat.

## Directory layout

```
specific_loci/
├── README.md                    # this file
├── ecoli_k12_mg1655/            # NC_000913.3 — 55 cards
├── saureus_n315/                # NC_002745.2 — 33 cards (HA-MRSA)
├── saureus_usa300/              # NC_007793.1 — 18 cards (CA-MRSA)
├── mtb_h37rv/                   # NC_000962.3 — 47 cards (PE/PPE-rich)
├── spneumo_tigr4/               # NC_003028.3 — 17 cards
├── kpneumo_hs11286/             # NC_016845.1 — 42 cards (KPC, integrons)
├── paeruginosa_pao1/            # NC_002516.2 — 21 cards
├── bsubtilis_168/               # NC_000964.3 — 20 cards (10 rrn, ICEBs1)
├── vcholerae_o1/                # NC_002505/06.1 — 19 cards (CTXphi, SI)
└── nmeningitidis_mc58/          # NC_003112.2 — 28 cards (opa, pilE/pilS)
```

Total: 300 locus cards.

## Sources at a glance

- NCBI RefSeq genome annotations (each species' RefSeq record)
- ISfinder database — https://www-is.biotoul.fr
- ICEberg 2.0 — doi:10.1093/nar/gky1123
- CARD (Comprehensive Antibiotic Resistance Database) — doi:10.1093/nar/gkz935
- PHASTER (PHAge Search Tool Enhanced Release) — doi:10.1093/nar/gkw387
- Blattner et al. 1997 (E. coli MG1655) — doi:10.1126/science.277.5331.1453
- Kuroda et al. 2001 (S. aureus N315) — doi:10.1016/S0140-6736(01)05321-1
- Diep et al. 2006 (USA300) — doi:10.1016/S0140-6736(06)68231-7
- Cole et al. 1998 (M. tb H37Rv) — doi:10.1038/31159
- Tettelin et al. 2001 (S. pneumoniae TIGR4) — doi:10.1126/science.1061217
- Liu et al. 2012 (K. pneumoniae HS11286) — doi:10.1128/JB.06487-11
- Stover et al. 2000 (P. aeruginosa PAO1) — doi:10.1038/35023079
- Kunst et al. 1997 (B. subtilis 168) — doi:10.1038/36786
- Heidelberg et al. 2000 (V. cholerae N16961) — doi:10.1038/35020000
- Tettelin et al. 2000 (N. meningitidis MC58) — doi:10.1126/science.287.5459.1809
- IWG-SCC SCCmec typing — doi:10.1128/AAC.00579-09

## Coordinate caveats

Many of the per-copy IS-element and rrn-operon start/end values are drawn from
the relevant RefSeq feature tables but rounded to integer coordinates for
readability. They are intended for *locus-class lookup* (does a read hit a
known multi-copy element?), not for nucleotide-precise mapping. Refine
against the live RefSeq GFF if a downstream tool needs exact bounds.
