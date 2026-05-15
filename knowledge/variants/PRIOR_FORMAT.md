# LLmap variant-prior file format (v1)

This is the binary-ish flat-text format that `llmap variant-ingest`
produces and `llmap align` consumes for Layer 4.

## On-disk layout

The format is a TSV with a fixed header line and per-bucket rows. The
file is bgzip-compressed and accompanied by a `.tbi`-style index for
chromosome-keyed random access.

```
# header (required first line)
#format llmap_priors_v1
#source <catalog_id>            # e.g. gnomad_sv_v4.1
#organism <organism>
#assembly <assembly>            # e.g. GRCh38, CHM13v2, mm39
#bucket_bp 100
#min_af 0.0
#n_buckets <int>
# data
chr	pos	n_snv_alt_alleles	snv_max_af	n_sv_overlapping	sv_max_af	sv_summary
chr1	10000	0	0.000	2	0.014	DEL:0.014,INS:0.001
chr1	10100	1	0.002	0	0.000	-
chr1	10200	3	0.310	1	0.003	INS:0.003
...
```

## Columns

| Column | Type | Meaning |
|--------|------|---------|
| `chr` | string | contig name as in reference FASTA |
| `pos` | int    | left edge of the 100-bp bucket (0-based) |
| `n_snv_alt_alleles` | int | count of distinct alt alleles in this bucket across all SNV sources |
| `snv_max_af` | float [0,1] | maximum population AF of any SNV in this bucket |
| `n_sv_overlapping` | int | count of distinct SVs that overlap this bucket (dedup by TYPE+coords) |
| `sv_max_af` | float [0,1] | maximum population AF of any overlapping SV |
| `sv_summary` | string | `TYPE:AF` items, comma-separated; `-` if empty |

`sv_summary` SV types follow the VCF v4.3 SV vocabulary (`DEL`, `INS`,
`DUP`, `INV`, `CNV`, `BND`).

## Bucket resolution

Default 100 bp. Rationale:

- Whole-genome human at 100 bp: ~31 M rows. After bgzip, ~250 MB per
  catalog. Manageable.
- 100 bp matches the typical short-read length and long-read soft-clip
  granularity that the mapper cares about.

Operators may select 10 bp (for high-resolution SV breakpoint priors at
the cost of ~10x file size) or 1 kb (for ultra-light footprint) at
ingest time via `--bucket-bp`.

## Random access

The runtime memory-maps the bgzip file. Per-chromosome offsets are
recorded in a sidecar `.idx` (a sorted `chr<TAB>byte_offset<TAB>n_rows`
table). Query at position `p` on `chr` is:

```
1. chr_offset, chr_rows = idx[chr]
2. bucket_idx = p / bucket_bp
3. record = mmap[chr_offset + bucket_idx * sizeof(row)]
```

The format is line-aligned so that even without the sidecar, binary
search on `pos` works.

## Multi-catalog overlay

`llmap align` accepts a comma-separated list of priors files:

```
--variant-priors dbsnp.priors,gnomad_sv.priors,hprc_panvariants.priors
```

At lookup time, each file is queried independently and the per-position
records are aggregated according to `schema.json:aggregation_rules`. The
overlay is computed lazily, only for positions the mapper actually
touches.

## Empty buckets

To keep file sizes small, buckets with `n_snv_alt_alleles == 0 &&
n_sv_overlapping == 0` are **omitted** from the file. The reader treats
a missing bucket as all-zero. On human dbSNP this typically prunes
30–40 % of buckets.

## Versioning

The `#format llmap_priors_v1` line is a hard requirement. The runtime
will refuse to load files with an unknown version. Adding new columns
will bump to `v2` and provide a migration path.

## Sanity-check tooling

`llmap variant-query --priors X.priors --chr chr14 --pos 105580000
--window 5000` dumps the relevant rows to stdout — useful for spot
checking ingest correctness against known polymorphic loci (e.g. the
IGH locus, the MHC, common CNVs).
