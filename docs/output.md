# Outputs

A `screening` run writes its main results directly into `--outdir`, plus a JSON
run log and (unless `--removetmp`) an archive of intermediates.

## Main outputs

| File | Contents |
|------|----------|
| `PREFIX.EEs.fa` | Endogenous element nucleotide sequences. |
| `PREFIX.EEs.tax.tsv` | Endogenous element taxonomy table (one row per element). |
| `PREFIX.EEs.gff3` | Endogenous element annotation in GFF3 format. |
| `PREFIX.EEs.flanks.fa` | Endogenous elements plus flanking regions (`--flank` nt each side). |
| `eefinder.log` | JSON run summary (see below). |

With `--clean_masked`, the mask-cleaned equivalents are also written:
`PREFIX.EEs.cleaned.fa` and `PREFIX.EEs.cleaned.tax.tsv`.

With `--overlap longest|targets`, the elements filtered out of the final results
are preserved under `tmp_outputs/` as `PREFIX.EEs.removed.{fa,tax.tsv}`.

## The taxonomy table (`PREFIX.EEs.tax.tsv`)

| Column | Meaning |
|--------|---------|
| `Element-ID` | `PREFIX\|CONTIG:START-END` identifier of the element. |
| `Sense` | Strand of the element (`+`/`-`). |
| `Protein-IDs` | Accession(s) of the best-matching reference protein(s). |
| `Protein-Products` | Product name(s) from the metadata CSV. |
| `Molecule_type` | Genome composition of the source (e.g. `ssRNA(+)`). |
| `Family` / `Genus` / `Species` | Taxonomic assignment from the metadata CSV. |
| `Host` | Host recorded for the reference. |
| `Overlaped_Element_ID` | Element(s) this one overlaps, if any. |
| `tag` | `overlaped` or `unique`. |
| `Average_pident` | Mean percent identity of the hits backing the element. |

## The run log (`eefinder.log`)

`eefinder.log` is a JSON document recording:

- `eefinder_version` — the installed EEfinder version;
- `arguments` — the resolved run arguments (including `translation_method`);
- `dependencies` — detected versions of bedtools, BLAST, DIAMOND, python, numpy
  and pandas, each flagged if it differs from the `env.yml` pin;
- per-step and total timing information.

```{tip}
When EEfinder is installed outside its source tree, point the dependency-drift
check at the reference file with `export EEFINDER_ENV_YML=/path/to/env.yml`.
```

## Intermediate files (`tmp_files/`)

Unless `--removetmp` is given, the intermediates are archived under `tmp_files/`.
Their names accrete suffixes as they pass through the pipeline, so you can trace
exactly which step produced each file:

```text
outdir/
├── eefinder.log
├── PREFIX.EEs.fa
├── PREFIX.EEs.tax.tsv
├── PREFIX.EEs.gff3
├── PREFIX.EEs.flanks.fa
└── tmp_files/
    ├── PREFIX.rn                     # prefixed headers
    ├── PREFIX.rn.fmt                 # length-filtered
    ├── PREFIX.rn.fmt.blastx          # similarity search (main)
    ├── PREFIX.rn.fmt.blastx.filtred  # redundant-hit filter
    ├── PREFIX.rn.fmt.blastx.filtred.bed
    ├── PREFIX.rn.fmt.blastx.filtred.bed.fasta          # putative EEs
    ├── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx   # host-bait search
    ├── ...
    └── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge...
```

With the prediction-based translation methods (`gv`/`rv`/`gv-rv`), the
predicted-protein coordinates TSVs also appear here (e.g.
`PREFIX.rn.fmt.pred.coords.tsv` for the main search and
`PREFIX.rn.fmt.blastx.filtred.bed.fasta.pred.coords.tsv` for the host-bait
search) — evidence that the translation method was applied to both searches.
