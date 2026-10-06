# Outputs

A `screening` run writes its main results directly into `--outdir`, plus a JSON
run log and (unless `--removetmp`) an archive of intermediates.

A single stage (`eefinder screening align`, …) writes only its own outputs,
under the stage names listed below; the `PREFIX.EEs.*` names and the archive are
produced by `all`.

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
| `element_id` | `CONTIG:START-END` identifier of the element (BED-style, 0-based half-open). The run prefix is **not** part of it — see the note below. |
| `sense` | Strand of the element, `pos` or `neg`. |
| `protein_ids` | `;`-separated `accession|pident` entries for the best-matching reference protein(s). |
| `protein_products` | Product name(s) from the metadata CSV. |
| `molecule_type` | Genome composition of the source (e.g. `ssRNA(+)`). |
| `taxonomy` | Lineage of the assignment — see [the lineage of a merged element](#the-lineage-of-a-merged-element). |
| `host` | Host recorded for the reference. |
| `overlaped_element_id` | Element(s) this one overlaps, if any. |
| `tag` | `overlaped` or `unique`. |
| `average_pident` | Mean percent identity of the hits backing the element. |

```{note}
The three places an element is named use two spellings. The taxonomy table's
`element_id` drops the run prefix (`ctg_1913:1754-2689`), while the `EEs.fa` /
`EEs.flanks.fa` headers and the GFF3 `ID` attribute carry it
(`Ae_aeg_Aag2_ctg_1913/ctg_1913:1754-2689`), so the sequences cross-reference the
genome the run was given. Prepend `{prefix}/` to an `element_id` to match a
FASTA record.
```

### The lineage of a merged element

An element is one interval, but `--merge_level` can build it from hits of
several taxa. Its `taxonomy` is then the **lowest common ancestor** of their
lineages: ranks are kept while every hit names the same taxon, and from the
first rank they disagree on down, the element carries no assignment. Two
species of one genus therefore resolve to the genus, two genera to the family.

The unresolved tail is left off rather than padded, so an element resolved only
to the family reads

```
r__Riboviria;k__Pararnavirae;p__Artverviricota;c__Revtraviricetes;o__Ortervirales;f__Caulimoviridae
```

A rank that *none* of the hits names is not a disagreement — viral lineages
often have no realm — so it stays `Unk` in the middle of the string and the
ranks below it still count. An element whose hits share nothing is
`Unclassified`. The rank vocabulary itself is described under
[the `Taxonomy` column](get-databases.md#the-taxonomy-column).

## The GFF3 annotation (`PREFIX.EEs.gff3`)

Each row of the taxonomy table becomes one feature, sorted by sequence id then
start. The `element_id` supplies the sequence id and the coordinates, converted
from the BED-style 0-based half-open range EEfinder carries internally to GFF3's
1-based inclusive convention. `average_pident` becomes the feature score, and
`--analysis` selects the feature type (`endogenous_viral_element` or
`endogenous_bacterial_element`).

Column 9 carries the taxonomy as attributes — reserved tags capitalised per the
spec, custom tags lower-case:

| Attribute | Source column |
|-----------|---------------|
| `ID` | `element_id`, prefixed with the run prefix |
| `Name` | the `s__` rank of `taxonomy` |
| `realm`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` | one attribute per rank of `taxonomy`; ranks the element does not reach are omitted |
| `molecule_type` | `molecule_type` |
| `product` | `protein_products` |
| `protein_ids` | `protein_ids`, re-separated with `,` (see below) |
| `host` | `host` |
| `overlap_status` | `tag` (`overlaped` / `unique`) |

The ranks are written as separate attributes rather than as the single lineage
string the TSV carries, so a genome browser can filter on one of them.

GFF3-reserved characters (`;`, `=`, `&`, `,`, tabs, newlines) are
percent-encoded in attribute values. The `;` that separates the `protein_ids`
entries in the TSV is one of them, so in the GFF3 those entries are joined with
`,` instead — the separator the spec reserves for a multi-valued attribute. A
comma occurring *inside* an entry is still escaped to `%2C`, so the two cannot
be confused.

## The cleaned outputs (`--clean_masked`)

`PREFIX.EEs.cleaned.fa` holds the subset of elements whose lower-case (soft-masked)
plus `N` content is at or below `--mask_per` percent, and
`PREFIX.EEs.cleaned.tax.tsv` is the taxonomy table restricted to those elements.
Both are a **subset** of the unfiltered pair, which is always written — the flag
adds a view rather than replacing one.

## The run log (`eefinder.log`)

`eefinder.log` is a JSON document recording:

- `eefinder_version` — the installed EEfinder version;
- `arguments` — the resolved run arguments (including `translation_method`);
- `dependencies` — detected versions of bedtools, BLAST, DIAMOND, the NCBI
  `datasets` CLI, python, numpy and pandas, each flagged if it differs from the
  `env.yml` pin;
- per-step and total timing information.

```{tip}
When EEfinder is installed outside its source tree, point the dependency-drift
check at the reference file with `export EEFINDER_ENV_YML=/path/to/env.yml`.
```

## Intermediate files (`tmp_files/`)

Unless `--removetmp` is given, the intermediates are archived under `tmp_files/`.
They are named after the stage that produced them, so the file a stage reads
and the file it writes are both obvious. `ScreeningPaths` in
[the Python API](api.md#file-names) is the single source of truth for these
names.

```text
outdir/
├── eefinder.log
├── PREFIX.EEs.fa
├── PREFIX.EEs.tax.tsv
├── PREFIX.EEs.gff3
├── PREFIX.EEs.flanks.fa
└── tmp_files/
    ├── PREFIX.cleaned_genome.fa          # prefixed headers + length-filtered
    ├── PREFIX.cleaned_genome.lengths.tsv
    ├── PREFIX.ee_hits.tsv                # similarity search (main)
    ├── PREFIX.ee_hits.filtered.tsv       # redundant-hit filter
    ├── PREFIX.ee_candidates.bed
    ├── PREFIX.ee_candidates.fa           # putative EEs
    ├── PREFIX.host_hits.tsv              # host-bait search
    ├── PREFIX.host_hits.filtered.tsv
    ├── PREFIX.ee_hits.validated.tsv      # survived the bait comparison
    ├── PREFIX.taxonomy_signature.csv
    ├── PREFIX.merge.annotated.bed
    ├── PREFIX.merge.merged.bed
    ├── PREFIX.merge.elements.bed
    ├── PREFIX.flanks.bed
    └── PREFIX.flanks.slop.bed
```

With the prediction-based translation methods (`gv`/`rv`/`gv-rv`), the
predicted-protein FASTAs and coordinate TSVs also appear here, named after the
query each search used — `PREFIX.cleaned_genome.fa.pred.coords.tsv` for the main
search and `PREFIX.ee_candidates.fa.pred.coords.tsv` for the host-bait one.
