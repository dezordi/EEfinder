# Running the pipeline (`screening`)

`screening` is the EE-finding pipeline. It takes a genome FASTA plus the protein
database, its metadata CSV, and the host-gene baits, and produces the endogenous
element sequences, taxonomy table, GFF3 annotation, and flanking regions.

## Example run (with the bundled `test_files/`)

The repository ships a small example dataset in
[`test_files/`](https://github.com/WallauBioinfo/EEfinder/tree/master/test_files)
so you can try the full pipeline end-to-end:

| File | Role | CLI flag |
|------|------|----------|
| `Ae_aeg_Aag2_ctg_1913.fasta` | Query genome contig (*Aedes aegypti* Aag2) | `-in` |
| `virus_subset.fa` | Viral protein database | `-db` |
| `virus_subset.csv` | Metadata for the viral database | `-mt` |
| `filter_subset.fa` | Host-gene bait proteins | `-bt` |

From the repository root, with the `EEfinder` environment active:

```bash
eefinder screening \
  -in test_files/Ae_aeg_Aag2_ctg_1913.fasta \
  -od results_test \
  -db test_files/virus_subset.fa \
  -mt test_files/virus_subset.csv \
  -bt test_files/filter_subset.fa \
  -ln 1000 \
  -id \
  -p 2 \
  -lm 100
```

- `-ln 1000` lowers the minimum contig length (the example contig is shorter than
  the 10000 nt default, which would otherwise filter it out).
- `-id` builds the BLAST/DIAMOND indexes for the databases (needed on the first
  run against a given database).
- `-p 2` uses two threads.
- `-lm 100` merges same-taxon elements within 100 nt.

## Pipeline steps

The `screening` command orchestrates these steps (each a small side-effect class;
files flow through disk with accreting suffixes `.rn`, `.fmt`, `.blastx`,
`.filtred`, `.bed`, `.tax`, …):

1. **InsertPrefix** — prefix every FASTA header (`>PREFIX/…`).
2. **RemoveShortSequences** — drop contigs below `--length`.
3. **MakeDB** — build BLAST or DIAMOND databases (`--index_databases`).
4. **SimilaritySearch** — the similarity search, run **twice** (main EE search +
   host-bait search). `--translation_method` controls both — see
   [Translation methods](translation-methods.md).
5. **FilterTable** — filter redundant hits by `qseqid`/range/sense.
6. **GetFasta** — extract putative EE sequences (bedtools).
7. **CompareResults** — drop EEs that hit host baits harder.
8. **GetTaxonomy** — join hits to the metadata CSV, build the taxonomy table.
9. **MergeBed** — merge truncated elements of the same genus/family
   (`--merge_level`).
10. **MaskClean** — optional soft-mask filter (`--clean_masked`).
11. **TagElements** — flag overlapping elements, add `Average_pident`.
12. **FilterOverlap** — resolve overlaps by the chosen strategy — see
    [Overlap resolution](overlap.md).
13. **WriteGFF3** — write the EE taxonomy table as a GFF3 annotation.
14. **GetLength + BedFlank + GetFasta** — extract flanking regions (`--flank`).

## Options reference

### Required inputs

| Option | Meaning |
|--------|---------|
| `-in/--genome_file` | Input genome FASTA (nucleotides). |
| `-od/--outdir` | Output directory. |
| `-db/--database` | Protein database FASTA (virus or bacteria). |
| `-mt/--dbmetadata` | Protein metadata CSV for `-db`. |
| `-bt/--hostgenesbaits` | Host-gene bait proteins FASTA. |

### Search & filtering

| Option | Default | Meaning |
|--------|---------|---------|
| `-md/--mode` | `blastx` | `blastx` or a DIAMOND sensitivity (`fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`). |
| `-tm/--translation_method` | `default` | `default`/`gv`/`rv`/`gv-rv` — see [Translation methods](translation-methods.md). |
| `-ln/--length` | `10000` | Minimum contig length for the search. |
| `-rj/--range_junction` | `100` | Range for junction of redundant hits — see [Custom arguments](custom-arguments.md). |
| `-p/--threads` | `1` | Threads for multi-threaded steps. |
| `-id/--index_databases` | off | Build the BLAST/DIAMOND indexes for the databases. |

### Element assembly

| Option | Default | Meaning |
|--------|---------|---------|
| `-lm/--limit` | `1` | Bases used to merge neighbouring same-taxon elements (bedtools merge). |
| `-ml/--merge_level` | `family` | Taxonomic level (`family`/`genus`) to merge by. |
| `-fl/--flank` | `10000` | Flanking-region length to extract. |
| `-ov/--overlap` | `keep` | `keep`/`longest`/`targets` — see [Overlap resolution](overlap.md). |
| `-tf/--target_families` | — | Family to KEEP (repeatable) with `--overlap targets`. |
| `-ntf/--non_target_families` | — | Family to DROP (repeatable) with `--overlap targets`. |

### Masking & output

| Option | Default | Meaning |
|--------|---------|---------|
| `-mp/--mask_per` | `50` | Lowercase-percentage threshold to call a region repetitive. |
| `-cm/--clean_masked` | off | Also emit mask-cleaned outputs (`*.cleaned.*`). |
| `-an/--analysis` | `virus` | GFF3 feature type (`virus` → `endogenous_viral_element`, `bacteria` → `endogenous_bacterial_element`). |
| `-pr/--prefix` | input filename | Prefix for output files and Element-IDs. |
| `-rm/--removetmp` | off | Delete intermediates instead of archiving them under `tmp_files/`. |
| `--debug` | off | Emit verbose debug logging (intermediate paths, per-step details). |

```{note}
The default `blastx` mode is the tested, reliable path. The DIAMOND modes
depend on the `diamond` build pinned in `env.yml` and can fail silently if that
build misbehaves — verify it if a DIAMOND run produces no hits.
```

See [Outputs](output.md) for the files produced and
[Custom arguments](custom-arguments.md) for the merge/junction behaviour with
worked examples.
