# Running the pipeline (`screening`)

`screening` is the EE-finding pipeline. It takes a genome FASTA plus the protein
database, its metadata CSV, and the host-gene baits, and produces the endogenous
element sequences, taxonomy table, GFF3 annotation, and flanking regions.

(example-run)=
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

A run goes through the twelve stages below — the same ones `eefinder.log`
records, in this order, with a timing for each. Every stage writes its result to
disk before the next picks it up, so any of them can be inspected afterwards
under `tmp_files/`.

1. **Prepare input data** — prefix every FASTA header (`>PREFIX/…`) and drop
   contigs shorter than `--length`.
2. **Index databases** — build the BLAST or DIAMOND indexes for `-db` and `-bt`
   (only with `--index_databases`).
3. **Similarity search** — search the genome against `-db`, then collapse
   redundant hits by query, coordinate range and strand (`--range_junction`).
   `--translation_method` controls how the genome is translated — see
   [Translation methods](translation-methods.md).
4. **Extraction of putative EEs** — cut the surviving candidate regions out of
   the genome.
5. **Filter step** — search those candidates against the host baits `-bt` and
   drop every candidate whose best bait hit outscores its best `-db` hit.
6. **Get basic taxonomy** — join the surviving hits to the metadata CSV. Its
   header is validated before the run starts (see
   [the required format](get-databases.md#the-metadata-csv-format)).
7. **Merge truncated elements** — join neighbouring fragments of the same taxon
   and strand (`--limit`, `--merge_level`) and re-extract the merged sequences.
8. **Clean EEs** — optional soft-mask filter (`--clean_masked`).
9. **Create final taxonomy** — build one row per element, flag overlapping
   elements and add `Average_pident`.
10. **Filter overlapping elements** — resolve overlaps by the chosen strategy —
    see [Overlap resolution](overlap.md).
11. **Generate GFF3 annotation** — write the taxonomy table as GFF3.
12. **Extract flanking regions** — `--flank` nt on each side of every element.

```{note}
The similarity search therefore runs **twice** — once for the genome against
`-db` (stage 3) and once for the candidates against `-bt` (stage 5) — always
with the same `--mode` and `--translation_method`.
```

(diamond-sensitivity)=
## The sensitivity trade-off

```{important}
DIAMOND is faster than BLAST, but **at a real cost in sensitivity**, and the
loss falls exactly where endogenous-element studies operate: on the most
divergent sequences.
```

DIAMOND seeds its alignments on a reduced amino-acid alphabet. Its authors
report that this does not compromise sensitivity, but the EEfinder benchmark
found a clear effect on the highly divergent sequences typical of EVE studies.
Screening the *Aedes aegypti* Aag2 genome (GCA_021653915) against a viral
protein database gave:

| Search mode | Elements recovered | Identity range | Runtime |
|-------------|--------------------|----------------|---------|
| `blastx` (BLAST) | **481** | 10.2–100 % | 56 h 14 min |
| `-md very-sensitive` (DIAMOND) | 225 | 16.2–100 % | 42 h 34 min |
| `-md fast` (DIAMOND) | 126 | 16.9–100 % | 8 h 6 min |

DIAMOND's most sensitive setting recovered **less than half** the elements BLAST
did, for a runtime saving of only ~25 %; `fast` mode recovered about a quarter of
them, but seven times faster. The identity ranges show where the difference
comes from — BLAST reaches down to 10.2 % identity, while both DIAMOND modes
floor out around 16 %, i.e. the oldest and most degraded integrations are the
ones that go missing.

**Use `blastx` (the default) whenever sensitivity matters**, which for EE
discovery is almost always. The DIAMOND modes are appropriate for rapid
exploratory screens, for very large sample sets, or when only comparatively
recent, well-conserved integrations are of interest — with the caveat that the
resulting element counts are not comparable to BLAST-based ones.

The full benchmark is reported in section 3.2 (*Benchmark of alignment tools*) of
[Dias, Dezordi & Wallau (2024)](https://doi.org/10.1016/j.csbj.2024.10.012)
([PMC11532726](https://pmc.ncbi.nlm.nih.gov/articles/PMC11532726/)).

## Options reference

### Required inputs

| Option | Meaning |
|--------|---------|
| `-in/--genome_file` | Input genome FASTA (nucleotides). |
| `-od/--outdir` | Output directory. |
| `-db/--database` | Protein database FASTA (virus or bacteria). |
| `-mt/--dbmetadata` | Protein metadata CSV for `-db`. Its seven columns are validated before the run starts — see [the required format](get-databases.md#the-metadata-csv-format). |
| `-bt/--hostgenesbaits` | Host-gene bait proteins FASTA. |

### Search & filtering

| Option | Default | Meaning |
|--------|---------|---------|
| `-md/--mode` | `blastx` | `blastx` or a DIAMOND sensitivity (`fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`). DIAMOND is faster but less sensitive — see [The sensitivity trade-off](#diamond-sensitivity). |
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
