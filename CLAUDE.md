# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What EEfinder is

EEfinder is a Python CLI that automates identification of **Endogenous Elements
(EEs)** — virus- or bacteria-derived sequences integrated into eukaryotic
genomes — via similarity search plus genomic-junction reasoning. Published in
*Computational and Structural Biotechnology Journal* (Dias, Dezordi & Wallau,
2024; https://doi.org/10.1016/j.csbj.2024.10.012).

Wiki: https://github.com/WallauBioinfo/EEfinder/wiki

## Architecture

The package lives in `eefinder/` (flat layout, no `src/`). Each processing step
is a small class whose `__init__` runs the work as a side effect (files in,
files out) — there is no shared in-memory pipeline object; steps communicate
through files on disk whose names accrete suffixes (`.rn`, `.fmt`, `.blastx`,
`.filtred`, `.bed`, `.tax`, ...).

The CLI (`eefinder/scripts/main.py`) is a `click` **group** (`cli`, the console
entry point) with two commands:

- **`screening`** — the EE-finding pipeline (everything below); the `screening`
  function was formerly the single top-level command.
- **`get-databases`** — a `click` **group** that downloads the RefSeq protein
  databases (`get_databases.py` `GetDatabases`) via the NCBI `datasets` CLI. It
  has one subcommand per database, each with group-specific defaults/options:
  `virus` (taxon default `10239`, `--exclude-uninformative` +
  `--standardize-proteins`) and `bacteria` (taxon default `2`,
  `--exclude-uninformative`, `--standardize-proteins` = generic cleaning only,
  no name map) produce a protein FASTA + metadata CSV (the
  `-db`/`-mt` inputs); `host` (taxon **required**) produces the `-bt` baits
  FASTA. The shared `-od`/`-pr`/`--refseq` options come from the
  `_common_download_options` decorator in `main.py`; `_run_get_databases` does
  the `datasets`-binary check and calls `GetDatabases`. The metadata CSV is
  rebuilt from the `protein.faa` headers (Accession/Protein/Species) joined with
  the `data_report.jsonl` taxonomy (Genus/Family/Molecule_type/Host). Protein
  standardization lives in `normalization.py` (`standardize_protein(name,
  mol_type, target)`), which dispatches per target: `virus` applies the bundled
  `data/viral_proteins.tsv` map, `bacteria`/`host` do generic cleaning only
  (extension points). All targets share the cleaning pipeline (bracket-tag
  removal, directive stripping, molecular-weight + misspelling normalisation,
  special-char removal, capitalisation, bare-`CDS`/`ORF` → `Unknown`).

The `screening` command orchestrates the steps in this order:

1. **prepare_data.py** `InsertPrefix` — prefix every FASTA header (`>PREFIX/…`).
2. **clean_data.py** `RemoveShortSequences` — drop contigs below `--length`.
3. **make_database.py** `MakeDB` — build BLAST or DIAMOND DBs (`--index_databases`).
4. **similarity_analysis.py** `SimilaritySearch` — blastx / diamond blastx.
5. **filter_table.py** `FilterTable` — filter redundant hits by `qseqid`/range/sense.
6. **bed.py** `GetFasta` — extract putative EE sequences (bedtools).
7. **compare_results.py** `CompareResults` — drop EEs that hit host baits harder.
8. **get_taxonomy.py** `GetTaxonomy` / `GetFinalTaxonomy` / `GetCleanedTaxonomy`
   — join hits to the metadata CSV and build the taxonomy table.
9. **bed.py** `GetAnnotBed` / `MergeBed` / `RemoveAnnotation` — merge truncated
   elements of the same genus/family (`--merge_level`).
10. **clean_data.py** `MaskClean` — optional soft-mask filter (`--clean_masked`).
11. **tag_elements.py** `TagElements` — flag overlapping elements, add
    `Average_pident`.
12. **overlap.py** `FilterOverlap` — unless `--overlap keep` (the default),
    filter elements tagged `overlaped` by the chosen strategy: `longest` (keep
    the longest of each cluster) or `targets` with **exactly one** of a keep-list
    (`--target_families`) or a drop-list (`--non_target_families`), resolved
    per overlap cluster. Removed elements are preserved under `tmp_outputs/`.
    Runs before the GFF3/flanking steps so they see the filtered results.
13. **gff.py** `WriteGFF3` — write the EE taxonomy table as a GFF3 annotation.
14. **get_length.py** `GetLength` + **bed.py** `GetBed`/`BedFlank`/`GetFasta` —
    extract flanking regions (`--flank`).

`utils.py` = path/timing helpers + the `StepInfo`/`RunArguments`/`RunInfo`
dataclasses (and `DownloadArguments`/`SequenceCounts`/`DownloadInfo` for the
`get-databases` log); `versions.py` = dependency-version detection + `env.yml`
comparison (reported in the log and warned about at startup); `log.py` = the
`eefinder` logger + `enable_debug()` (the `--debug` flag on both commands lowers
it to DEBUG; `logger.debug(...)` calls throughout are silent otherwise). The run
finishes by renaming intermediates to `PREFIX.EEs.*` and writing `eefinder.log`
(JSON: `eefinder_version`, `arguments`, `dependencies`, timing, and per-step
info). `get-databases` similarly writes `{outdir}/{prefix}.log` (`DownloadInfo`:
version, arguments, per-phase steps, timing, and a `sequence_counts` block —
`downloaded`/`excluded_uninformative`/`dropped_standardization`/`kept`).

### Inputs / outputs

- **Inputs:** genome FASTA (`-in`), protein DB FASTA (`-db`) + metadata CSV
  (`-mt`, columns `Accession,Species,Genus,Family,Molecule_type,Protein,Host`),
  host-gene baits FASTA (`-bt`).
- **Outputs:** `PREFIX.EEs.fa`, `PREFIX.EEs.tax.tsv`, `PREFIX.EEs.gff3`,
  `PREFIX.EEs.flanks.fa` (+ `.cleaned.*` when `--clean_masked`), plus
  `eefinder.log` and, unless `--removetmp`, a `tmp_files/` archive of
  intermediates. With `--overlap longest|targets`, a `tmp_outputs/` directory
  holds the `PREFIX.EEs.removed.*` elements filtered out of the final results.

## Environment & tooling

External binaries are **required at runtime**: `blastx`/`makeblastdb` (BLAST),
`diamond`, `bedtools` (for `screening`) and `datasets` (NCBI datasets CLI, for
`get-databases`). They are not pip-installable — use conda/micromamba.

```bash
micromamba env create -f env.yml      # or: conda env create -f env.yml
micromamba activate EEfinder
pip install .                         # or `pip install -e .` for development
pip install -r requirements-dev.txt   # pytest + black
```

- Python runtime deps: biopython, pandas (<2), numpy (<2), click. The version
  is read via stdlib `importlib.metadata` (no `pkg_resources`/`setuptools`
  runtime dependency).
- Build metadata is in `setup.py`; `pyproject.toml` holds **only** black +
  pytest config (do not add a `[build-system]` there without also migrating
  `setup.py`).

## Testing

```bash
pytest                      # full suite
pytest -m "not integration" # unit tests only (no external binaries needed)
pytest -m integration       # end-to-end CLI runs against test_files/
```

- Unit tests use synthetic inputs in `tmp_path` — no binaries required.
- Integration tests shell out to the `eefinder` console script and
  auto-**skip** when `blastx`/`makeblastdb`/`bedtools` are absent.
- See `docs/testing.md` for details. `test_files/` holds the example inputs.

## Style

- Format with **black** (line length 88): `black eefinder tests`. The package
  and tests are black-clean; keep them that way.
- Target Python is pinned to **3.9** (`env.yml`). Avoid syntax that only parses
  on 3.10+ (e.g. parenthesised multi-item `with` statements); `py_compile` under
  3.9 should stay green.

## Conventions & gotchas

- **Side-effect classes:** instantiating a step class runs it. Don't expect
  return values; check the output file.
- **Filename chaining:** downstream steps hard-code the accreted suffix of the
  upstream file. Changing an output name means updating every consumer in
  `main.py`.
- **The default `blastx` mode is the reliable path.** The DIAMOND modes can
  fail silently because the subprocess stderr is routed to `DEVNULL`; verify the
  `diamond` build (env pins `diamond=2.2.3`) if a DIAMOND run produces no hits.
- Keep changes minimal and focused; update `CHANGELOG.md` each session.

## Changelog

Record notable changes per session in `CHANGELOG.md` (Keep a Changelog format).
