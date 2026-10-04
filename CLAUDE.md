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
is organised at **two levels**:

- **Steps** (package root: `prepare_data.py`, `filter_table.py`, `bed.py`, ...)
  are small classes whose `__init__` runs the work as a side effect (files in,
  files out) and return nothing. Each names its output by appending a suffix to
  its input (`.filtred`, `.merge`, `.fmt`, ...). They are the implementation.
- **Stages** (`eefinder/stages/`) are the public layer: one class per pipeline
  stage, taking explicit input paths, running on `.run()`, and returning a
  dataclass naming every file produced. `eefinder/stages/paths.py`
  (`ScreeningPaths`) is the single source of truth for every file name — stages
  rename each step's suffixed output to the canonical name, so the stage
  interface does not depend on step naming. **New code composes stages.**

`models.py` holds the output schemas (`HIT_TABLE`, `FILTERED_HIT_TABLE`,
`ELEMENT_TABLE`, ...). A stage declares `models = {output_attr: model}` and
`Stage.run` validates each one after `_execute`, so a stage enforces its
guarantees rather than documenting them.

`search_methods.py` makes hit production pluggable: a `SearchMethod` writes a
table satisfying `HIT_TABLE`, and `SearchMethod.run` checks it.
`resolve_search_method(mode, translation_method)` maps the two CLI options onto
the registry (`blastx`, `diamond`, `predicted`); `SequenceAlignment` and
`PutativeElementsFilter` accept a `search_method=` override. Adding a method
means subclassing, registering, and satisfying `HIT_TABLE` — nothing downstream
changes.

The CLI entry point is the `cli` group in `eefinder/scripts/main.py`, which
holds only the group itself and registers the two command groups — each in its
own module under `eefinder/scripts/`:

- **`screening`** (`eefinder/scripts/screening.py`) — a `click` **group** of
  nine subcommands, one per stage plus `all`. It uses a `DefaultGroup` subclass
  whose `parse_args` prepends `all` when the first argument is not a
  subcommand, so `eefinder screening <options>` keeps working. The subcommands
  are `prepare`, `clean`, `align`, `filter`, `taxonomy`, `merge`,
  `postprocess`, `flanks`, `all`; options are shared through the
  `_outdir_options`/`_database_options`/`_search_options`/
  `_postprocess_options`/`_common_options` decorators. A single stage defaults
  its file inputs to the `ScreeningPaths` names for `-od`/`-pr`. Only `all`
  renames results to `PREFIX.EEs.*`, archives intermediates and writes
  `eefinder.log`. Every subcommand takes `--versions-yml`/`--versions-key`
  (`versions.write_versions_yml`) and `--debug`.
- **`get-databases`** (`eefinder/scripts/get_databases.py`) — a `click`
  **group** that downloads the RefSeq protein databases (the implementation is
  `eefinder/get_databases.py` `GetDatabases`; note the two modules share a name
  in different packages) via the NCBI `datasets` CLI. It
  has one subcommand per database, each with group-specific defaults/options:
  `virus` (taxon default `10239`, `--exclude-uninformative` +
  `--standardize-proteins`) and `bacteria` (taxon default `2`,
  `--exclude-uninformative`, `--standardize-proteins` = generic cleaning only,
  no name map) produce a protein FASTA + metadata CSV (the
  `-db`/`-mt` inputs); `host` (taxon **required**) produces the `-bt` baits
  FASTA. The shared `-od`/`-pr`/`--refseq` options come from the
  `_common_download_options` decorator in `scripts/get_databases.py`;
  `_run_get_databases` does
  the `datasets`-binary check and calls `GetDatabases`. The metadata CSV is
  rebuilt from the `protein.faa` headers (Accession/Protein/Species) joined with
  the `data_report.jsonl` taxonomy (Genus/Family/Molecule_type/Host). Protein
  standardization lives in `normalization.py` (`standardize_protein(name,
  mol_type, target)`), which dispatches per target: `virus` applies the bundled
  `data/viral_proteins.tsv` map, `bacteria`/`host` do generic cleaning only
  (extension points). All targets share the cleaning pipeline (bracket-tag
  removal, directive stripping, molecular-weight + misspelling normalisation,
  special-char removal, capitalisation, bare-`CDS`/`ORF` → `Unknown`).
  `--exclude-taxon` (repeatable; defaults to `2697049`/SARS-CoV-2 on `virus`,
  cleared with `--exclude-taxon none`) leaves a branch out **at fetch time**:
  `taxon_exclusion.py` `expand_taxon_excluding` walks the lineage from the requested
  taxon to the excluded one and keeps every child except the one on the path, so
  the branch is never requested. The resulting taxa go to `datasets` via
  `--inputfile` (max 100, so `batch_taxa` splits longer lists into several
  downloads, each extracted into its own `part_N/` and merged by
  `find_data_reports`). Records attached directly to a rank on the path are
  unreachable this way (measured: 0.005%).


`progress.py` = terminal progress reporting: it mirrors a subprocess's own
progress display (the `datasets` CLI draws one) instead of capturing it, wraps
`click.progressbar` for the in-process loops, and implements the download
retry/stall detection (`run_with_retries`; a transfer is "stalled" only when
neither new output nor growth of the output file is seen within the timeout).
Everything degrades to the previous silent behaviour off a terminal or under
`EEFINDER_NO_PROGRESS=1`.

`utils.py` = `-mt` header validation (`check_metadata_columns` /
`check_metadata_file` — the metadata CSV is read by column position downstream,
so a missing column is an error at CLI start-up and extra/reordered columns are
warned about and fixed in memory) + path/timing helpers + the
`StepInfo`/`RunArguments`/`RunInfo` dataclasses (and `DownloadArguments`/`SequenceCounts`/`DownloadInfo` for the
`get-databases` log); `versions.py` = dependency-version detection + `env.yml`
comparison, plus `report_run_context` (the startup banner, called by the
`screening all` subcommand); `log.py` = the
`eefinder` logger + `enable_debug()` (the `--debug` flag on both commands lowers
it to DEBUG; `logger.debug(...)` calls throughout are silent otherwise). The run
finishes by renaming intermediates to `PREFIX.EEs.*` and writing `eefinder.log`
(JSON: `eefinder_version`, `arguments`, `dependencies`, timing, and per-step
info). `--released-before` bounds a build to a date for reproducibility: native in the
`datasets` client for bacteria/host, applied by EEfinder per **organism** for
viruses (their subcommand has no such flag, and proteins cannot be traced to a
genome record). `get-databases` also writes `{prefix}.tracking.tsv` (one row per downloaded
accession: product before/after standardisation, whether it was removed and why,
and the cd-hit cluster + representative for duplicates, plus each organism's
earliest release date) and `{prefix}.clstr`,
and deletes the downloaded zip and the extracted `*_ncbi/` directory unless
`--keep-download` is passed. It similarly writes `{outdir}/{prefix}.log` (`DownloadInfo`:
version, arguments, per-phase steps, timing, and a `sequence_counts` block —
`downloaded`/`excluded_uninformative`/`clustered_identical`/
`dropped_standardization`/`kept`). Unless `--no-cluster`, `get-databases` runs a
`cd-hit` 100%-identity/100%-coverage dedup on the protein FASTA
(`cluster_identical_proteins`, reusing `translation.cluster_proteins`) before
building the metadata CSV, so the CSV only describes the retained
representatives.

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

External binaries are **required at runtime**: `blastx`/`blastp`/`makeblastdb`
(BLAST), `diamond`, `bedtools` (for `screening`) and `datasets` (NCBI datasets
CLI, for `get-databases`); `cd-hit` is needed for `--translation_method gv-rv`
and for the `get-databases` `--cluster` step (100%/100% duplicate collapse, on
by default). They are not pip-installable — use conda/micromamba.
The `gv`/`rv`/`gv-rv` methods also need the pip packages `pyrodigal-gv` and
`pyrodigal-rv` (pinned in `env.yml`).

```bash
micromamba env create -f env.yml      # or: conda env create -f env.yml
micromamba activate EEfinder
pip install .                         # or `pip install -e .` for development
pip install ".[dev]"                  # + pytest + black (or requirements-dev.txt)
```

- Python runtime deps (declared in `pyproject.toml`): click, biopython
  (**`<1.86`** — `Bio.Blast.Applications` was removed there and
  `make_database.py`/`similarity_analysis.py` use it), pandas (<3), numpy (<3). `pip install .` pulls them; the external binaries
  still come from `env.yml`.
- Build metadata lives in `pyproject.toml` (**hatchling** backend, `[project]`
  table with runtime deps + `dev` extra + the `eefinder` console script). There
  is no `setup.py` / `MANIFEST.in`. The same file also holds the black and pytest
  config. The version is set in `[project].version`; `eefinder/__init__.py` reads
  it back at runtime via stdlib `importlib.metadata` (no `pkg_resources` /
  `setuptools` runtime dependency).

## Testing

```bash
pytest                      # full suite
pytest -m "not integration" # unit tests only (no external binaries needed)
pytest -m integration       # end-to-end CLI runs against test_files/
```

- Unit tests use synthetic inputs in `tmp_path` — no binaries required.
  `test_models.py`, `test_search_methods.py` and `test_stages.py` cover the
  output models, the search-method registry and the stage base class.
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

- **Side-effect step classes:** instantiating a *step* class runs it and
  returns nothing. A *stage* returns its outputs from `.run()` — use stages.
- **Filename chaining is confined to the step classes.** Stages rename their
  outputs to `ScreeningPaths` names, so a new file name is one edit in
  `stages/paths.py`. Never hard-code an intermediate name elsewhere.
- **The `PREFIX.EEs.*` names are published** and must not be renamed.
- **The default `blastx` mode is the reliable path.** The DIAMOND modes can
  fail silently because the subprocess stderr is routed to `DEVNULL`; verify the
  `diamond` build (env pins `diamond=2.0.15`) if a DIAMOND run produces no hits.
- **The integration golden files are the safety net.** Any refactor of the
  pipeline must leave `PREFIX.EEs.{fa,tax.tsv,gff3,flanks.fa}` byte-identical;
  `pytest -m integration` asserts it. Running `all` and running the eight
  stages in sequence must agree.
- Keep changes minimal and focused; update `CHANGELOG.md` each session.

## Changelog

Record notable changes per session in a local `CHANGELOG.md`.
