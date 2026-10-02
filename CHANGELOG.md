# Changelog

All notable changes to EEfinder are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/).


## [2.0.0]

Major release, relative to 1.1.3. The CLI became a command group and gained a
command that builds its own reference databases; the similarity search gained
alternative translation methods and overlap resolution; and the project gained
GFF3 output, an auditable run log, a test suite and CI.

Each feature is documented in full on
[the documentation site](https://eefinder.readthedocs.io) — this section records
what changed, not how each option works.

### Added
- **`get-databases`, a command that builds the reference inputs** (`-db`, `-mt`,
  `-bt`) from NCBI RefSeq via the `datasets` CLI, replacing the manual downloads
  through the NCBI Virus web UI. One subcommand per database (`virus`,
  `bacteria`, `host`), each producing the protein FASTA and — for `virus` and
  `bacteria` — the matching metadata CSV, rebuilt from the downloaded protein
  headers joined with the taxonomy report. A build can collapse exact duplicate
  proteins, drop uninformative records, rewrite synonymous protein names to
  canonical ones, leave a branch of the taxonomy out of the download entirely,
  and be bounded to a release date so it can be reproduced later. Each run writes
  a JSON log with its sequence tallies and a tracking table recording the fate of
  every downloaded accession, including the records that never reached the
  database. Whole-RefSeq transfers are multi-GB and fail intermittently, so
  downloads report progress and recover from reported errors, corrupt archives
  and silent hangs on their own. See
  [Acquiring databases](docs/get-databases.md).

  > **Note on a default that drops data:** on `virus`, SARS-CoV-2 is left out of
  > the download by default — it is 61% of all viral records in GenBank and of no
  > interest to an endogenous-element search. Pass `--exclude-taxon none` to
  > include it.

- **Alternative translation methods for the similarity search**
  (`-tm/--translation_method`). Besides the six-frame `blastx` default, the
  genome can be searched through proteins predicted up front with pyrodigal-gv
  and/or pyrodigal-rv, aligned with `blastp`. Amino-acid coordinates are traced
  back to contig nucleotides, so the output files have the same shape in every
  mode. The method applies to both searches in a run. See
  [Translation methods](docs/translation-methods.md).
- **Overlap resolution** (`-ov/--overlap`). Elements that overlap while being
  assigned to different families could previously only be reported as such;
  they can now be resolved by keeping the longest of each cluster, or by a
  keep-list or drop-list of families. Filtered-out elements are preserved under
  `tmp_outputs/` rather than discarded. See
  [Overlap resolution](docs/overlap.md).
- **GFF3 annotation output** — `PREFIX.EEs.gff3`, the taxonomy table as a
  genome-browser-loadable annotation, with `-an/--analysis` selecting the feature
  type (`endogenous_viral_element` / `endogenous_bacterial_element`).
- **An auditable run log.** `eefinder.log` records the resolved arguments, the
  parameters used at each step, per-step timing, and the detected versions of
  every dependency flagged against the `env.yml` pins, with a startup warning on
  drift. `--debug` adds verbose tracing, and `kill -USR1 <pid>` dumps a traceback
  of every thread so a run that appears stuck can be diagnosed.
- **A `pytest` suite and CI.** Unit tests for every data-processing step (no
  external binaries needed) plus scenario-driven end-to-end runs against
  `test_files/`, including byte-for-byte golden comparisons of the main outputs.
  A GitHub Actions workflow runs the binary-free tests on every pull request.
  See [Testing](docs/testing.md).
- **`accessory_scripts/update_ictv_table.py`** regenerates the bundled ICTV
  genome-composition table (the source of `Molecule_type`) from the current ICTV
  Virus Metadata Resource, reporting what a new Master Species List release would
  change and flagging the families whose rename would silently empty
  `Molecule_type`. The table's format and update procedure are now documented in
  `eefinder/data/README.md`, which previously only recorded where the file came
  from — and pointed at an ICTV page that has no downloadable table.
- **Documentation pages for all of the above**, plus a `CHANGELOG.md` (this
  file). The Read-the-Docs site itself arrived in 1.1.2.

### Changed
- **The CLI is a command group and the pipeline moved under `screening`.** It is
  invoked as **`eefinder screening <options>`** instead of `eefinder <options>`;
  the options themselves are unchanged. **(Breaking.)**
- **`--merge_level` now defaults to `family`** instead of `genus`, so truncated
  fragments of one integration are merged more readily by default.
  **(Breaking: changes results.)**
- **The genome is prepared in a single pass.** Prefixing the headers and dropping
  short contigs were two steps chained through an intermediate file, writing the
  whole genome to disk twice — 3.2 GB of intermediates for a 1.6 GB genome, the
  first copy never read again. The output is byte-for-byte what the two steps
  produced.
- **The README install section is reinstated.** 1.1.3 removed it in favour of the
  documentation site; it is back because `README.md` is the PyPI landing page,
  and now also covers the Bioconda and container routes.
- **The installation docs were brought forward to the command group** — the
  container examples run `eefinder screening …`, and the Bioconda recipe note
  lists the new binary dependencies.
- **Packaging.** The console script points at the command group; `pyproject.toml`
  gains the `dev` extra (pytest + black) and the black and pytest configuration.
  The dependency bounds are 1.1.2's, unchanged, and verified against
  pandas 2.3 / numpy 2.2.
- **The pipeline was refactored** for type hints, NumPy-style docstrings and
  dataclasses for the run-log structures, without changing the default-path
  outputs — asserted byte-for-byte by the golden integration runs.

### Removed
- **The top-level pipeline invocation `eefinder <options>`.** In 1.1.3 the entry
  point was a single command that ran the pipeline directly; it is now the
  `screening` subcommand of the group. **(Breaking.)**

### Fixed
- **`--clean_masked` produced an empty cleaned taxonomy table.** The cleaned
  FASTA record IDs keep the `PREFIX/` that the taxonomy table's `Element-ID`
  drops, so the id comparison never matched and `*.EEs.cleaned.tax.tsv` came out
  with only a header. IDs are now compared with the prefix removed.
- **The run log reported the wrong `merge_level`.** It was read from the same
  argument position as `length`, so it always echoed `--length`.
- **The run log's end timestamps were malformed** by a stray `%` in the format
  string (`%Y-%m-%d %H:%M%:%S`).

## [1.1.3]

Released from the official repository
([PR #43](https://github.com/WallauBioinfo/EEfinder/pull/43)) while 2.0.0 was in
development, and the release 2.0.0 builds on.

### Added
- **The `-mt` metadata header is validated** (`utils.check_metadata_columns` /
  `utils.check_metadata_file`). The taxonomy steps address the joined metadata
  **by column position** (`_SPECIES_COL` … `_HOST_COL` in `get_taxonomy.py`), so
  an extra, missing or reordered column in `-mt` used to produce a silently wrong
  taxonomy table. `check_metadata_file` reads the header alone at CLI start-up,
  so a bad file fails before the similarity search rather than hours into the
  run, and `GetTaxonomy` reselects the seven expected columns in the canonical
  order before merging. A missing column is a `ValueError` naming it; extra
  columns and a different column order are warnings — the extras are dropped and
  the order fixed in memory, and the file on disk is never rewritten. Verified
  end to end on the 2.0.0 tree: a shuffled eight-column metadata file produces a
  taxonomy table byte-for-byte identical to the canonical one.
- **Bioconda and container installation documented** (`docs/installation.md`).
  `bioconda::eefinder` installs EEfinder together with the binaries it drives,
  and the BioContainers image needs nothing on the host; Docker and
  Singularity/Apptainer invocations are given, including the mount/user caveat
  around `-id` (the index files are written next to the database FASTA). The
  developer guide gains a section on the Bioconda recipe, which lives in
  `bioconda-recipes`, not here: the autobump bot follows PyPI versions, but
  dependency changes have to be carried over by hand.

### Removed
- **The install instructions were dropped from `README.md`** in favour of the
  Read-the-Docs page. 2.0.0 reinstates them (see above).

## [1.1.2]

Packaging and documentation release ([PR #40](https://github.com/WallauBioinfo/EEfinder/pull/40)).
No change to the pipeline or its results.

### Added
- **A Read-the-Docs documentation site** (`docs/`, `.readthedocs.yaml`) built with
  Sphinx + MyST and the `sphinx_rtd_theme`, replacing the GitHub wiki as the
  reference documentation: installation, databases, running, outputs, custom
  arguments, accessory scripts and a developer guide.
- **A `LICENSE` file** (MIT), with `license`/`license-files` declared in the
  package metadata.
- **PyPI publishing** — `eefinder` is installable with `pip install eefinder`,
  released by a `publish.yml` GitHub Actions workflow using Trusted Publishing.

### Changed
- **The build moved from `setup.py` to `pyproject.toml`** with the **hatchling**
  backend. `setup.py` was removed. The `[project]` table carries the runtime
  dependencies (`click`, `biopython>=1.79,<1.86`, `pandas<3`, `numpy<3`),
  `requires-python = ">=3.9,<3.12"`, maintainers, URLs, classifiers and the sdist
  manifest.
- **`biopython` is capped at `<1.86`**, which removed `Bio.Blast.Applications` —
  still used by `make_database.py` and `similarity_analysis.py`.
- **The version is read with stdlib `importlib.metadata`** instead of
  `pkg_resources`, dropping the runtime dependency on `setuptools`.
- The README points at the documentation site instead of the wiki, and documents
  the PyPI and from-source installation routes.

## [1.1.1]

Small fix release ([PR #39](https://github.com/WallauBioinfo/EEfinder/pull/39)).

### Fixed
- **EEfinder logged on the root logger.** `logging.getLogger("")` returns the
  root logger, so the messages carried no identifiable name and could not be
  configured or silenced independently of everything else in the process. The
  logger is now the named `eefinder` one.
- **The closing message still said the study was unpublished.** It now prints the
  *Computational and Structural Biotechnology Journal* citation
  ([doi:10.1016/j.csbj.2024.10.012](https://doi.org/10.1016/j.csbj.2024.10.012)).

## [1.1.0]

Output and naming release ([PR #38](https://github.com/WallauBioinfo/EEfinder/pull/38)).

### Added
- **`Average_pident` in `PREFIX.EEs.tax.tsv`** — the mean percent identity of the
  hits supporting each element. The per-hit identity is carried through the merge
  step (appended to each accession in `Protein-IDs` as `accession|pident`) and
  averaged per element, so the strength of the evidence behind an assignment is
  visible in the output table rather than only in the intermediates.

### Changed
- **`-bt` was renamed `--baits` → `--hostgenesbaits`**, to say what the file must
  contain: the **host's own** genes, not an arbitrary filter set. The short flag
  is unchanged. **(Breaking for anyone using the long form.)**
- The README was reduced to install/usage and pointed at the GitHub wiki; the
  frozen dependency table (every transitive conda package) was dropped.

### Fixed
- **pandas `FutureWarning`s** from chained-assignment `fillna(..., inplace=True)`
  and from assigning a string into an inferred-dtype column (`sense`), which
  would have become errors in pandas 3.
- Removed a large commented-out block in `filter_table.py` that duplicated the
  redundant-hit logic; the explanation now lives in
  [Custom arguments](docs/custom-arguments.md).

## [1.0.0]

Documentation and typing release ([PR #34](https://github.com/WallauBioinfo/EEfinder/pull/34)).
The version reaching 1.0.0 went with the README moving from "Alpha dev version"
to "Beta dev version"; the pipeline and its outputs are those of 0.3.1.

### Added
- **Docstrings and type hints across every module.** Each step class gained a
  class-level docstring naming its inputs and the CLI argument each one comes
  from (the `Keyword arguments:` style; the NumPy-style rewrite came in 2.0.0),
  and every method gained parameter and return annotations.

### Changed
- Docstrings moved from the methods to the classes, so the documented contract is
  visible where the class is instantiated (instantiating a step **runs** it).
- The codebase was formatted with black.

### Fixed
- **Three numeric arguments were interpolated into shell commands unvalidated.**
  `bedtools merge -d` (`--limit`), `diamond makedb --threads` and
  `diamond blastx -p` (`--threads`) are now cast with `int()` before being
  formatted into the command line, so a value arriving as a string cannot
  produce a malformed command.

[2.0.0]: https://github.com/WallauBioinfo/EEfinder/compare/v1.1.3...v2.0.0
[1.1.3]: https://github.com/WallauBioinfo/EEfinder/compare/v1.1.2...v1.1.3
[1.1.2]: https://github.com/WallauBioinfo/EEfinder/compare/v1.1.1...v1.1.2
[1.1.1]: https://github.com/WallauBioinfo/EEfinder/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/WallauBioinfo/EEfinder/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/WallauBioinfo/EEfinder/compare/v0.3.1...v1.0.0
