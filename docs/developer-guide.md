# Developer guide

## Package layout

EEfinder uses a flat layout — the package lives in `eefinder/`, no `src/`. Each
processing step is a small **side-effect class** whose `__init__` runs the work
(files in, files out); there is no shared in-memory pipeline object. Steps
communicate through files on disk whose names accrete suffixes (`.rn`, `.fmt`,
`.blastx`, `.filtred`, `.bed`, `.tax`, …).

| Module | Responsibility |
|--------|----------------|
| `scripts/main.py` | The `click` command group (`screening`, `get-databases`). |
| `prepare_data.py` | `InsertPrefix` — prefix FASTA headers. |
| `clean_data.py` | `RemoveShortSequences`, `MaskClean`. |
| `make_database.py` | `MakeDB` — build BLAST/DIAMOND databases. |
| `similarity_analysis.py` | `SimilaritySearch` — the search (run twice). |
| `translation.py` | Protein prediction + coordinate traceback (`gv`/`rv`/`gv-rv`). |
| `filter_table.py` | `FilterTable` — redundant-hit filtering. |
| `bed.py` | `GetFasta`, `MergeBed`, `BedFlank`, … (bedtools wrappers). |
| `compare_results.py` | `CompareResults` — host-bait filtering. |
| `get_taxonomy.py` | `GetTaxonomy` / `GetFinalTaxonomy` / `GetCleanedTaxonomy`. |
| `tag_elements.py` | `TagElements` — overlap flags + `Average_pident`. |
| `overlap.py` | `FilterOverlap` — overlap resolution strategies. |
| `gff.py` | `WriteGFF3`. |
| `get_databases.py` | `GetDatabases` — the `get-databases` implementation. |
| `normalization.py` | `standardize_protein` — per-target protein-name cleaning. |
| `utils.py` | Path/timing helpers + the run-info dataclasses. |
| `versions.py` | Dependency-version detection + `env.yml` comparison. |
| `log.py` | The `eefinder` logger + `enable_debug()`. |

## Conventions & gotchas

- **Side-effect classes:** instantiating a step class runs it. Don't expect
  return values — check the output file.
- **Filename chaining:** downstream steps hard-code the accreted suffix of the
  upstream file. Changing an output name means updating every consumer in
  `main.py`.
- **The default `blastx` mode is the reliable path.** The DIAMOND modes can fail
  silently because the subprocess stderr is routed to `DEVNULL`; verify the
  `diamond` build (env pins `diamond=2.2.3`) if a DIAMOND run produces no hits.
- **Debug logging:** `--debug` (on both commands) lowers the `eefinder` logger to
  DEBUG via `log.enable_debug()`; the `logger.debug(...)` calls throughout are
  silent otherwise.

## Style

- Format with **black** (line length 88): `black eefinder tests`. The package and
  tests are black-clean; keep them that way.
- Target Python is pinned to **3.9** (`env.yml`). Avoid syntax that only parses
  on 3.10+ (e.g. parenthesised multi-item `with` statements); use
  `from __future__ import annotations` for newer typing syntax.
- Public functions/classes get NumPy-style docstrings.

## Build metadata

Build metadata lives entirely in `pyproject.toml`, built with the **hatchling**
backend. The `[project]` table declares the runtime dependencies, the `dev`
extra (`pytest` + `black`), and the `eefinder` console script; the same file
holds the black and pytest configuration. There is no `setup.py` or
`MANIFEST.in`. The version is set once in `[project].version` and read back at
runtime via stdlib `importlib.metadata` (no `pkg_resources` / `setuptools`
runtime dependency).

```bash
pip install .            # runtime install
pip install -e ".[dev]"  # editable install with pytest + black
python -m build          # build the wheel/sdist (needs the `build` package)
```

## Continuous integration

`.github/workflows/tests.yml` runs the binary-free unit tests
(`pytest -m "not integration"`) on every pull request to `master`. It installs
only the pip runtime dependencies (kept in sync with `env.yml`) and pytest — no
BLAST/DIAMOND/bedtools — so the tool-gated tests skip cleanly. See
[Testing](testing.md) for the full suite.

## Building the docs locally

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
# open docs/_build/html/index.html
```
