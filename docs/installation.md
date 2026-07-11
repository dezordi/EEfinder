# Installation

EEfinder is a Python package that drives several external bioinformatics
binaries. Those binaries are **not** pip-installable, so the supported route is a
conda/micromamba environment built from the bundled `env.yml`.

## Requirements

| Dependency | Role | Provided by |
|------------|------|-------------|
| Python 3.9 | runtime | `env.yml` |
| BLAST (`blastx`/`blastp`/`makeblastdb`) | similarity search + database build | `env.yml` (`blast`) |
| DIAMOND (`diamond`) | fast alternative to BLAST | `env.yml` (`diamond`) |
| bedtools | sequence extraction / merging | `env.yml` (`bedtools`) |
| NCBI datasets CLI (`datasets`) | database download (`get-databases`) | `env.yml` (`ncbi-datasets-cli`) |
| cd-hit | dedup for `--translation_method gv-rv` | `env.yml` (`cd-hit`) |
| pyrodigal-gv / pyrodigal-rv | protein prediction (`gv`/`rv`/`gv-rv`) | `env.yml` (pip) |
| biopython, pandas (<2), numpy (<2), click | Python runtime deps | `env.yml` (pip) |

```{note}
BLAST 2.17.0 is available on `linux-64` but **not** on macOS Apple Silicon
(`osx-arm64`), which tops out at 2.16.0. On Apple Silicon either pin
`blast=2.16.0` in `env.yml` or build on a `linux-64` / `osx-64` machine.
```

## Install with micromamba (recommended)

```bash
git clone https://github.com/WallauBioinfo/EEfinder.git
cd EEfinder

micromamba env create -f env.yml     # or: conda env create -f env.yml
micromamba activate EEfinder

pip install .                        # or `pip install -e .` for development
```

## Verify the installation

```bash
eefinder --version
# eefinder, version 2.0.0

eefinder --help
# Usage: eefinder [OPTIONS] COMMAND [ARGS]...
#   screening       Run the EEfinder screening pipeline on a genome.
#   get-databases   Download RefSeq protein databases (and metadata) ...
```

## Development install

To run the test suite and format the code, add the development dependencies:

```bash
pip install -r requirements-dev.txt   # pytest + black
```

See the [Developer guide](developer-guide.md) and [Testing](testing.md) pages
for details.

```{tip}
When EEfinder is installed **outside** its source tree, the run log's
dependency-drift check needs to find the reference `env.yml`. Point it there
with `export EEFINDER_ENV_YML=/path/to/env.yml`.
```
