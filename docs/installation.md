# Installation

EEfinder is a Python package that drives several external bioinformatics
binaries. Those binaries are **not** pip-installable, so a `pip install` alone is
never a complete installation — they have to come from somewhere else: the
[Bioconda package](#install-from-bioconda), a
[container image](#run-with-a-container) (nothing installed on the host), or a
conda/micromamba environment (the bundled `env.yml`, or the one-liner under
[Install from PyPI](#install-from-pypi)).

## Requirements

| Dependency | Role | Provided by |
|------------|------|-------------|
| Python 3.9 | runtime | `env.yml` |
| BLAST (`blastx`/`blastp`/`makeblastdb`) | similarity search + database build | `env.yml` (`blast`) |
| DIAMOND (`diamond`) | fast alternative to BLAST | `env.yml` (`diamond`) |
| bedtools | sequence extraction / merging | `env.yml` (`bedtools`) |
| NCBI datasets CLI (`datasets`) | database download (`get-databases`) | `env.yml` (`ncbi-datasets-cli`, pinned to 18.36.0) |
| cd-hit | dedup for `--translation_method gv-rv` **and** the `get-databases` `--cluster` step (on by default) | `env.yml` (`cd-hit`) |
| pyrodigal-gv / pyrodigal-rv | protein prediction (`gv`/`rv`/`gv-rv`) | `env.yml` (pip) |
| biopython (<1.86), pandas (<3), numpy (<3), click | Python runtime deps | `env.yml` (pip) |

`env.yml` pins the exact versions EEfinder is developed and tested against. The
Bioconda package and the container image carry all of the above, including the
Python interpreter.

(install-from-bioconda)=
## Install from Bioconda

The [Bioconda package](https://anaconda.org/bioconda/eefinder) is the
recommended route: it declares the external binaries as dependencies, so one
command installs EEfinder *and* the tools it drives.

```bash
conda create -n EEfinder -c conda-forge -c bioconda eefinder
conda activate EEfinder
```

Or into an existing environment:

```bash
conda install -c conda-forge -c bioconda eefinder
```

```{tip}
`mamba` and `micromamba` work as drop-in replacements for `conda` in either
command and resolve considerably faster.
```

(run-with-a-container)=
## Run with a container

The Bioconda package is repackaged automatically as a
[BioContainers image](https://quay.io/repository/biocontainers/eefinder), which
bundles EEfinder with the binaries it drives. Nothing is installed on the host —
useful on machines where you cannot create conda environments, and for
reproducible or HPC runs.

```{important}
BioContainers images are **always version-tagged and there is no `latest`
tag** — a tag has to be given explicitly. Pick one from the
[tag list](https://quay.io/repository/biocontainers/eefinder?tab=tags); the
examples below write it as `<tag>` (of the form
`<version>--<bioconda build string>`, e.g. `2.0.0--pyhdfd78af_0`). The
`get-databases` command needs a tag of version 2.0.0 or later.
```

### Docker

The container sees only what you mount into it, so the genome, the databases and
the output directory all have to live under a mounted path. Mounting the current
directory as the working directory is the simplest arrangement:

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  quay.io/biocontainers/eefinder:<tag> \
  eefinder --version
```

A full run from the repository root, over the bundled `test_files/` (the same
example as in [Running the pipeline](#example-run)):

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  quay.io/biocontainers/eefinder:<tag> \
  eefinder screening \
    -in test_files/Ae_aeg_Aag2_ctg_1913.fasta \
    -od results_test \
    -db test_files/virus_subset.fa \
    -mt test_files/virus_subset.csv \
    -bt test_files/filter_subset.fa \
    -ln 1000 -id -p 2 -lm 100
```

```{note}
`-u "$(id -u):$(id -g)"` makes the results belong to you instead of `root`. It
matters for more than tidiness: `-id` writes the BLAST/DIAMOND index files
*next to the database FASTA*, so that directory must be mounted writable by the
same user — index once outside the container, or keep the databases under the
mounted tree.
```

### Singularity / Apptainer

Pull the image once into a local `.sif` file:

```bash
singularity pull eefinder.sif \
  docker://quay.io/biocontainers/eefinder:<tag>
```

Or download the prebuilt Singularity image from the Galaxy depot, which skips
the Docker-to-SIF conversion:

```bash
curl -Lo eefinder.sif \
  https://depot.galaxyproject.org/singularity/eefinder:<tag>
```

Then run it. Singularity (and Apptainer, its current name — the commands are
interchangeable) bind-mounts the current directory and runs as you, so no `-v`
or `-u` equivalent is needed for data under `$PWD`:

```bash
singularity exec eefinder.sif eefinder --version

singularity exec eefinder.sif eefinder screening \
  -in test_files/Ae_aeg_Aag2_ctg_1913.fasta \
  -od results_test \
  -db test_files/virus_subset.fa \
  -mt test_files/virus_subset.csv \
  -bt test_files/filter_subset.fa \
  -ln 1000 -id -p 2 -lm 100
```

```{tip}
Paths outside the current directory need an explicit bind — databases on a
shared filesystem, for instance: `singularity exec -B /scratch/dbs:/dbs
eefinder.sif eefinder screening -db /dbs/viral.fa ...`.
```

(install-from-pypi)=
## Install from PyPI

The [PyPI package](https://pypi.org/project/eefinder/) ships the Python code
only. Install the external binaries first, then EEfinder itself:

```bash
# 1. the binaries (conda-forge / bioconda)
conda create -n EEfinder -c conda-forge -c bioconda \
  "python>=3.9,<3.12" "blast>=2.5" "diamond>=2.0.15" "bedtools>=2.27" \
  ncbi-datasets-cli cd-hit
conda activate EEfinder

# 2. the package
pip install eefinder
```

This route resolves the binaries to whatever versions are current. To get the
exact versions EEfinder is tested against, install from source with `env.yml` as
below.

## Install from source

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

With a container, prefix the same command — `docker run --rm
quay.io/biocontainers/eefinder:<tag> eefinder --version`, or `singularity exec
eefinder.sif eefinder --version` — and expect the version of the *image tag*,
which can lag the PyPI and Bioconda releases.

## Development install

To run the test suite and format the code, add the development dependencies:

```bash
pip install ".[dev]"                  # pytest + black
# or, with the exact versions CI uses:
pip install -r requirements-dev.txt
```

See the [Developer guide](developer-guide.md) and [Testing](testing.md) pages
for details.

```{tip}
When EEfinder is installed **outside** its source tree, the run log's
dependency-drift check needs to find the reference `env.yml`. Point it there
with `export EEFINDER_ENV_YML=/path/to/env.yml`.
```
