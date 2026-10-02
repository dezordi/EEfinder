# EEfinder

EEfinder is a Python CLI and package that automates the identification of
**Endogenous Elements (EEs)** — virus- or bacteria-derived sequences integrated
into eukaryotic genomes — by combining similarity search with genomic-junction
reasoning.

#### Install

EEfinder drives BLAST, DIAMOND, bedtools, the NCBI `datasets` CLI and cd-hit,
which are not pip-installable, so they have to come from conda or a container.

##### From Bioconda

The shortest route — the package declares the external binaries, so one command
installs EEfinder and everything it drives.

```bash
conda create -n EEfinder -c conda-forge -c bioconda eefinder
conda activate EEfinder
```

##### From PyPI

Install the binaries first, then the package:

```bash
conda create -n EEfinder -c conda-forge -c bioconda \
  "python>=3.9,<3.12" "blast>=2.5" "diamond>=2.0.15" "bedtools>=2.27" \
  ncbi-datasets-cli cd-hit
conda activate EEfinder

pip install eefinder
```

##### From source

Cloning also gives you `env.yml`, which pins the exact versions EEfinder is
developed and tested against (BLAST 2.17.0, DIAMOND 2.2.3, bedtools 2.31.1,
ncbi-datasets-cli 18.36.0, cd-hit 4.8.1), plus the example data in `test_files/`.

```bash
git clone https://github.com/WallauBioinfo/EEfinder.git
cd EEfinder
conda env create -f env.yml
conda activate EEfinder
pip install .
```

##### With a container

The Bioconda package is repackaged as a
[BioContainers image](https://quay.io/repository/biocontainers/eefinder), so
nothing is installed on the host. The images are always version-tagged (there is
no `latest`) — pick a tag from the repository and see
[the installation docs](https://eefinder.readthedocs.io/en/latest/installation.html)
for the Docker and Singularity/Apptainer invocations.

#### Check tool

```bash
eefinder --version

#eefinder, version 2.0.0
```

#### Usage

EEfinder is a command group with two commands — `get-databases` to acquire the
reference data and `screening` to run the pipeline:

```bash
# download the reference databases
eefinder get-databases virus -od db/ -pr virus
eefinder get-databases host -tx "Aedes aegypti" -od db/ -pr host

# run the pipeline
eefinder screening -in genome.fasta -od results/ \
  -db db/virus.fa -mt db/virus.csv -bt db/host.fa
```

#### Documentation

Full documentation is hosted on Read the Docs: **https://eefinder.readthedocs.io**

#### Cite us

If you use EEfinder in your research, please cite https://www.sciencedirect.com/science/article/pii/S2001037024003325:

```
Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEFinder: A general-purpose tool for identification of bacterial and viral endogenized elements in eukaryotic genomes. Computational and Structural Biotechnology Journal. https://doi.org/10.1016/j.csbj.2024.10.012
```
