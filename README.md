# EEfinder

EEfinder is a tool/python package that automatizes several tasks related to identification of Endogenous Elements present on Eukaryotic Genomes.

#### Install

EEfinder relies on BLAST, DIAMOND, bedtools and the NCBI datasets CLI, installed through a conda environment (see `env.yml`, which pins BLAST 2.17.0, DIAMOND 2.2.3, bedtools 2.31.1 and ncbi-datasets-cli 18.1.0). Note that BLAST 2.17.0 is not available for macOS on Apple Silicon (osx-arm64); on that platform pin `blast=2.16.0` instead.

```bash
git clone https://github.com/WallauBioinfo/EEfinder.git
cd EEfinder
conda env create -f env.yml
conda activate EEfinder
pip install .
```

#### Check tool

```bash
eefinder --version

#eefinder, version 1.1.1
```

For more information, check [EEfinder Wiki here](https://github.com/WallauBioinfo/EEfinder/wiki)

EEfinder exposes two subcommands:

| Command | Purpose |
|---------|---------|
| `eefinder screening` | Run the EE-finding pipeline on a genome. |
| `eefinder get-databases` | Download the RefSeq protein databases it needs. |

#### Acquiring databases (`get-databases`)

`get-databases` automates the manual NCBI RefSeq downloads described in the wiki
using the [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/) CLI. It is a
command group with one subcommand per database — `virus`, `bacteria` and `host`
— because each group has its own defaults and options. `virus`/`bacteria` fetch
protein sequences and build the metadata CSV
(`Accession,Species,Genus,Family,Molecule_type,Protein,Host`) consumed by
`screening`; `host` fetches the `-bt` baits FASTA only:

```bash
# Viral protein database + metadata (the -db / -mt inputs)
eefinder get-databases virus -tx Flaviviridae -od db/ -pr virus

# Host proteins used as -bt baits (no metadata CSV)
eefinder get-databases host -tx "Aedes aegypti" -od db/ -pr host
```

Options shared by every subcommand:

| Option | Meaning |
|--------|---------|
| `-tx/--taxon` | NCBI taxon name or tax id (e.g. `Flaviviridae`, `10239`). |
| `-od/--outdir` | Output directory. |
| `-pr/--prefix` | Output basename (default: the dataset type → `virus.fa`/`virus.csv`). |
| `--refseq/--all-sequences` | Restrict to RefSeq (default) or fetch everything. |

`virus` and `bacteria` also take `--exclude-uninformative` and
`--standardize-proteins` (see below; `bacteria` gets the generic cleaning only,
since there is no bacterial name map yet). `-tx` defaults to `10239` (Viruses)
for `virus` and `2` (Bacteria) for `bacteria`, so `eefinder get-databases virus
-od db/` downloads the whole RefSeq viral protein set; for `host` it is required.
Every subcommand also accepts `--debug` for verbose logging.

Each run writes a JSON log `{outdir}/{prefix}.log` (e.g. `virus.log`) recording
the arguments, per-phase timing, and a `sequence_counts` block — how many
sequences were `downloaded`, `excluded_uninformative`, `dropped_standardization`
and `kept`. How the metadata CSV columns are filled:

- `Accession`, `Protein`, `Species` — from the protein FASTA headers.
- `Genus`, `Family` — inferred from the NCBI taxonomy lineage by ICTV name
  suffix (`-viridae` / single-word `-virus`).
- `Host` — from the datasets taxonomy report.
- `Molecule_type` — looked up by family from a bundled
  [ICTV genome-composition table](https://ictv.global/virus-properties)
  (`eefinder/data/`), since the datasets report does not carry it.

`--exclude-uninformative` (on `virus`/`bacteria`, on by default) drops
`hypothetical protein` and `uncharacterized protein` records from the downloaded
database.

`--standardize-proteins` (a `virus`/`bacteria` option, on by default) rewrites
the CSV `Protein` column. Every target shares a generic cleaning pass: leaked
NCBI `[key=value]` tags (`[organism=…]`, `[gbkey=CDS]`, …) and the special
characters `:,/\?!` plus quotes are removed; molecular-weight tokens are
normalised (`33 kDa`, `33-kDa`, `33K-like protein` → `33 kDa protein`); common
misspellings are fixed so variants converge (`membran` → `membrane`,
`polyprotien` → `polyprotein`); a leading `CDS:`/`ORF:` directive is stripped
(records that are only a bare `CDS`/`ORF` directive are **dropped** from both the
CSV and the FASTA); and the leading letter is capitalised (`nucleoprotein (N)` →
`Nucleoprotein (N)`). For `virus`, canonical names from the bundled map
(`eefinder/data/viral_proteins.tsv`) are additionally applied, collapsing
synonyms per `Molecule_type` scope (e.g. every RdRp spelling, including compound
names like `P2-RdRp` or `CP/RdRp fusion`, → `RdRp`). `bacteria` has no name map
yet, so it gets the generic cleaning only. The map is a first draft meant to be
extended.

#### Example run (with the bundled `test_files/`)

The repository ships a small example dataset in [`test_files/`](test_files/) so
you can try the full pipeline end-to-end:

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

- `-ln 1000` lowers the minimum contig length (the example contig is shorter
  than the 10000 nt default, which would otherwise filter it out).
- `-id` builds the BLAST/DIAMOND indexes for the databases (needed on the first
  run against a given database).
- Add `-cm` to also emit the mask-cleaned outputs.
- `-an/--analysis` (`virus` default, or `bacteria`) sets the GFF3 feature type
  (`endogenous_viral_element` / `endogenous_bacterial_element`).

##### Resolving overlaping elements (`--overlap` / `-ov`)

Elements that overlap another element assigned to a *different* family are
tagged `overlap_status=overlaped`. `--overlap` decides what to do with them:

| Value | Behaviour |
|-------|-----------|
| `keep` (default) | Keep every element (no filtering). |
| `longest` | Among overlaping elements, keep the longest and drop the shorter ones. |
| `targets` | Resolve each overlap cluster by a family list — a **keep-list** (`--target_families`) or a **drop-list** (`--non_target_families`). |

With `--overlap targets` you must provide **exactly one** of two repeatable
family lists:

- `--target_families` / `-tf` — **keep-list**: in a cluster that contains a
  target-family member, keep the target-family elements and drop the rest; a
  cluster with no target member is kept as-is.
- `--non_target_families` / `-ntf` — **drop-list**: drop the listed families
  from each cluster, unless *every* member of the cluster is a listed family (in
  which case the cluster is kept, never wiped).

```bash
# keep-list: keep only these families in an overlap
eefinder screening ... -ov targets -tf Flaviviridae -tf Caulimoviridae

# drop-list: remove these families from an overlap
eefinder screening ... -ov targets -ntf Retroviridae
```

Filtering is applied to every final result (`.EEs.fa`, `.EEs.tax.tsv`,
`.EEs.gff3`, `.EEs.flanks.fa`, and the `--clean_masked` variants). The
filtered-out elements are **not** discarded — they are written to a
`tmp_outputs/` directory (`PREFIX.EEs.removed.{fa,tax.tsv}`) for inspection.

Results are written to `results_test/`:

```
Ae_aeg_Aag2_ctg_1913.EEs.fa         # EE nucleotide sequences
Ae_aeg_Aag2_ctg_1913.EEs.tax.tsv    # EE taxonomy table
Ae_aeg_Aag2_ctg_1913.EEs.gff3       # EE annotation in GFF3 format
Ae_aeg_Aag2_ctg_1913.EEs.flanks.fa  # EEs + flanking regions
eefinder.log                        # JSON run summary
```

`eefinder.log` also records the EEfinder version and the versions of the
dependencies used (bedtools, BLAST, DIAMOND, python, numpy, pandas), flagging
any that differ from the pins in `env.yml`. When EEfinder is installed outside
its source tree, point it at the reference file with
`export EEFINDER_ENV_YML=/path/to/env.yml`.

> Note: the default `blastx` mode is the tested, reliable path. The DIAMOND
> modes (`-md fast`, etc.) depend on the `diamond` build pinned in `env.yml`.

##### Translation method (`--translation_method` / `-tm`)

By default the similarity search translates the genome in the six reading frames
(`blastx` / `diamond blastx`). `--translation_method` swaps this for up-front
protein prediction, and applies to **both** searches in a run (the main EE
search and the host-bait search) so they never diverge:

| Value | Behaviour |
|-------|-----------|
| `default` | Six-frame `blastx` / `diamond blastx` (current behaviour). |
| `gv` | Predict proteins with [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv), align with `blastp` / `diamond blastp`. |
| `rv` | Predict proteins with [pyrodigal-rv](https://github.com/LanderDC/pyrodigal-rv), align with `blastp` / `diamond blastp`. |
| `gv-rv` | Run both predictors, drop redundancy with `cd-hit` (100%/100%), then align. |

```bash
eefinder screening ... -tm gv
```

The prediction modes map each protein hit's amino-acid coordinates back to
nucleotide coordinates on the source contig, so the output schema (and every
downstream step) is identical to `default`. They need the extra dependencies
`pyrodigal-gv`, `pyrodigal-rv` (pip) and — for `gv-rv` — `cd-hit` (all pinned in
`env.yml`).

#### Tests

EEfinder has a `pytest` suite (unit tests for the data-processing steps plus
end-to-end runs of the CLI against `test_files/`). Install the development
dependencies and run:

```bash
pip install -r requirements-dev.txt
pytest                       # full suite
pytest -m "not integration"  # unit tests only (no external binaries needed)
```

See [docs/testing.md](docs/testing.md) for details.

#### Cite us

If you use EEfinder in your research, please cite https://www.sciencedirect.com/science/article/pii/S2001037024003325:

```
Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEFinder: A general-purpose tool for identification of bacterial and viral endogenized elements in eukaryotic genomes. Computational and Structural Biotechnology Journal. https://doi.org/10.1016/j.csbj.2024.10.012
```