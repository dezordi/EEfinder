# Acquiring databases (`get-databases`)

`screening` needs three reference inputs:

| Input | CLI flag | Produced by |
|-------|----------|-------------|
| Protein database FASTA | `-db` | `get-databases virus` / `bacteria` |
| Protein metadata CSV | `-mt` | `get-databases virus` / `bacteria` |
| Host-gene bait FASTA | `-bt` | `get-databases host` |

`get-databases` automates the manual NCBI RefSeq downloads (previously done by
hand through the NCBI Virus web UI) using the
[NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/) CLI. It is a **command
group** with one subcommand per database, because each group has its own
defaults and options.

```{note}
`get-databases` requires the `datasets` binary (`ncbi-datasets-cli`) on `PATH`.
Version **≥ 18.1** is required: earlier builds crash on viral downloads because
NCBI added the "acellular root" taxonomy rank above Viruses (taxid 10239) in
2025.
```

## Subcommands

```bash
# Viral protein database + metadata (the -db / -mt inputs)
eefinder get-databases virus -tx Flaviviridae -od db/ -pr virus

# Bacterial protein database + metadata
eefinder get-databases bacteria -tx 2 -od db/ -pr bacteria

# Host proteins used as -bt baits (no metadata CSV)
eefinder get-databases host -tx "Aedes aegypti" -od db/ -pr host
```

### Shared options

| Option | Meaning |
|--------|---------|
| `-tx/--taxon` | NCBI taxon name or tax id (e.g. `Flaviviridae`, `10239`). |
| `-od/--outdir` | Output directory. |
| `-pr/--prefix` | Output basename (default: the dataset type → `virus.fa` / `virus.csv`). |
| `--refseq/--all-sequences` | Restrict to RefSeq (default) or fetch everything. |
| `--cluster/--no-cluster` | Collapse 100%-identical / 100%-coverage duplicate proteins with `cd-hit` before writing the database (on by default). |
| `--debug` | Verbose logging (resolved arguments, the `datasets` command, sequence tallies). |

`-tx` defaults to `10239` (Viruses) for `virus` and `2` (Bacteria) for
`bacteria`, so `eefinder get-databases virus -od db/` downloads the whole RefSeq
viral protein set. For `host`, `-tx` is **required**.

### `virus` / `bacteria`-only options

Both accept `--exclude-uninformative` and `--standardize-proteins` (both on by
default). `host` produces only the baits FASTA and takes neither.

## How the metadata CSV is built

The CSV consumed by `screening -mt` has the columns
`Accession,Species,Genus,Family,Molecule_type,Protein,Host`. `get-databases`
rebuilds it from the downloaded `protein.faa` headers joined with the datasets
`data_report.jsonl` taxonomy:

- `Accession`, `Protein`, `Species` — from the protein FASTA headers.
- `Genus`, `Family` — inferred from the NCBI taxonomy lineage by ICTV name
  suffix (`-viridae` / single-word `-virus`).
- `Host` — from the datasets taxonomy report.
- `Molecule_type` — looked up by family from a bundled
  [ICTV genome-composition table](https://ictv.global/virus-properties)
  (`eefinder/data/`), since the datasets report does not carry it.

## Deduplication (`--cluster`)

Before the metadata CSV is built, `get-databases` collapses **exact duplicate**
proteins with `cd-hit` at 100% identity **and** 100% coverage
(`-c 1.0 -aL 1.0 -aS 1.0`) — so only sequences that are identical over their
entire length are merged to a single representative. This is a lossless
deduplication that shrinks the database (and speeds up the `screening` search)
without discarding any distinct sequence. The metadata CSV is then built from the
deduplicated FASTA, so it describes only the retained representatives.

Clustering is **on by default** for every subcommand (`virus`, `bacteria`,
`host`) and needs `cd-hit` on `PATH`; pass `--no-cluster` to skip it. The number
of duplicates removed is reported in the run log as `clustered_identical`.

## Cleaning and standardisation

### `--exclude-uninformative`

Drops `hypothetical protein` and `uncharacterized protein` records from the
downloaded database (on `virus`/`bacteria`, on by default).

### `--standardize-proteins`

Rewrites the CSV `Protein` column so synonymous names converge (which makes the
final taxonomy table far easier to read and aggregate). Every target shares a
**generic cleaning pass**:

- leaked NCBI `[key=value]` tags (`[organism=…]`, `[gbkey=CDS]`, …) are removed,
  including one level of nesting;
- the special characters `:,/\?!` plus quotes are removed;
- molecular-weight tokens are normalised (`33 kDa`, `33-kDa`,
  `33K-like protein` → `33 kDa protein`);
- common misspellings are fixed so variants converge (`membran` → `membrane`,
  `polyprotien` → `polyprotein`, `glycop` → `glycoprotein`);
- a leading `CDS:` / `ORF:` directive is stripped (records that are **only** a
  bare `CDS`/`ORF` directive are dropped from both the CSV and the FASTA);
- the leading letter is capitalised (`nucleoprotein (N)` → `Nucleoprotein (N)`).

For `virus`, canonical names from the bundled map
(`eefinder/data/viral_proteins.tsv`) are additionally applied, collapsing
synonyms per `Molecule_type` scope — e.g. every RdRp spelling (including compound
names like `P2-RdRp` or `CP/RdRp fusion`) → `RdRp`; the various Capsid spellings
→ `Capsid Protein`. `bacteria` has no name map yet, so it gets the generic
cleaning only; `host` is not standardised. The map is a first draft meant to be
extended — see `eefinder/data/README.md`.

## The `get-databases` run log

Each run writes a JSON log `{outdir}/{prefix}.log` (e.g. `virus.log`) mirroring
the `screening` `eefinder.log` structure: version, resolved arguments, per-phase
timing, and a **`sequence_counts`** block reporting how many sequences were:

| Key | Meaning |
|-----|---------|
| `downloaded` | records fetched from NCBI |
| `excluded_uninformative` | dropped by `--exclude-uninformative` |
| `clustered_identical` | removed as 100%/100% duplicates by `--cluster` |
| `dropped_standardization` | dropped by `--standardize-proteins` (e.g. bare `CDS`, hypothetical) |
| `kept` | records written to the final FASTA/CSV |

## Appendix: manual download

The RefSeq protein database and metadata CSV can still be obtained by hand from
[NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) (select the RefSeq
protein set, then download **Protein** as FASTA and the accompanying **CSV**,
choosing "Download all records"). The bacterial `bac_retriever.py` accessory
script and the manual host-protein selection described in the
[wiki](https://github.com/WallauBioinfo/EEfinder/wiki) remain valid, but
`get-databases` is the recommended, reproducible path.
