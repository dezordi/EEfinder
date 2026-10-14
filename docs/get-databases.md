# Acquiring databases (`get-databases`)

`screening` needs three reference inputs:

| Input | CLI flag | Produced by |
|-------|----------|-------------|
| Protein database FASTA | `-db` | `get-databases virus` / `bacteria` |
| Protein metadata CSV | `-mt` | `get-databases virus` / `bacteria` |
| Host-gene bait FASTA | `-bt` | `get-databases host` |

`get-databases` automates the manual NCBI RefSeq downloads (previously done by
hand through the NCBI Virus web UI) using the
[NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/) CLI.

```{note}
`get-databases` requires the `datasets` binary (`ncbi-datasets-cli`) on `PATH`.
Version **≥ 18.1** is required: earlier builds crash on viral downloads because
NCBI added the "acellular root" taxonomy rank above Viruses (taxid 10239) in
2025. `env.yml` pins **18.36.0**.
```

(databases-define-search)=
## The databases define the search space

These three files are not just inputs — they are what makes the run a search for
one thing rather than another. EEfinder has no built-in reference data and no
notion of which taxa are "interesting": it reports similarity to whatever proteins
you give it, labelled with whatever taxonomy your metadata table carries.

Three consequences are worth stating explicitly:

- **The donor lineage is your choice.** A RefSeq viral protein set screens for
  endogenous viral elements; a bacterial protein set screens for endogenous
  bacterial elements. Nothing in the pipeline is specific to either, so any
  donor lineage whose proteins you can assemble into a FASTA can be screened
  for on the same code path.
- **Scope follows the database.** Restricting `-db` to a single family — or a
  single genus — restricts the whole analysis to it, at a fraction of the
  runtime. Conversely, a hit can only ever be reported for a lineage that is
  represented in the database: absence of an element is evidence of absence
  *from your reference set*, not from the genome.
- **The taxonomy is yours too.** The `taxonomy` lineage, `molecule_type` and
  `host` in the output tables come from the `Taxonomy`, `Molecule_type` and
  `Host` columns of `-mt`. EEfinder never consults an external taxonomy: an
  element reported from a single reference carries that reference's lineage
  verbatim, and one merged from several carries their lowest common ancestor
  (see [the lineage of a merged element](output.md#the-lineage-of-a-merged-element)).
  The rank names, spellings and any `Unk` ranks in your results are therefore
  the ones in your table (see [the `Taxonomy` column](#taxonomy-column)).

The `-bt` bait set defines the opposite side of the comparison: it is the
background that candidate elements must *fail* to match better than they match
`-db`. It should therefore represent the host's own proteome (see
[Host-gene baits](#host-gene-baits)).

(host-gene-baits)=
## Host-gene baits (`-bt`)

The bait set is a collection of **host proteins** used to discard false
positives: a putative EE that matches a host gene with a higher bitscore than it
matches the viral/bacterial reference is removed (see
[Running the pipeline](screening.md)).

Use the NCBI proteins of the host species, or of the closest available
relative. No metadata table is needed for this input — only the FASTA.
`test_files/filter_subset.fa` shows the expected shape (RefSeq *Aedes
albopictus* proteins used against an *Aedes aegypti* genome).

## Subcommands

```bash
# Viral protein database + metadata (the -db / -mt inputs)
eefinder get-databases virus -od db/ -pr virus

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
| `-md/--mode` | Search the database is for: its index is built beside the FASTA — `blastx` (`makeblastdb`, default) or a DIAMOND sensitivity (`diamond makedb`). See [Indexing](#indexing). |
| `-p/--threads` | Threads for the `cd-hit` deduplication and `diamond makedb` (default 1). The download itself is network-bound and unaffected. |
| `--refseq/--all-sequences` | Restrict to RefSeq (default) or fetch everything. |
| `--exclude-taxon` | Leave a branch of the taxonomy out of the download entirely; repeatable. Defaults to SARS-CoV-2 for `virus`. See [Excluding a virus from the download](#excluding-a-virus-from-the-download). |
| `--include-unranked/--skip-unranked` | Also download the subtrees with no taxon at `--split-level`, instead of only reporting them. Their lineage still distinguishes them. default = skip. Not on `host`. |
| `--split-level` | Split the download into one request per `family` (default) or `genus`, instead of one request for the whole taxon. `none` restores the single request. Subtrees with no taxon at that rank are skipped and logged. Not on `host`. |
| `--cluster/--no-cluster` | Collapse 100%-identical / 100%-coverage duplicate proteins with `cd-hit` before writing the database (on by default). |
| `--released-before` | Only include data released on or before this date (`YYYY-MM-DD`), so a build can be reproduced later. |
| `--keep-download` | Keep the downloaded zip and the extracted `*_ncbi/` directory (deleted by default). |
| `--attempts` | How many times to try the download before giving up (default 3). |
| `--stall-timeout` | Seconds without any sign of life after which a download is treated as hung, killed and retried (default 180; `0` waits forever). |
| `--debug` | Verbose logging (resolved arguments, the `datasets` command, sequence tallies). |

`-tx` defaults to `10239` (Viruses) for `virus` and `2` (Bacteria) for
`bacteria`, so `eefinder get-databases virus -od db/` downloads the whole RefSeq
viral protein set. For `host`, `-tx` is **required**.

### Per-subcommand options

`--exclude-uninformative` is accepted by **all three** subcommands (on by
default). `--standardize-proteins` exists only on `virus` and `bacteria`, since
`host` produces the baits FASTA and no metadata CSV to standardise.

## Excluding a virus from the download

SARS-CoV-2 is **9,210,670 of the 15,027,633 viral records in GenBank — 61%** of
everything under taxon `10239`, and 99.2% of *Coronaviridae* (numbers from October 11th,
2026).

`--exclude-taxon` leaves it out. It is **not** a post-download filter: the
branch is never requested from NCBI, so nothing about it is transferred.

```bash
# the default: SARS-CoV-2 is left out
eefinder get-databases virus -od db/ --all-sequences

# leave out something else as well (repeatable, tax id or scientific name)
eefinder get-databases virus -od db/ --exclude-taxon 2697049 \
                                     --exclude-taxon "Influenza A virus"

# opt out and download the branch after all
eefinder get-databases virus -od db/ --exclude-taxon none
```

The option is available on `bacteria` and `host` too, where it has no default.

### How it works

The `datasets` CLI has no exclusion flag — its `virus genome taxon` subcommand
can only be told what to *include*. So EEfinder inverts the request. It walks
the lineage from the taxon being downloaded to the taxon being excluded and, at
each step, keeps every child except the one that leads to the exclusion:

```
Viruses (10239)          18 children -> keep 17
  Riboviria                7 children -> keep 6
    Orthornavirae          9 children -> keep 8
      ...
        Sarbecovirus       2 children -> keep 1
          B. pandemicum    2 children -> keep 1   <- SARS-CoV-2 dropped here
```

The union of those siblings is exactly "all viruses except SARS-CoV-2": **63
taxa**, which fit in a single `datasets --inputfile` call (the CLI accepts 100).
Longer lists are split across several downloads and extracted one package per
subdirectory. `--debug` prints the per-rank breakdown, and the run log records
the expansion under `arguments.exclude_taxa`.

The lineage is resolved live rather than hard-coded, which matters: SARS-CoV-2's
parent taxon was recently renamed to *Betacoronavirus pandemicum* by the ICTV.

### The one caveat

Records attached **directly** to a rank on that path — `"Betacoronavirus sp."`
and similar — belong to no child, so this method does not reach them. Measured
against NCBI's own totals that is **273 records out of 5.8 million (0.005%)**,
all of them unclassified entries whose genus and family EEfinder could not have
resolved anyway.

## Progress and interrupted downloads

A whole NCBI download takes a long time, so `get-databases` reports what it
is doing. The `datasets` CLI draws its own transfer progress
(`Downloading: virus.zip  65.5kB 125kB/s`, then `Validating package files`) and
EEfinder passes it straight through instead of swallowing it; the steps EEfinder
performs itself — extracting the archive, merging the `protein.faa` files and
building the metadata table — get their own bars:

```
Downloading: db/virus.zip   1.2GB 8.4MB/s
Validating package files [====================================] 100% 5/5
Extracting         [####################################]  312/312
Merging FASTAs     [####################################]  312/312
Building metadata  [###########-------------------------]  1.2M/4.1M
```

Progress is shown **only on a terminal**: piped into a file or a workflow engine
the redraw sequences would be noise, so the output falls back to being captured
silently. `EEFINDER_NO_PROGRESS=1` turns it off on a terminal too.

### Splitting the download (`--split-level`)

A single request for a broad taxon is the transfer that fails and hangs: asking
for every virus at once moves the whole of RefSeq viral in one package.
`--split-level` asks for **one family at a time** instead (the default), and the
packages are concatenated before the database is built — clustering,
standardisation and the metadata CSV all see the merged set, so the result is
the same as a single download.

```bash
# the default: one request per viral family
eefinder get-databases virus -od db/

# one request per genus instead (more, smaller requests)
eefinder get-databases virus -od db/ --split-level genus

# back to a single request
eefinder get-databases virus -od db/ --split-level none
```

By default only existent taxa (family or genus) are reported. `--include-unranked`
downloads them as well: their full lineage still tells them apart, so including
them does not cluster unclassified viruses together.

```bash
eefinder get-databases virus -od db/ --include-unranked
```

`--split-level` composes with `--exclude-taxon`: the split is planned first, and
an excluded branch is pruned out of whichever split taxon contains it.

### Retries and stalled transfers

Large NCBI downloads fail intermittently, and `get-databases` handles the usual
failure modes on its own: an error the tool reports, a corrupt archive (verified
after every download), an HTTP/2 stream reset (retried over HTTP/1.1, which
completes the same transfer), and a silent hang where the connection stays open
and nothing more arrives.

A transfer counts as **stalled** only when neither new output from `datasets`
nor growth of the archive is seen for `--stall-timeout` seconds (default 180),
so a slow-but-working transfer is never mistaken for a hung one. The attempt is
then killed, the partial archive discarded, and the download retried — up to
`--attempts` times (default 3) with a growing backoff.

A long wait is announced rather than silent, so a slow transfer can be told
apart from a dead one:

```
INFO:eefinder:Still waiting on datasets: no output and no new data for 40s
(giving up and retrying at 180s)
```

Failures that are **permanent** — a misspelled taxon, an unknown flag — are not
retried, since the next attempt would produce the same message:

```
$ eefinder get-databases virus -tx NotARealTaxon -od db/
Failed to download databases: datasets download failed after 1 attempt(s)
(exit 1): Error: The taxonomy name 'NotARealTaxon' is not a recognized virus taxon.
```

## Outputs

```
db/
├── virus.fa              the protein database (screening -db)
├── virus.csv             its metadata (screening -mt)
├── virus.fa.p*           its BLAST index (virus.fa.dmnd with a DIAMOND --mode)
├── virus.tracking.tsv    the treatment of every downloaded accession
├── virus.clstr           the cd-hit clusters: members and representatives
└── virus.log             the JSON run summary
```

The downloaded zip and the `*_ncbi/` directory it is extracted into are
**deleted** once they have been read: for a whole-RefSeq download they are
several GB duplicating what is already in the outputs. `--keep-download` keeps
them.

A download split by `--exclude-taxon` writes one package per batch
(`virus.part1.zip`, …) with a `virus.partN.taxa.txt` listing the taxa it asked
for, extracted into `virus_ncbi/part_N/`. They are cleaned up on the same terms
as the single-package files.

### `tracking.tsv`

One row per sequence **as downloaded** — including the ones that did not make it
into the database, which appear nowhere else:

| Column | Meaning |
|--------|---------|
| `Accession`, `Species` | the record as NCBI delivered it |
| `Protein_downloaded` | the product name in the FASTA header |
| `Protein_final` | the name after `--standardize-proteins` |
| `Name_changed` | `yes` when standardisation renamed it (`RNA-dependent RNA polymerase` → `RdRp`) |
| `Status` | `kept` or `removed` |
| `Reason` | why it was removed: `released_after_cutoff`, `uninformative_product`, `identical_duplicate`, `product_standardized_to_unknown` |
| `Organism_release_date` | earliest release date among the organism's genome records — recorded for every accession, kept or removed, so a dated build can be checked rather than trusted |
| `Cluster`, `Cluster_representative` | the `cd-hit` cluster it landed in and the accession that stands for it — so a dropped duplicate points at the sequence that replaced it |

`Reason` also has a catch-all, `absent_from_final_database`: the table is
reconciled against the finished FASTA, so a record that disappeared without any
step reporting it is still accounted for. It does happen — `cd-hit` silently
discards sequences of up to 10 residues.

```
Accession              Species                   Protein_downloaded              Protein_final  Name_changed  Status   Reason                 Organism_release_date  Cluster  Cluster_representative
AAA17535.1:1-535       Soybean dwarf virus       RNA-dependent RNA polymerase    RdRp           yes           kept                            1990-02-10             2148     AAA17535.1:1-535
AAA02620.1:1-167       Human betaherpesvirus 6A  glycoprotein                    Glycoprotein   yes           kept                            1991-03-15             80639    AAA02620.1:1-167
NP_047212.1:1-1433     Bunyamwera virus          M polyprotein                   M polyprotein  no            kept                            1991-07-30             15       NP_047212.1:1-1433
YP_009117083.1:1-2239  Maprik virus              RNA-dependent RNA polymerase                                 removed  released_after_cutoff  2015-01-31             3012     YP_009117083.1:1-2239
YP_009177301.1:1-210   Suffolk virus             hypothetical protein                                         removed  uninformative_product  2015-11-03             41       YP_009177301.1:1-210
AAA02569.1:1-109       Human coronavirus OC43    12.9 kDa nonstructural protein                               removed  identical_duplicate    1993-03-03             56645    XXN76581.1:1-109
```

`Protein_final` and `Name_changed` are filled only for sequences that reach the
metadata table; a sequence removed before that step leaves them empty.

## Reproducible builds (`--released-before`)

NCBI databases grows continuously, so "the viral database" means nothing without a date.
`--released-before YYYY-MM-DD` restricts a build to data released on or before
that day, and the date is recorded in `{prefix}.log`, so a database can be
rebuilt as it stood then:

```bash
eefinder get-databases virus -od db/ -pr virus_2020 --released-before 2020-12-31
```

**How the cutoff is applied differs by dataset, and the difference matters:**

| Dataset | Applied by | Granularity |
|---------|-----------|-------------|
| `bacteria`, `host` | the NCBI client's own `--released-before` | per assembly (exact) |
| `virus` | EEfinder, after the download | **per organism** |

The `datasets` client implements the cutoff for genome downloads but its **virus**
subcommand offers only `--released-after`, so for viruses EEfinder applies it
itself, using the release dates in the download's own `data_report.jsonl`.

That filtering is per **organism**, not per genome record, because the protein
FASTA carries no link back to the record a protein came from (mapping them by the
report's `proteinCount` was tested and does not hold).

So an organism is included when its **earliest** record predates the cutoff —
i.e. when the virus already existed then — and all of its proteins come with it.
This matters for segmented viruses, whose segments can be released years apart
(La Crosse virus: 2002 and 2023). The alternative, dropping an organism because
one segment came later, would erase viruses that demonstrably existed at the
cutoff.

Every excluded accession appears in `tracking.tsv` with
`Reason = released_after_cutoff` and its `Organism_release_date`, so the cutoff
can be verified line by line. For `bacteria`/`host` the excluded records are never
downloaded, so they cannot appear there — the parameter in `{prefix}.log` is the
record in that case.

## The metadata CSV format

`-mt` is the table that turns a protein accession into a taxonomic assignment.
It must have a header row with these five columns, which EEfinder matches **by
name**:

```text
Accession,Taxonomy,Molecule_type,Protein,Host
YP_009664712.1,r__Riboviria;k__Orthornavirae;p__Negarnaviricota;c__Monjiviricetes;o__Mononegavirales;f__Rhabdoviridae;g__Tibrovirus;s__Bas-Congo tibrovirus,ssRNA(-),N protein,Homo sapiens
YP_009665181.1,r__Riboviria;k__Pararnavirae;p__Artverviricota;c__Revtraviricetes;o__Ortervirales;f__Retroviridae;g__Gammaretrovirus;s__Chick syncytial virus,ssRNA-RT,polymerase,
```

(taxonomy-column)=
### The `Taxonomy` column

The whole lineage lives in one field, one slot per rank, always all eight:

```text
r__{realm};k__{kingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species}
```

A rank NCBI has no taxon at is written `Unk` rather than left empty, so the
field count never varies.

Family- and genus-based options (`--merge_level`, `--target_families`,
`--non_target_families`) read the rank they need out of this string. When that
rank is `Unk`, the **whole lineage** is used as the taxon instead -- so two
unclassified viruses from different orders are never merged into one element,
which is what a single `Unknown` bucket used to do.

The `taxonomy` column of a screening run is the same vocabulary but not
necessarily the same shape: a merged element is assigned the lowest common
ancestor of its hits, and the ranks it does not reach are left off instead of
padded with `Unk` (see [the lineage of a merged
element](output.md#the-lineage-of-a-merged-element)).

The `Accession` values must match the first token of the corresponding FASTA
headers in `-db` (e.g. `>YP_009664712.1 N protein [Bas-Congo tibrovirus]`).
`test_files/virus_subset.csv` is a working miniature example.

**Every protein in `-db` must have a row in `-mt`.** The run stops before the
search if any does not, naming the offending accessions: the taxonomy join is a
left join, so an undescribed protein would otherwise reach the output with an
empty lineage. Extra metadata rows are fine -- a database may be a
subset of the table that describes it.

The header is checked before the analysis starts, and again when the taxonomy is
assembled:

| Header | Behaviour |
|--------|-----------|
| The five columns, in the order above | The run proceeds silently. |
| The five columns, in a different order | A warning is logged and the columns are reordered in memory. The metadata file itself is never rewritten. |
| Extra columns beyond the five | A warning is logged and they are ignored. |
| Any of the five missing | The run stops with an error naming the missing column(s). |

```{note}
The names must be spelled exactly as above — `Molecule_type`, not
`Molecule type`. A wrong column *order* is fixed for you; a missing column
cannot be, which is why it is an error rather than a warning.
```

```{note}
To build a database by hand instead of with `get-databases`, any protein FASTA
works, as long as every accession in it has a row in a CSV with the five
columns above (`Accession,Taxonomy,Molecule_type,Protein,Host`).
```

## How the metadata CSV is built

`get-databases` writes the five columns described [above](#the-metadata-csv-format)
by joining the downloaded `protein.faa` headers with the datasets
`data_report.jsonl`, matched on the organism name:

- `Accession`, `Protein` — from the protein FASTA header (`Protein` after
  [`--standardize-proteins`](#standardize-proteins)).
- `Taxonomy` — the report gives each organism's lineage as a list of tax ids
  **without ranks**, so EEfinder collects every tax id that appears in any
  lineage and resolves their ranks with `datasets summary taxonomy` (once per
  tax id, in batches, rather than once per organism).
  Each lineage is then written as the eight-rank string, `Unk` wherever NCBI
  has no taxon at that rank. `s__` is the organism name as NCBI delivered it,
  which can sit below species rank (a strain or isolate).
- `Host` — the host organism of the organism's first genome record in the
  report.
- `Molecule_type` — looked up by the `f__` rank in a bundled
  [ICTV genome-composition table](https://ictv.global/report/information/virus-properties)
  (`eefinder/data/ictv_genome_composition.tsv`), since the report does not
  carry it. A family missing from the table, an `Unk` family and every
  bacterial family leave it empty.

Two fallbacks keep a build from failing on incomplete taxonomy:

- If the rank lookup fails, a warning is logged and ranks are inferred from
  ICTV name suffixes instead (`-viridae` for the family, a single-word
  `-virus` for the genus). Only `f__`, `g__` and `s__` are filled; every other
  rank is `Unk`.
- A protein whose organism has no entry in the report gets a lineage with only
  `s__` filled, and an empty `Host`.

## Deduplication (`--cluster`)

Each downloaded package (family or genus) is deduplicated **on its own, before the packages are
merged**, so every `cd-hit` run stays small instead of one run over the whole
database. The partitions are processed one at a time; the parallelism comes from
`cd-hit` itself, via `-p/--threads`.

`get-databases` collapses **exact duplicate**
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
downloaded database. Available on every subcommand, on by default.

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

For `virus`, the rules in `eefinder/data/protein_rules.yaml` are additionally
applied. They hold three kinds of adjustment:

- **typos** — whole-word spelling corrections, applied everywhere;
- **rewrites** — regex fixes to the *shape* of a name, such as dropping the
  redundant `protein` of `ORF3 protein` or the `unknown protein` tail of
  `14 kDa protein unknown protein`;
- **proteins** — synonym → canonical name, e.g. every RdRp spelling (including
  compound names like `P2-RdRp` or `CP/RdRp fusion`) → `RdRp`.

A synonym rule is scoped by `Molecule_type` **and** by taxon, so it never
rewrites proteins of another genome group, and a family that names a protein
against its group's convention can override the group-wide rule. `bacteria` gets
the generic cleaning only; `host` is not standardised.

The rules are a first draft meant to be extended, and adding a term never
requires a code change — `eefinder/data/README.md` documents the file, its
scopes and how the competing rules are ranked.

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

## Indexing

The databases must be indexed once for the search engine you intend to use.
`get-databases` builds the index of what it downloads, for the search given by
`-md`/`--mode` (`blastx` by default; any DIAMOND sensitivity builds the same
`.dmnd`), so `screening` finds it beside the FASTA and reuses it. `screening`
builds the BLAST (`makeblastdb`) or DIAMOND (`diamond makedb`) index of `-db`
and `-bt` when it is missing, e.g. for a database built with another `--mode` or
by hand; pass `-id`/`--index_databases` to rebuild it (e.g. after replacing a
FASTA in place).
