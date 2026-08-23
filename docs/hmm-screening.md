# Profile-HMM screening (`-md hmmer`)

Endogenous elements are frequently old, degraded and only remotely similar to
extant viral proteins — exactly the regime where pairwise search loses elements.
`screening -md hmmer` searches **profile HMMs** (HMMER3) instead, which pool the
information of a whole protein family into position-specific scores.

The engine slots into the existing pipeline: once the search emits the usual
tabular hit table, every downstream step (filtering, BED extraction, taxonomy,
merging, mask cleaning, overlap resolution, GFF3, flanks) runs unchanged.

```{note}
`-md hmmer` needs `hmmsearch` on `PATH` (`hmmer` in `env.yml`), a viral profile
database (`get-databases virus-hmm`) and, by default, a host profile database
(`get-databases host-hmm`).
```

## Quick start

```bash
# 1. viral profiles: built from a RefSeq download of one clade, or downloaded
#    from VOGDB (--source vogdb), or both merged into one database
eefinder get-databases virus-hmm -tx Chuviridae -od db_hmm/ -pr virus_hmm -p 8

# 2. host profiles — build with the SAME options (see "Comparability" below)
eefinder get-databases host-hmm -tx "Aedes aegypti" -od db_hmm/ -pr host_hmm -p 8

# 3. screen
eefinder screening -in genome.fa -od out/ -md hmmer \
    -dbh db_hmm/virus_hmm.hmm -mth db_hmm/virus_hmm.csv \
    -bth db_hmm/host_hmm.hmm -p 8
```

## How the search works

```
contigs → six-frame translation → hmmsearch → --domtblout → coordinate traceback
        → {query}.blastx → FilterTable → … (unchanged pipeline)
```

1. **Six-frame translation.** `hmmsearch` needs protein queries. The default
   `--translation_method` becomes `sixframe`: all six frames are translated,
   split at stop codons, and every segment of ≥20 aa is kept. This is the honest
   analogue of `blastx` — an ORF predictor (`gv`/`rv`/`gv-rv`, still selectable)
   would miss pseudogenised elements. Each segment records its contig
   coordinates.
2. **`hmmsearch`** runs once per profile database with `--domtblout`.
3. **Conversion + traceback.** The domain table is rewritten into the 12
   `outfmt 6` columns and the amino-acid coordinates are mapped back to contig
   nucleotides, so the emitted `{query}.blastx` is indistinguishable in shape from
   a `blastx` result.

### Column semantics in HMM mode

| Column | Meaning |
|--------|---------|
| `sseqid` | the **profile id** (joins the `Accession` column of `-mth`) |
| `pident` | `acc × 100`, the domain's mean **posterior probability** — *not* percent identity, which a profile hit does not have. It also feeds `Average_pident` in the taxonomy table |
| `length` | aligned length in amino acids (as in `blastx`) |
| `sstart`/`send` | profile coordinates (`hmm from`/`hmm to`) |
| `evalue` | the independent per-domain E-value (`i-Evalue`) |
| `bitscore` | the per-domain score |

```{warning}
`Average_pident` in a `-md hmmer` run is a posterior-accuracy value (typically
80–100), not a sequence identity. Do not compare it with the identities of a
`blastx` run.
```

### Reproducible E-values

`hmmsearch` E-values scale with the number of searched sequences, which varies
with genome size — the same element would get different E-values in two genomes.
EEfinder therefore fixes `-Z`/`--domZ` (`--hmm_z`, default `1000000`) and records
it in `eefinder.log`.

## Discarding host genes: the symmetric comparison

The classic rule (“keep the element when the viral hit out-scores the host hit”)
relies on both numbers coming from the same program. **A profile-HMM bit score and
a BLAST bit score are not on a common scale**: a deep curated profile out-scores a
pairwise alignment as a matter of course, so mixing them would systematically
keep candidates that are really host genes. EEfinder never does that.

Instead, `-md hmmer` compares **like with like**: the candidate regions are
translated once and searched against *both* profile databases, and the decision
uses **bit score per aligned residue**:

```
keep  if  density(viral) > density(host) × (1 + --hmm_score_margin)
```

Why this form:

- **bit scores, not E-values** — a bit score does not depend on database size or
  `-Z`; an E-value does, and the two databases differ in size by design;
- **density, not raw score** — raw scores grow with alignment length, so a long
  weak host match would beat a short strong viral one;
- **same builder on both sides** — see below.

| `--hmm_host_filter` | What it does |
|---------------------|--------------|
| `hmm` (default) | the symmetric comparison above; needs `-bth` |
| `blastp` | falls back to the `-bt` baits, dropping an element only when a host hit is **significant** (E ≤ `--hmm_evalue`) **and** covers ≥ `--host_min_coverage` of it — a significance/coverage test, never a cross-engine score comparison |
| `none` | no host check at all; warns, and host genes homologous to viral proteins will be reported as elements |

Running `-md hmmer` without `-bth` is an **error**, not a silent downgrade to a
weaker filter.

### Comparability of the two databases

The comparison only means something when both profile sets were produced the same
way. Every database records its builder settings in its `{prefix}.log`, and
`screening` refuses to run the symmetric filter when they differ, naming the
setting:

```
virus_hmm.hmm and host_hmm.hmm were built with different settings
(cluster_identity), so their profile scores are not comparable.
Rebuild one of them, or pass --allow_builder_mismatch.
```

### Host-calibrated thresholds (`calibrate-hmm` + `--hmm_use_ga`)

An orthogonal, cheap refinement that composes with the comparison above. It asks,
once per database and host, **what score each viral profile achieves on the
host's own proteins**, and stores that score (plus a margin) as the profile's
`GA` gathering threshold:

```bash
eefinder calibrate-hmm -dbh db_hmm/virus_hmm.hmm -bt db_hmm/host.faa -p 8
eefinder screening ... -md hmmer --hmm_use_ga ...
```

A hit then only counts when it beats what that profile can achieve on host genes
— both sides of *that* comparison being the same profile scoring real sequences.
It catches the "this profile has cellular homologs" case before the per-region
comparison happens: in practice the raised thresholds land on exactly the
retroelement-related profiles (reverse transcriptase, gag-pol, ORF B).

| Option | Default | Meaning |
|--------|---------|---------|
| `--margin` | `0.1` | relative margin added to each observed host score |
| `--quantile` | `1.0` | which point of the host-score distribution to take; `0.99` ignores a single outlier, useful when the proteome may itself contain an endogenised element |
| `--hmm_evalue` / `--hmm_z` | `1e-5` / `1e6` | the significance the **floor** corresponds to |

Two things worth knowing:

- profiles that never score on the host do **not** get a zero threshold — they get
  the bit score equivalent to `--hmm_evalue` at `--hmm_z`, because `--cut_ga`
  *replaces* the E-value cutoff rather than adding to it, and a zero floor would
  make a calibrated database more permissive than an uncalibrated one;
- calibration is **host-specific** (a database calibrated for *Aedes* is not valid
  for *Drosophila*), and `--hmm_use_ga` applies to the `-dbh` databases only —
  never to the host database, which is the thing being compared against, not
  calibrated.

Calibrating against the curated `-bt` bait set is safer than against a whole
proteome, which may already contain endogenised elements that would raise
thresholds and cost true positives.

### Residual bias, and what to do about it

Viral profiles come from whole RefSeq families across many species; a single host
proteome yields many shallow clusters, and a deeper profile scores a true homolog
higher. Untreated, this biases the comparison toward *keeping* candidates. Two
mitigations:

- build the host database from **several related proteomes** (repeat `-tx`), which
  deepens host families the same way the viral side is deep;
- require the viral side to win by a margin (`--hmm_score_margin 0.1`–`0.2`);
- calibrate the viral database (above), which handles the host-homologous profiles
  independently of the per-region comparison.

## Hybrid mode: profiles discover, BLAST names (`--taxonomy_refine blast`)

A profile is cross-taxon by construction, so its taxonomy is only as precise as
the group it was built from — for a prebuilt set, frequently a whole family or
worse. Hybrid mode keeps the profile search as the **discovery** step and takes
the **taxonomy** from a BLAST search of the surviving candidate regions against
the sequence database:

```bash
eefinder screening -in genome.fa -od out/ -md hmmer \
    -dbh db_hmm/virus_hmm.hmm -mth db_hmm/virus_hmm.csv \
    -bth db_hmm/host_hmm.hmm \
    --taxonomy_refine blast -db db/virus.fa -mt db/virus.csv
```

For every element with a sequence-level hit, the reported protein becomes that
hit's accession and `Average_pident` becomes its **real** percent identity.
Elements only the profiles could place keep their profile id and posterior
accuracy — that population is the reason for using profiles in the first place,
and the run log reports how many of each there were. The pre-refinement
assignments are kept as `*.concat.nr.unrefined` under `tmp_files/`.

The search is cheap: it runs over the candidate regions only, not the genome.

## Which source names an element (`--source_priority`)

When several `-dbh` databases hit the same region, `FilterTable` keeps the
highest-scoring hit — and profile depth drives bit score, so a deep prebuilt
profile usually out-scores a shallow taxon-scoped one. The label then comes from
the source with the coarser assignment, exactly the one the other source exists
to improve.

```bash
eefinder screening ... --source_priority NCBIREFSEQ,VOGDB --source_priority_margin 0.05
```

When a preferred source has a hit within the margin of the best score, its hit
names the element. **Detection is unaffected** — only the label changes. With
`--taxonomy_refine blast` the question is moot, since the taxonomy comes from the
sequence database either way.

## Options

| Option | Default | Purpose |
|--------|---------|---------|
| `-dbh/--hmmdatabase` | — | viral profile database; **repeatable** (one `-mth` each) |
| `-mth/--hmmmetadata` | — | its metadata CSV, in the same order as `-dbh` |
| `-bth/--hostgeneshmm` | — | host profile database; required for the default filter |
| `--hmm_host_filter` | `hmm` | `hmm` / `blastp` / `none` |
| `--hmm_score_margin` | `0.0` | margin the viral side must win by |
| `--host_min_coverage` | `0.5` | coverage term of the `blastp` fallback |
| `--hmm_evalue` | `1e-5` | `-E`/`--domE` |
| `--hmm_z` | `1000000` | fixed `-Z`/`--domZ` |
| `--hmm_min_coverage` | `0.0` | minimum fraction of the profile a hit must cover |
| `--hmm_use_ga` | off | use each profile's `GA` threshold (`--cut_ga`); requires every profile to carry one |
| `--hmm_sensitive` | off | pass `--max` to `hmmsearch` (much more sensitive, much slower) |
| `--allow_builder_mismatch` | off | run the symmetric filter on databases built differently |
| `--taxonomy_refine` | `none` | `blast` takes each element's taxonomy from a sequence-level search (needs `-db`/`-mt`) |
| `--source_priority` | — | comma-separated source namespaces preferred when scores are close |
| `--source_priority_margin` | `0.05` | relative margin within which `--source_priority` decides |
| `--translation_method` | `sixframe` | `sixframe` or a prediction method (`gv`/`rv`/`gv-rv`) |

## Searching several databases

`-dbh`/`-mth` are repeatable, so profile sets from different sources can be
searched together; each keeps its own thresholds and calibration state, and the
hit tables are concatenated before filtering:

```bash
eefinder screening -in genome.fa -od out/ -md hmmer \
    -dbh db_hmm/chuviridae.hmm -mth db_hmm/chuviridae.csv \
    -dbh db_hmm/flaviviridae.hmm -mth db_hmm/flaviviridae.csv \
    -bth db_hmm/host_hmm.hmm
```

Profile ids are namespaced per source (`NCBIREFSEQ__…`), so the merged metadata
never collides; `screening` errors out if two databases share an `Accession`.
Redundant hits are not a problem: `FilterTable` already keeps one hit per
`(contig, --range_junction window, strand)`, so the same region hit by two
profiles yields **one** element.

## Runtime notes

- `hmmsearch` cost is roughly linear in the number of profiles, so two databases
  cost about the sum of both.
- The host search runs only over the extracted candidate regions, so it is cheap;
  the one-off cost is building the host database.
- `--hmm_sensitive` (`--max`) is roughly an order of magnitude slower.

## Choosing between the engines

| Configuration | Command | When |
|---------------|---------|------|
| classic | `-md blastx` | baseline; reproduces the published results |
| profile HMM | `-md hmmer` | remote homology / degraded elements, with the symmetric host filter |
| hybrid | `-md hmmer --taxonomy_refine blast` | profile sensitivity with sequence-level taxonomy — the recommended configuration when a sequence database is available |
| DIAMOND | `-md very-sensitive` | large genomes where BLAST is too slow |
