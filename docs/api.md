# Python API

Every stage of the `screening` logic is a class in `eefinder.stages`. They
all have the same shape: build one with its inputs, call `run()`, read the
returned outputs.

```python
from eefinder.stages import DataCleaning

result = DataCleaning(
    genome_file="genome.fa",
    outdir="results",
    prefix="Ae_aeg_Aag2",
    length=1000,
).run()

result.cleaned_genome      # 'results/Ae_aeg_Aag2.cleaned_genome.fa'
result.kept, result.total  # 1, 2
result.step_info           # timing record, as it appears in eefinder.log
```

`run()` returns a dataclass, not a path string, so a caller never has to
reconstruct a file name. `result.files()` gives every produced path as a dict,
and `result.existing()` narrows that to the ones on disk.


## Running the whole pipeline

```python
from eefinder.stages import ScreeningPipeline

result = ScreeningPipeline(
    genome_file="genome.fa",
    database="db/virus.fa",
    dbmetadata="db/virus.csv",
    hostgenesbaits="db/host.fa",
    outdir="results",
    prefix="Ae_aeg_Aag2",
    index_databases=True,
).run()

result.final_tax     # 'results/Ae_aeg_Aag2.EEs.tax.tsv'
result.step_infos    # one StepInfo per stage, in execution order
```

## Chaining stages yourself

Each stage takes the previous one's outputs, so composing them is explicit:

```python
from eefinder.stages import (
    DataCleaning, SequenceAlignment, PutativeElementsFilter,
    TaxonomyAssignment, MergeFragmentedElements,
)

common = dict(outdir="results", prefix="Ae_aeg_Aag2")

cleaned = DataCleaning(genome_file="genome.fa", length=1000, **common).run()
aligned = SequenceAlignment(
    genome=cleaned.cleaned_genome, database="db/virus.fa", threads=8, **common
).run()
filtered = PutativeElementsFilter(
    candidates=aligned.ee_candidates,
    ee_hits_filtered=aligned.ee_hits_filtered,
    hostgenesbaits="db/host.fa",
    threads=8,
    **common,
).run()
assigned = TaxonomyAssignment(
    ee_hits_validated=filtered.ee_hits_validated,
    dbmetadata="db/virus.csv",
    **common,
).run()
merged = MergeFragmentedElements(
    taxonomy_signature=assigned.taxonomy_signature,
    genome=cleaned.cleaned_genome,
    limit=100,
    **common,
).run()
```

## The stages

| Class | Subcommand | Inputs | Outputs |
|-------|-----------|--------|---------|
| `PrepareInputs` | `prepare` | `database`, `dbmetadata`, `hostgenesbaits`, `index_databases` | `database`, `dbmetadata`, `hostgenesbaits`, `indexed` (raises unless every `-db` protein has a `-mt` row) |
| `DataCleaning` | `clean` | `genome_file`, `length` | `cleaned_genome`, `kept`, `total` |
| `SequenceAlignment` | `align` | `genome`, `database`, `mode`, `translation_method`, `range_junction`, `search_method` | `ee_hits`, `ee_hits_filtered`, `ee_candidates_bed`, `ee_candidates` |
| `PutativeElementsFilter` | `filter` | `candidates`, `ee_hits_filtered`, `hostgenesbaits` | `host_hits`, `host_hits_filtered`, `ee_hits_validated` |
| `TaxonomyAssignment` | `taxonomy` | `ee_hits_validated`, `dbmetadata` | `taxonomy_signature` |
| `MergeFragmentedElements` | `merge` | `taxonomy_signature`, `genome`, `limit`, `merge_level` | `merge_bed`, `merge_bed_merged`, `merge_elements_bed`, `ee_elements`, `ee_elements_tax` |
| `PostProcessing` | `postprocess` | `ee_elements`, `ee_elements_tax`, `clean_masked`, `mask_per`, `overlap`, `analysis` | `ee_elements_gff3`, `ee_elements_cleaned*`, `removed_*` |
| `FlanksExtraction` | `flanks` | `ee_elements`, `genome`, `flank` | `genome_lengths`, `flanks_bed`, `flanks_slop_bed`, `flanks` |
| `ScreeningPipeline` | `all` | all of the above | `final_fasta`, `final_tax`, `final_gff3`, `final_flanks`, `step_infos` |

Every stage also takes `outdir` and `prefix`. The output directory is created
if it does not exist, so a stage can be run on its own.

## File names

`ScreeningPaths` is the single source of truth for every path a run uses:

```python
from eefinder.stages import ScreeningPaths

paths = ScreeningPaths(outdir="results", prefix="Ae_aeg_Aag2")
paths.cleaned_genome       # 'results/Ae_aeg_Aag2.cleaned_genome.fa'
paths.taxonomy_signature   # 'results/Ae_aeg_Aag2.taxonomy_signature.csv'
paths.final_tax            # 'results/Ae_aeg_Aag2.EEs.tax.tsv'
paths.intermediates()      # every intermediate, in pipeline order
```

A stage run on its own defaults its inputs to these names, which is why
`eefinder screening align -od results -pr Ae_aeg_Aag2 -db db/virus.fa` needs no
`--genome`.

(output-models)=
## Output models

A model is a schema plus a check. Each stage declares the models its outputs
satisfy, and `run()` enforces them — so a stage states its guarantees rather
than documenting them.

```python
from eefinder.models import HIT_TABLE, ELEMENT_TABLE, ModelError

HIT_TABLE.columns          # ('qseqid', 'sseqid', 'pident', ...)
HIT_TABLE.validate("hits.tsv")      # raises ModelError if malformed
ELEMENT_TABLE.validate("x.tax.tsv")
```

| Model | Shape |
|-------|-------|
| `HIT_TABLE` | BLAST `outfmt 6`, headerless, nucleotide query coordinates |
| `FILTERED_HIT_TABLE` | `outfmt 6` plus `sense`, `bed_name`, `tag` |
| `CANDIDATE_BED` | headerless BED3 |
| `TAXONOMY_SIGNATURE` | per-hit signature CSV |
| `ELEMENT_TABLE` | per-element taxonomy table |
| `CANDIDATE_FASTA`, `ELEMENT_FASTA`, `GENOME_FASTA` | FASTA |

## Search methods

How hits are produced is pluggable. `HIT_TABLE` is the interface: a method
takes a nucleotide query and a protein database and writes that table.

```python
from eefinder.search_methods import SEARCH_METHODS, resolve_search_method

sorted(SEARCH_METHODS)                              # ['blastx', 'diamond', 'predicted']
resolve_search_method("blastx", "default").name     # 'blastx'
resolve_search_method("very-sensitive", "default")  # DiamondSearch
resolve_search_method("blastx", "gv-rv").name       # 'predicted'
```

`mode` and `translation_method` are resolved to one of these.

### Choosing a translation method

`translation_method` selects how the proteins are obtained, and `mode` which
aligner runs. Both are plain stage arguments:

```python
from eefinder.stages import SequenceAlignment

SequenceAlignment(
    genome="results/run.cleaned_genome.fa",
    database="db/virus.fa",
    outdir="results",
    prefix="run",
    translation_method="gv",   # 'default', 'gv', 'rv' or 'gv-rv'
    mode="blastx",             # 'blastx', or a DIAMOND sensitivity
    threads=8,
).run()
```

| `translation_method` | Proteins come from | Aligner |
|---------------------|--------------------|---------|
| `default` | six-frame translation of the query | `blastx` / `diamond blastx` |
| `gv` | pyrodigal-gv prediction | `blastp` / `diamond blastp` |
| `rv` | pyrodigal-rv prediction | `blastp` / `diamond blastp` |
| `gv-rv` | both predictors, deduplicated with cd-hit | `blastp` / `diamond blastp` |

The valid values are `eefinder.translation.TRANSLATION_METHODS`. The prediction
methods need `pyrodigal-gv`/`pyrodigal-rv`, and `gv-rv` also needs `cd-hit`;
a missing binary raises before the search starts.

```{important}
Pass the **same** `translation_method` to `SequenceAlignment` and to
`PutativeElementsFilter`. The two searches of a run are meant to be comparable,
and the CLI guarantees this by threading one value into both — composing the
stages by hand puts that back on the caller.
```

## Dependency versions

```python
from eefinder.versions import write_versions_yml

write_versions_yml("versions.yml", "EEFINDER_ALIGN")
```

```yaml
"EEFINDER_ALIGN":
    eefinder: 2.0.0
    bedtools: 2.27.1
    blast: 2.5.0
    diamond: 2.0.15
```

Every `screening` subcommand takes `--versions-yml PATH` and `--versions-key`
to write the same file from the command line.
