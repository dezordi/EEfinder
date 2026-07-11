# EEfinder Documentation

**EEfinder** is a Python CLI and package that automates the identification of
**Endogenous Elements (EEs)** — virus- or bacteria-derived sequences integrated
into eukaryotic genomes — by combining similarity search with genomic-junction
reasoning.

It was published in the *Computational and Structural Biotechnology Journal*
(Dias, Dezordi & Wallau, 2024,
[doi:10.1016/j.csbj.2024.10.012](https://doi.org/10.1016/j.csbj.2024.10.012)).

## Main features

- **Two CLI commands** — `screening` (the EE-finding pipeline) and
  `get-databases` (automated RefSeq database download via the NCBI `datasets`
  CLI).
- **Similarity search + junction reasoning** — filters redundant hits, merges
  fragmented elements of the same taxon, and removes false positives that match
  host genes more strongly than the viral/bacterial reference.
- **Selectable translation methods** — the classic six-frame `blastx`/`diamond
  blastx`, or up-front protein prediction with
  [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) /
  [pyrodigal-rv](https://github.com/LanderDC/pyrodigal-rv) (`gv`/`rv`/`gv-rv`).
- **Overlap resolution** — keep, longest-wins, or target-family strategies for
  elements assigned to different families in the same region.
- **Reproducible run logs** — a JSON `eefinder.log` recording the version,
  resolved arguments, dependency versions (flagging drift from `env.yml`), and
  per-step timing.

## Documentation contents

```{toctree}
:maxdepth: 2

installation
get-databases
screening
translation-methods
overlap
custom-arguments
output
testing
developer-guide
```

## Quick links

- [GitHub repository](https://github.com/WallauBioinfo/EEfinder)
- [Issues / bugs](https://github.com/WallauBioinfo/EEfinder/issues)

## Cite us

If you use EEfinder in your research, please cite:

> Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEfinder: A
> general-purpose tool for identification of bacterial and viral endogenized
> elements in eukaryotic genomes. *Computational and Structural Biotechnology
> Journal*. https://doi.org/10.1016/j.csbj.2024.10.012
