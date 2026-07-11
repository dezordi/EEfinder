# EEfinder

> [!WARNING]
> **This is a fork of the [official EEfinder repository](https://github.com/WallauBioinfo/EEfinder).**
> It carries many experimental features still under testing that are not part of
> the upstream release. Use it **at your own risk** — behaviour and outputs may
> change without notice. For research use, prefer the official repository.

EEfinder is a Python CLI and package that automates the identification of
**Endogenous Elements (EEs)** — virus- or bacteria-derived sequences integrated
into eukaryotic genomes — by combining similarity search with genomic-junction
reasoning.

## Documentation

- [Installation](docs/installation.md) — environment and dependencies.
- [Acquiring databases](docs/get-databases.md) — `get-databases` (download, standardisation, 100%/100% clustering).
- [Running the pipeline](docs/screening.md) — `screening` steps and the full option reference.
- [Translation methods](docs/translation-methods.md) — `--translation_method` (`default`/`gv`/`rv`/`gv-rv`).
- [Overlap resolution](docs/overlap.md) — `--overlap` (`keep`/`longest`/`targets`).
- [Custom arguments](docs/custom-arguments.md) — merge length, merge level, range junction, flanks.
- [Outputs](docs/output.md) — result files, taxonomy table and run log.
- [Testing](docs/testing.md) and [Developer guide](docs/developer-guide.md).

## Cite us

If you use EEfinder in your research, please cite:

> Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEfinder: A
> general-purpose tool for identification of bacterial and viral endogenized
> elements in eukaryotic genomes. *Computational and Structural Biotechnology
> Journal.* https://doi.org/10.1016/j.csbj.2024.10.012
