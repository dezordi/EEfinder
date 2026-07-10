"""Download EEfinder similarity-search databases via the NCBI datasets CLI.

EEfinder needs, as inputs to the ``screening`` command:

* a **protein database** FASTA (``-db``) plus a **metadata CSV** (``-mt``) with the
  columns ``Accession,Species,Genus,Family,Molecule_type,Protein,Host``; and
* a **host-gene baits** FASTA (``-bt``).

This module automates acquiring them from NCBI RefSeq — the manual procedure
described in the wiki (https://github.com/WallauBioinfo/EEfinder/wiki) — using
the `NCBI datasets <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/>`_ command
line tool (``ncbi-datasets-cli`` in ``env.yml``). Three dataset types are
supported:

* ``virus``    -- RefSeq viral proteins + metadata (the ``-db``/``-mt`` inputs);
* ``bacteria`` -- RefSeq bacterial proteins + metadata (bacteria screening mode);
* ``host``     -- RefSeq host proteins used as ``-bt`` baits (no metadata CSV).

The metadata CSV is rebuilt from two sources bundled in the datasets download:
the ``protein.faa`` headers (``Accession``, ``Protein`` product and ``Species``)
and the ``data_report.jsonl`` taxonomy report (``Genus`` and ``Family`` inferred
from the unranked lineage by ICTV name suffix, plus ``Host``), joined by
organism name. ``Molecule_type`` is not in the datasets report, so it is filled
from the bundled ICTV genome-composition table (``data/``) keyed by family.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import time
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import NamedTuple

import pandas as pd
from Bio import SeqIO

from eefinder import __version__
from eefinder.log import logger
from eefinder.utils import (
    check_outdir,
    DownloadArguments,
    DownloadInfo,
    SequenceCounts,
    StepInfo,
)
from eefinder.normalization import standardize_protein, strip_bracket_tags

#: Dataset types this module can download.
DATASET_CHOICES = ("virus", "bacteria", "host")

#: Dataset types that also produce a metadata CSV.
_METADATA_DATASETS = ("virus", "bacteria")

#: Default NCBI taxon per dataset when ``--taxon`` is omitted. ``virus`` and
#: ``bacteria`` default to their whole-database roots (10239 = Viruses,
#: 2 = Bacteria); ``host`` has no default.
DEFAULT_TAXA = {"virus": "10239", "bacteria": "2"}

#: Product substrings whose proteins carry no taxonomic signal; optionally
#: dropped from a download via ``--exclude-uninformative``.
UNINFORMATIVE_PRODUCTS = ("hypothetical protein", "uncharacterized protein")

#: Columns of the metadata CSV consumed by the ``screening`` command (``-mt``).
METADATA_COLUMNS = [
    "Accession",
    "Species",
    "Genus",
    "Family",
    "Molecule_type",
    "Protein",
    "Host",
]

#: NCBI datasets CLI binary (from ``ncbi-datasets-cli``).
DATASETS_BINARY = "datasets"

#: Bundled ICTV family -> genome-composition table (used for ``Molecule_type``,
#: which the NCBI datasets report does not provide). Sourced from
#: https://ictv.global/virus-properties.
_ICTV_GENOME_TABLE = (
    Path(__file__).resolve().parent / "data" / ("ictv_genome_composition.tsv")
)


def _load_genome_composition() -> "dict[str, str]":
    """Load the ICTV ``family -> genome composition`` map from the data file."""
    table: dict[str, str] = {}
    if _ICTV_GENOME_TABLE.is_file():
        with open(_ICTV_GENOME_TABLE) as handle:
            next(handle, None)  # skip the header row
            for line in handle:
                family, _, genome = line.rstrip("\n").partition("\t")
                if family:
                    table[family] = genome
    return table


#: ICTV genome composition keyed by virus family.
GENOME_COMPOSITION = _load_genome_composition()


def molecule_type_for_family(family: str) -> str:
    """Return the ICTV genome composition (``Molecule_type``) for a family.

    Empty string when the family is unknown (e.g. an unclassified virus or any
    bacterial family, which the ICTV virus table does not cover).
    """
    return GENOME_COMPOSITION.get(family, "")


@dataclass
class ProteinHeader:
    """The three fields parsed out of a protein FASTA header."""

    accession: str
    product: str
    organism: str


@dataclass
class TaxonomyRecord:
    """Per-organism taxonomy pulled from a datasets ``data_report.jsonl``."""

    species: str
    genus: str
    family: str
    mol_type: str
    host: str


def parse_protein_header(header: str) -> ProteinHeader:
    """Parse a protein FASTA header into accession/product/organism.

    Handles both header styles EEfinder databases come in:

    * NCBI *datasets* CDS proteins ---
      ``YP_013613119.1:1-301 P1 [organism=Paris potyvirus 4] [isolate=YLJ]``:
      the organism is the ``[organism=...]`` tag and the product is the text
      before the first bracket group.
    * NCBI Virus RefSeq proteins ---
      ``YP_009664712.1 N protein [Bas-Congo tibrovirus]``: the organism is the
      trailing bare ``[...]`` group.

    In both cases the accession is the first whitespace-delimited token (kept
    whole, including any ``:start-end`` suffix, so it matches the FASTA id
    reported by BLAST).

    Parameters
    ----------
    header : str
        A FASTA header line, with or without a leading ``>``.

    Returns
    -------
    ProteinHeader
        ``organism``/``product`` are empty strings when the header lacks them.
    """
    header = header.lstrip(">").strip()
    accession, _, remainder = header.partition(" ")
    remainder = remainder.strip()

    organism = ""
    # The organism value may contain one level of nested "[...]" (a strain tag),
    # e.g. "[organism=Maize streak virus - A[South Africa]]".
    match = re.search(r"\[organism=((?:[^\[\]]|\[[^\[\]]*\])*)\]", remainder)
    if match:
        # datasets CDS format: product is everything before the first bracket.
        organism = match.group(1).strip()
        first_bracket = remainder.find("[")
        product = remainder[:first_bracket].strip()
    else:
        # RefSeq format: organism is the trailing bare "[...]" group, if any.
        trailing = re.search(r"\[([^\[\]]*)\]\s*$", remainder)
        if trailing:
            organism = trailing.group(1).strip()
            product = remainder[: trailing.start()].strip()
        else:
            product = remainder
    # Defensively drop any leaked "[key=value]" tag (e.g. "[organism=...]") the
    # branch above may have left behind for unusual header layouts.
    product = re.sub(r"\s+", " ", strip_bracket_tags(product)).strip()
    return ProteinHeader(accession=accession, product=product, organism=organism)


def _genus_family_from_lineage(lineage: "list[dict]") -> "tuple[str, str]":
    """Infer ``(genus, family)`` from an unranked datasets ``lineage``.

    The datasets ``data_report.jsonl`` lists the taxonomic lineage as
    ``{"name", "taxId"}`` entries **without** rank information, so ranks are
    inferred from the ICTV virus naming suffixes: families end in ``-viridae``
    and a genus is a single-word name ending in ``-virus``. The most specific
    (last) match of each is used.

    Returns empty strings for names that do not follow the ICTV conventions
    (e.g. unclassified viruses or phages).
    """
    genus = family = ""
    for node in lineage:
        name = node.get("name", "")
        lower = name.lower()
        if lower.endswith("viridae"):
            family = name
        elif " " not in name and lower.endswith("virus"):
            genus = name
    return genus, family


def parse_taxonomy_report(report_path: str) -> dict[str, TaxonomyRecord]:
    """Build an ``organism -> TaxonomyRecord`` map from a ``data_report.jsonl``.

    The datasets virus data report is a JSON-lines file with one record per
    genome. ``Species`` is the ``virus.organismName``, ``Genus``/``Family`` are
    inferred from the (unranked) ``virus.lineage`` via
    :func:`_genus_family_from_lineage`, and ``Host`` is the top-level
    ``host.organismName``. ``Molecule_type`` is left empty: the datasets report
    does not carry it. Fields are read defensively; the first record wins per
    organism.

    Parameters
    ----------
    report_path : str
        Path to the ``data_report.jsonl`` extracted from the datasets download.

    Returns
    -------
    dict[str, TaxonomyRecord]
        Keyed by organism (species) name.
    """
    records: dict[str, TaxonomyRecord] = {}
    with open(report_path) as report:
        for line in report:
            line = line.strip()
            if not line:
                continue
            data = json.loads(line)
            virus = data.get("virus") or data.get("organism") or {}
            organism = virus.get("organismName", "")
            if not organism or organism in records:
                continue

            genus, family = _genus_family_from_lineage(virus.get("lineage", []))
            host = (data.get("host") or {}).get("organismName", "")
            records[organism] = TaxonomyRecord(
                species=organism,
                genus=genus,
                family=family,
                mol_type="",  # not present in the datasets report
                host=host,
            )
    return records


def build_metadata_frame(
    protein_fasta: str,
    taxonomy: dict[str, TaxonomyRecord],
    standardize: bool = False,
    dataset: str = "virus",
) -> pd.DataFrame:
    """Join protein headers with taxonomy into the screening metadata table.

    Parameters
    ----------
    protein_fasta : str
        Concatenated protein FASTA (its headers supply Accession/Protein/Species).
    taxonomy : dict[str, TaxonomyRecord]
        Organism -> taxonomy map from :func:`parse_taxonomy_report`.
    standardize : bool
        When ``True``, rewrite the ``Protein`` column with
        :func:`standardize_protein` (canonical names via the target-specific
        logic, special characters removed, leading letter capitalised). Records
        that standardise to ``"Unknown"`` (bare ``CDS``/``ORF`` directives) are
        **dropped** from the table.
    dataset : str
        The database target (``"virus"``/``"bacteria"``), selecting which
        standardisation logic :func:`standardize_protein` applies.

    Returns
    -------
    pandas.DataFrame
        One row per (kept) protein, with :data:`METADATA_COLUMNS`.
    """
    rows = []
    for record in SeqIO.parse(protein_fasta, "fasta"):
        header = parse_protein_header(record.description)
        tax = taxonomy.get(
            header.organism,
            TaxonomyRecord(header.organism, "", "", "", ""),
        )
        # Molecule_type is absent from the datasets report, so it is looked up
        # from the ICTV genome-composition table by family.
        mol_type = molecule_type_for_family(tax.family)
        if standardize:
            protein = standardize_protein(header.product, mol_type, target=dataset)
            if protein == "Unknown":
                continue  # bare CDS/ORF: drop the record entirely
        else:
            protein = header.product
        rows.append(
            {
                "Accession": header.accession,
                "Species": header.organism,
                "Genus": tax.genus,
                "Family": tax.family,
                "Molecule_type": mol_type,
                "Protein": protein,
                "Host": tax.host,
            }
        )
    return pd.DataFrame(rows, columns=METADATA_COLUMNS)


def build_download_command(
    dataset: str,
    taxon: str,
    zip_path: str,
    refseq: bool = True,
    datasets_bin: str = DATASETS_BINARY,
) -> str:
    """Compose the ``datasets download`` command for a dataset type.

    Parameters
    ----------
    dataset : str
        One of :data:`DATASET_CHOICES`.
    taxon : str
        NCBI taxon name or tax id to download (e.g. ``Flaviviridae``).
    zip_path : str
        Destination path for the downloaded ``ncbi_dataset`` zip.
    refseq : bool
        Restrict to RefSeq sequences (recommended).
    datasets_bin : str
        Path/name of the ``datasets`` executable.

    Returns
    -------
    str
        The shell command (one token per space) to run.
    """
    if dataset not in DATASET_CHOICES:
        raise ValueError(f"Unknown dataset type: {dataset!r}")

    if dataset == "virus":
        command = (
            f"{datasets_bin} download virus genome taxon {shlex.quote(taxon)} "
            f"--include protein "
            f"--filename {zip_path}"
        )
        if refseq:
            command += " --refseq"
    else:  # bacteria / host: genome download, keep only the proteins
        command = (
            f"{datasets_bin} download genome taxon {shlex.quote(taxon)} "
            f"--include protein "
            f"--filename {zip_path}"
        )
        if refseq:
            command += " --assembly-source RefSeq"
    return command


class ConcatCounts(NamedTuple):
    """Result of :func:`concat_protein_faas`: file and sequence tallies."""

    files: int
    total: int
    written: int


def concat_protein_faas(
    extract_dir: str,
    output_fasta: str,
    exclude_products: "tuple[str, ...]" = (),
) -> ConcatCounts:
    """Concatenate every ``protein.faa`` under ``extract_dir`` into one FASTA.

    A datasets download stores proteins under ``ncbi_dataset/data`` (one file
    for viruses, one per assembly for genomes), so all ``protein.faa`` files are
    merged in sorted order.

    Parameters
    ----------
    extract_dir : str
        Directory the datasets zip was extracted into.
    output_fasta : str
        Destination FASTA path.
    exclude_products : tuple[str, ...]
        Lower-cased product substrings whose records are dropped (e.g.
        :data:`UNINFORMATIVE_PRODUCTS`). Empty (the default) keeps everything.

    Returns
    -------
    ConcatCounts
        ``files`` merged, ``total`` sequences seen, and ``written`` sequences
        kept (``total - written`` were dropped as uninformative).

    Raises
    ------
    FileNotFoundError
        If no ``protein.faa`` file is present in the download.
    """
    faas = sorted(Path(extract_dir).rglob("protein.faa"))
    if not faas:
        raise FileNotFoundError(f"no protein.faa found under {extract_dir}")
    exclude = tuple(term.lower() for term in exclude_products)
    total = written = 0
    keep = True
    with open(output_fasta, "w") as out:
        for faa in faas:
            with open(faa) as handle:
                for line in handle:
                    if line.startswith(">"):
                        total += 1
                        if exclude:
                            product = parse_protein_header(line).product.lower()
                            keep = not any(term in product for term in exclude)
                        else:
                            keep = True
                        if keep:
                            written += 1
                    if keep:
                        out.write(line)
    return ConcatCounts(files=len(faas), total=total, written=written)


def find_data_report(extract_dir: str) -> str | None:
    """Return the first ``data_report.jsonl`` under ``extract_dir``, if any."""
    reports = sorted(Path(extract_dir).rglob("data_report.jsonl"))
    return str(reports[0]) if reports else None


def filter_fasta_by_ids(fasta_path: str, keep_ids: "set[str]") -> int:
    """Rewrite ``fasta_path`` in place, keeping only records in ``keep_ids``.

    A record is kept when the first whitespace-delimited token of its header
    (i.e. its accession/id) is in ``keep_ids``. Used to drop from the FASTA the
    same records dropped from the metadata CSV, keeping the two in sync.

    Parameters
    ----------
    fasta_path : str
        FASTA to filter (overwritten via a temporary file).
    keep_ids : set[str]
        Header ids (first token, without ``>``) to retain.

    Returns
    -------
    int
        Number of records dropped.
    """
    tmp_path = f"{fasta_path}.tmp"
    dropped = 0
    keep = True
    with open(fasta_path) as fasta_in, open(tmp_path, "w") as fasta_out:
        for line in fasta_in:
            if line.startswith(">"):
                record_id = line[1:].split(None, 1)[0]
                keep = record_id in keep_ids
                if not keep:
                    dropped += 1
            if keep:
                fasta_out.write(line)
    os.replace(tmp_path, fasta_path)
    return dropped


class GetDatabases:
    """Download and assemble an EEfinder database via the NCBI datasets CLI.

    Runs on instantiation: downloads the datasets zip, extracts it, writes the
    concatenated protein FASTA (``{outdir}/{prefix}.fa``) and, for the ``virus``
    and ``bacteria`` datasets, the metadata CSV (``{outdir}/{prefix}.csv``). A
    JSON run summary with kept/dropped sequence counts is written to
    ``{outdir}/{prefix}.log`` and exposed as ``self.sequence_counts``.

    Parameters
    ----------
    dataset : str
        One of :data:`DATASET_CHOICES`.
    taxon : str
        NCBI taxon name or tax id to download.
    outdir : str
        Output directory (created if missing).
    prefix : str, optional
        Basename for the output files; defaults to ``dataset``.
    refseq : bool
        Restrict to RefSeq sequences (recommended).
    exclude_uninformative : bool
        Drop ``hypothetical protein`` / ``uncharacterized protein`` records
        (:data:`UNINFORMATIVE_PRODUCTS`) from the FASTA and CSV.
    standardize_proteins : bool
        Rewrite the CSV ``Protein`` column to canonical names via the bundled
        viral protein map (:func:`standardize_protein`).
    datasets_bin : str
        Path/name of the ``datasets`` executable.
    """

    def __init__(
        self,
        dataset: str,
        taxon: str,
        outdir: str,
        prefix: "str | None" = None,
        refseq: bool = True,
        exclude_uninformative: bool = False,
        standardize_proteins: bool = False,
        datasets_bin: str = DATASETS_BINARY,
    ) -> None:
        if dataset not in DATASET_CHOICES:
            raise ValueError(f"Unknown dataset type: {dataset!r}")
        self.dataset = dataset
        self.taxon = taxon
        self.outdir = check_outdir(outdir)
        self.prefix = prefix or dataset
        self.refseq = refseq
        self.exclude_uninformative = exclude_uninformative
        self.standardize_proteins = standardize_proteins
        self.datasets_bin = datasets_bin

        self.get_databases()

    def get_databases(self) -> None:
        """Download the dataset and write the FASTA (and CSV where applicable)."""
        run_start = time.time()
        steps: list[StepInfo] = []
        zip_path = f"{self.outdir}/{self.prefix}.zip"
        extract_dir = f"{self.outdir}/{self.prefix}_ncbi"
        fasta_out = f"{self.outdir}/{self.prefix}.fa"
        logger.debug(
            f"Paths: zip={zip_path} extract_dir={extract_dir} fasta={fasta_out}"
        )

        logger.info(f"Downloading {self.dataset} proteins for taxon '{self.taxon}'")
        start = time.time()
        self._download(zip_path)
        steps.append(
            StepInfo.from_times(
                "Download",
                start,
                time.time(),
                f"Downloaded {self.dataset} proteins for taxon '{self.taxon}' "
                f"(refseq={self.refseq}) to {zip_path}.",
            )
        )

        logger.info("Extracting the datasets archive")
        start = time.time()
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(extract_dir)
        logger.debug(f"Extracted archive into {extract_dir}")
        steps.append(
            StepInfo.from_times(
                "Extract",
                start,
                time.time(),
                f"Extracted {zip_path} into {extract_dir}.",
            )
        )

        exclude = UNINFORMATIVE_PRODUCTS if self.exclude_uninformative else ()
        start = time.time()
        counts = concat_protein_faas(extract_dir, fasta_out, exclude_products=exclude)
        excluded_uninformative = counts.total - counts.written
        merged = f"{counts.files} protein.faa file(s) merged"
        if exclude:
            merged += "; dropped hypothetical/uncharacterized proteins"
        logger.info(f"Wrote {fasta_out} ({merged})")
        logger.debug(
            f"Sequences: downloaded={counts.total} "
            f"excluded_uninformative={excluded_uninformative} "
            f"written={counts.written}"
        )
        steps.append(
            StepInfo.from_times(
                "Write protein FASTA",
                start,
                time.time(),
                f"Merged {counts.files} protein.faa file(s): {counts.total} "
                f"sequences, {excluded_uninformative} dropped as uninformative, "
                f"{counts.written} written to {fasta_out}.",
            )
        )

        dropped_standardization = 0
        if self.dataset in _METADATA_DATASETS:
            start = time.time()
            report = find_data_report(extract_dir)
            logger.debug(f"data_report.jsonl: {report}")
            taxonomy = parse_taxonomy_report(report) if report else {}
            logger.debug(f"Parsed taxonomy for {len(taxonomy)} organism(s)")
            if not report:
                logger.warning(
                    "No data_report.jsonl in the download; the metadata CSV will "
                    "have empty Genus/Family/Molecule_type/Host columns."
                )
            frame = build_metadata_frame(
                fasta_out,
                taxonomy,
                standardize=self.standardize_proteins,
                dataset=self.dataset,
            )
            if self.standardize_proteins:
                # Keep the FASTA in sync: drop the records that standardisation
                # removed from the table (bare CDS/ORF and hypothetical proteins).
                dropped_standardization = filter_fasta_by_ids(
                    fasta_out, set(frame["Accession"])
                )
                if dropped_standardization:
                    logger.info(
                        f"Dropped {dropped_standardization} unknown "
                        f"(CDS/ORF/hypothetical) protein(s) from {fasta_out}"
                    )
            csv_out = f"{self.outdir}/{self.prefix}.csv"
            frame.to_csv(csv_out, index=False)
            logger.info(f"Wrote {csv_out} ({len(frame)} records)")
            steps.append(
                StepInfo.from_times(
                    "Write metadata CSV",
                    start,
                    time.time(),
                    f"Wrote {len(frame)} records to {csv_out} "
                    f"(standardize={self.standardize_proteins}); dropped "
                    f"{dropped_standardization} sequence(s) from the FASTA.",
                )
            )

        kept = counts.written - dropped_standardization
        self.sequence_counts = SequenceCounts(
            downloaded=counts.total,
            excluded_uninformative=excluded_uninformative,
            dropped_standardization=dropped_standardization,
            kept=kept,
        )
        logger.info(
            f"Sequences: {kept} kept, "
            f"{excluded_uninformative + dropped_standardization} dropped "
            f"(of {counts.total} downloaded)"
        )
        self._write_log(run_start, time.time(), steps)

    def _write_log(
        self, start_time: float, end_time: float, steps: list[StepInfo]
    ) -> None:
        """Write the ``{prefix}.log`` JSON run summary (kept/dropped counts)."""
        info = DownloadInfo.from_run(
            eefinder_version=__version__,
            arguments=DownloadArguments(
                dataset=self.dataset,
                taxon=self.taxon,
                outdir=self.outdir,
                prefix=self.prefix,
                refseq=self.refseq,
                exclude_uninformative=self.exclude_uninformative,
                standardize_proteins=self.standardize_proteins,
            ),
            sequence_counts=self.sequence_counts,
            start_time=start_time,
            end_time=end_time,
            steps_information=steps,
        )
        log_path = f"{self.outdir}/{self.prefix}.log"
        logger.debug(f"Writing download summary to {log_path}")
        with open(log_path, "w") as json_out:
            json.dump(asdict(info), json_out, indent=4)
        logger.info(f"Wrote {log_path}")

    def _download(self, zip_path: str) -> None:
        """Run ``datasets download``, raising on a non-zero exit."""
        command = build_download_command(
            self.dataset,
            self.taxon,
            zip_path,
            refseq=self.refseq,
            datasets_bin=self.datasets_bin,
        )
        logger.debug(f"datasets command: {command}")
        result = subprocess.run(shlex.split(command), capture_output=True, text=True)
        if result.returncode != 0:
            message = result.stderr.strip() or result.stdout.strip()
            if "ACELLULAR_ROOT" in message or "RankType" in message:
                message += (
                    " -- this is an outdated NCBI datasets CLI: it predates the "
                    "'acellular root' taxonomy rank added above Viruses. Upgrade "
                    "to datasets >= 18.1 (see env.yml)."
                )
            raise RuntimeError(
                f"datasets download failed (exit {result.returncode}): {message}"
            )
