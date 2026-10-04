"""Filesystem and run-bookkeeping helpers for the EEfinder pipeline.

The timing/summary structures emitted at the end of a run are modelled as
dataclasses (:class:`StepInfo`, :class:`RunArguments`, :class:`RunInfo`) so the
run log has a single, typed source of truth. Use :func:`dataclasses.asdict` to
serialise them to JSON.
"""

from __future__ import annotations
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING
import pandas as pd
from eefinder.log import logger

if TYPE_CHECKING:
    from eefinder.versions import DependencyVersion, SystemInfo

#: Timestamp format used throughout the run log.
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"

#: Columns the ``-mt`` metadata CSV must provide, in the order EEfinder uses
#: them downstream.
EXPECTED_METADATA_COLUMNS = [
    "Accession",
    "Species",
    "Genus",
    "Family",
    "Molecule_type",
    "Protein",
    "Host",
]


def check_outdir(outdir: str) -> str:
    """Create the output directory and return its normalised path.

    Parameters
    ----------
    outdir : str
        Desired output directory, with or without a trailing slash.

    Returns
    -------
    str
        The directory path without a trailing slash. The directory (and any
        missing parents) is created if it does not already exist.

    Example
    -------
    >>> check_outdir("results/")  # doctest: +SKIP
    'results'
    """
    outdir = re.sub("/$", "", outdir)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    return outdir


def check_metadata_columns(columns: list) -> list:
    """Validate a metadata header and return the columns EEfinder needs.

    Parameters
    ----------
    columns : list
        Column names of the ``-mt`` metadata table, in file order.

    Returns
    -------
    list
        :data:`EXPECTED_METADATA_COLUMNS`, i.e. the expected columns in the
        order the taxonomy steps read them. Use it to select (and thereby
        reorder) the metadata columns of a ``DataFrame``.

    Raises
    ------
    ValueError
        If any expected column is missing.

    Notes
    -----
    Extra columns and a different column order are warnings, not errors: the
    extras are dropped and the order is fixed in memory. The metadata file on
    disk is never modified.
    """
    missing_columns = [
        column for column in EXPECTED_METADATA_COLUMNS if column not in columns
    ]
    if missing_columns:
        raise ValueError(
            f"the metadata file does not have the column(s): {', '.join(missing_columns)}. "
            f"The metadata file must have the columns: {', '.join(EXPECTED_METADATA_COLUMNS)}."
        )

    extra_columns = [
        column for column in columns if column not in EXPECTED_METADATA_COLUMNS
    ]
    if extra_columns:
        logger.warning(
            f"The metadata file has extra column(s): {', '.join(extra_columns)}. "
            "They will be ignored."
        )

    expected_columns_order = [
        column for column in columns if column in EXPECTED_METADATA_COLUMNS
    ]
    if expected_columns_order != EXPECTED_METADATA_COLUMNS:
        logger.warning(
            f"The metadata file columns are not in the expected order: {', '.join(EXPECTED_METADATA_COLUMNS)}. "
            "They will be reordered in memory, the metadata file is not modified."
        )

    return EXPECTED_METADATA_COLUMNS


def check_metadata_file(metadata_file: str) -> list:
    """Validate the header of the ``-mt`` metadata CSV.

    Only the header row is read, so a malformed metadata file stops the run
    before the similarity search starts rather than at the taxonomy step.

    Parameters
    ----------
    metadata_file : str
        CSV table with the protein taxonomy metadata, parsed from ``-mt``.

    Returns
    -------
    list
        See :func:`check_metadata_columns`.
    """
    header = pd.read_csv(metadata_file, nrows=0).columns.tolist()

    return check_metadata_columns(header)


def _format_timestamp(epoch: float) -> str:
    """Render a POSIX timestamp using :data:`TIME_FORMAT`."""
    return datetime.fromtimestamp(epoch).strftime(TIME_FORMAT)


def _elapsed_minutes(start_time: float, end_time: float) -> str:
    """Return the elapsed wall-clock time in minutes, formatted to 4 decimals."""
    return f"{(end_time - start_time) / 60:.4f}"


@dataclass
class StepInfo:
    """Timing and description record for a single pipeline step."""

    step: str
    start_time: str
    end_time: str
    total_time_minutes: str
    message: str

    @classmethod
    def from_times(
        cls, step: str, start_time: float, end_time: float, message: str
    ) -> "StepInfo":
        """Build a :class:`StepInfo` from POSIX start/end timestamps.

        Parameters
        ----------
        step : str
            Human-readable name of the pipeline step.
        start_time, end_time : float
            POSIX timestamps (e.g. from :func:`time.time`) bounding the step.
        message : str
            Description of what the step produced.

        Returns
        -------
        StepInfo
        """
        return cls(
            step=step,
            start_time=_format_timestamp(start_time),
            end_time=_format_timestamp(end_time),
            total_time_minutes=_elapsed_minutes(start_time, end_time),
            message=message,
        )


@dataclass
class RunArguments:
    """The full set of resolved command-line arguments for a run."""

    genome_file: str
    prefix: str
    outdir: str
    database: str
    dbmetadata: str
    baits: str
    mode: str
    length: int
    flank: int
    limit: int
    range_junction: int
    mask_per: int
    clean_masked: bool
    threads: int
    removetmp: bool
    index_databases: bool
    merge_level: str
    analysis: str
    overlap: str
    target_families: list
    non_target_families: list
    translation_method: str


@dataclass
class DownloadArguments:
    """The resolved arguments of a ``get-databases`` run."""

    dataset: str
    taxon: str
    outdir: str
    prefix: str
    refseq: bool
    exclude_uninformative: bool
    standardize_proteins: bool
    cluster: bool
    #: Download robustness settings (retries and stall detection).
    attempts: int = 1
    stall_timeout: float = 0.0
    keep_download: bool = False
    released_before: str = ""
    #: Taxa left out of the download entirely (never requested from NCBI).
    exclude_taxa: str = ""
    split_level: str = ""


@dataclass
class SequenceCounts:
    """How many sequences a ``get-databases`` run kept vs dropped."""

    downloaded: int
    excluded_uninformative: int
    clustered_identical: int
    dropped_standardization: int
    kept: int


@dataclass
class DownloadInfo:
    """Top-level ``get-databases`` summary serialised to ``{prefix}.log``."""

    eefinder_version: str
    arguments: DownloadArguments
    sequence_counts: SequenceCounts
    start_time: str
    end_time: str
    total_time_minutes: str
    steps_information: list[StepInfo] = field(default_factory=list)
    skipped_taxa: list = field(default_factory=list)
    failed_taxa: list = field(default_factory=list)

    @classmethod
    def from_run(
        cls,
        eefinder_version: str,
        arguments: DownloadArguments,
        sequence_counts: SequenceCounts,
        start_time: float,
        end_time: float,
        steps_information: list[StepInfo],
        skipped_taxa: list = None,
        failed_taxa: list = None,
    ) -> "DownloadInfo":
        """Assemble the download summary from metadata, timestamps and steps.

        Parameters
        ----------
        eefinder_version : str
            Version of EEfinder that produced the download.
        arguments : DownloadArguments
            Resolved arguments used for the download.
        sequence_counts : SequenceCounts
            Kept vs dropped sequence tallies.
        start_time, end_time : float
            POSIX timestamps bounding the whole download.
        steps_information : list[StepInfo]
            One :class:`StepInfo` per download phase.
        skipped_taxa : list, optional
            Taxa left out of a split download, each with its reason.
        failed_taxa : list, optional
            Taxa whose download failed after every attempt.

        Returns
        -------
        DownloadInfo
        """
        return cls(
            eefinder_version=eefinder_version,
            arguments=arguments,
            sequence_counts=sequence_counts,
            start_time=_format_timestamp(start_time),
            end_time=_format_timestamp(end_time),
            total_time_minutes=_elapsed_minutes(start_time, end_time),
            steps_information=steps_information,
            skipped_taxa=list(skipped_taxa or ()),
            failed_taxa=list(failed_taxa or ()),
        )


@dataclass
class RunInfo:
    """Top-level run summary serialised to ``eefinder.log``."""

    eefinder_version: str
    system: "SystemInfo"
    arguments: RunArguments
    dependencies: list["DependencyVersion"]
    start_time: str
    end_time: str
    total_time_minutes: str
    steps_information: list[StepInfo] = field(default_factory=list)

    @classmethod
    def from_run(
        cls,
        eefinder_version: str,
        system: "SystemInfo",
        arguments: RunArguments,
        dependencies: list["DependencyVersion"],
        start_time: float,
        end_time: float,
        steps_information: list[StepInfo],
    ) -> "RunInfo":
        """Assemble the run summary from run metadata, timestamps and step info.

        Parameters
        ----------
        eefinder_version : str
            Version of EEfinder that produced the run.
        system : SystemInfo
            Operating system and host context of the run.
        arguments : RunArguments
            Resolved arguments used for the run.
        dependencies : list[DependencyVersion]
            Detected vs env.yml-pinned versions of the dependencies.
        start_time, end_time : float
            POSIX timestamps bounding the whole run.
        steps_information : list[StepInfo]
            One :class:`StepInfo` per executed pipeline step.

        Returns
        -------
        RunInfo
        """
        return cls(
            eefinder_version=eefinder_version,
            system=system,
            arguments=arguments,
            dependencies=dependencies,
            start_time=_format_timestamp(start_time),
            end_time=_format_timestamp(end_time),
            total_time_minutes=_elapsed_minutes(start_time, end_time),
            steps_information=steps_information,
        )
