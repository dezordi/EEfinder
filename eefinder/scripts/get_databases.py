#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""The ``eefinder get-databases`` command group.

One subcommand per database: ``virus`` and ``bacteria`` produce a protein FASTA
and its metadata CSV, ``host`` produces the baits FASTA.
"""

import click
import shutil
import sys
from eefinder.log import logger, enable_debug
from eefinder.get_databases import (
    GetDatabases,
    DATASETS_BINARY,
    CDHIT_BINARY,
    DEFAULT_TAXA,
)
from eefinder.taxon_exclusion import DEFAULT_VIRUS_EXCLUSIONS, NO_EXCLUSION
from eefinder.progress import DEFAULT_ATTEMPTS, DEFAULT_STALL_TIMEOUT
from eefinder import __version__


def _common_download_options(func):
    """Attach the -od/-pr/--refseq/--debug options shared by every command."""
    func = click.option(
        "--stall-timeout",
        help="Seconds without any sign of life (no output from datasets, no "
        "growth of the archive) after which a download is treated as hung, "
        f"killed and retried. 0 waits forever. default = {DEFAULT_STALL_TIMEOUT}",
        type=float,
        default=DEFAULT_STALL_TIMEOUT,
    )(func)
    func = click.option(
        "--released-before",
        help="Only include data released on or before this date (YYYY-MM-DD), so "
        "a database can be rebuilt as it stood then. For bacteria/host the NCBI "
        "client applies it per assembly; for viruses its subcommand has no such "
        "flag, so EEfinder applies it per organism, from the release dates in the "
        "download's own report.",
        default=None,
    )(func)
    func = click.option(
        "--keep-download/--remove-download",
        help="Keep the downloaded zip and the extracted ncbi_dataset directory. "
        "They are deleted by default, since everything in them is already in the "
        "FASTA, the metadata CSV and the tracking table. default = remove",
        default=False,
    )(func)
    func = click.option(
        "--attempts",
        help="How many times to try the download before giving up. NCBI "
        f"transfers fail and hang intermittently. default = {DEFAULT_ATTEMPTS}",
        type=int,
        default=DEFAULT_ATTEMPTS,
    )(func)
    func = click.option(
        "--debug",
        help="Emit verbose debug logging (download command, extraction, "
        "per-record standardization details). default = off",
        is_flag=True,
    )(func)
    func = click.option(
        "--cluster/--no-cluster",
        help="Collapse 100%-identical / 100%-coverage duplicate proteins with "
        "cd-hit before writing the database. default = cluster",
        default=True,
    )(func)
    func = click.option(
        "--exclude-taxon",
        "exclude_taxa",
        help="Tax id or scientific name to leave out of the download; repeatable. "
        "The branch is never requested from NCBI: the taxon being downloaded is "
        "expanded into the taxa that cover it minus this one. 'virus' defaults "
        f"to {DEFAULT_VIRUS_EXCLUSIONS[0]} (SARS-CoV-2), which is 61% of every "
        f"viral record in GenBank; pass '--exclude-taxon {NO_EXCLUSION}' to "
        "switch that off.",
        multiple=True,
    )(func)
    func = click.option(
        "--refseq/--all-sequences",
        help="Restrict the download to RefSeq (default) or fetch all sequences.",
        default=True,
    )(func)
    func = click.option(
        "-pr",
        "--prefix",
        help="Basename for the output files, default = the dataset type "
        "(e.g. virus.fa/virus.csv).",
        default=None,
    )(func)
    func = click.option(
        "-od",
        "--outdir",
        help="Path and dir to store the downloaded database.",
        required=True,
    )(func)
    return func


def _resolve_exclusions(exclude_taxa, defaults=()):
    """Apply the per-dataset default exclusions and the 'none' escape hatch.

    Passing ``--exclude-taxon none`` clears the defaults, so a build that really
    wants every branch can still ask for one.
    """
    given = tuple(t.strip() for t in (exclude_taxa or ()) if t.strip())
    if not given:
        return tuple(defaults)
    if any(t.lower() == NO_EXCLUSION for t in given):
        kept = tuple(t for t in given if t.lower() != NO_EXCLUSION)
        return kept
    return given


def _run_get_databases(
    dataset,
    taxon,
    outdir,
    prefix,
    refseq,
    exclude_uninformative,
    standardize_proteins,
    cluster=True,
    debug=False,
    attempts=DEFAULT_ATTEMPTS,
    stall_timeout=DEFAULT_STALL_TIMEOUT,
    keep_download=False,
    released_before=None,
    exclude_taxa=(),
):
    """Check for the datasets binary and run :class:`GetDatabases`."""
    if debug:
        enable_debug()
    logger.debug(
        f"get-databases {dataset} arguments: taxon={taxon!r} outdir={outdir!r} "
        f"prefix={prefix!r} refseq={refseq} "
        f"exclude_uninformative={exclude_uninformative} "
        f"standardize_proteins={standardize_proteins} cluster={cluster} "
        f"attempts={attempts} stall_timeout={stall_timeout} "
        f"keep_download={keep_download} released_before={released_before!r} "
        f"exclude_taxa={exclude_taxa!r}"
    )
    if shutil.which(DATASETS_BINARY) is None:
        click.secho(
            f"'{DATASETS_BINARY}' was not found on PATH. Install the NCBI "
            "datasets CLI (conda package 'ncbi-datasets-cli', pinned in env.yml).",
            err=True,
            fg="red",
        )
        sys.exit(1)
    if cluster and shutil.which(CDHIT_BINARY) is None:
        click.secho(
            f"'{CDHIT_BINARY}' was not found on PATH. Install it (conda package "
            "'cd-hit', pinned in env.yml) or pass --no-cluster.",
            err=True,
            fg="red",
        )
        sys.exit(1)
    try:
        GetDatabases(
            dataset=dataset,
            taxon=taxon,
            outdir=outdir,
            prefix=prefix,
            refseq=refseq,
            exclude_uninformative=exclude_uninformative,
            standardize_proteins=standardize_proteins,
            cluster=cluster,
            attempts=attempts,
            stall_timeout=stall_timeout or None,
            keep_download=keep_download,
            released_before=released_before,
            exclude_taxa=exclude_taxa,
        )
    except Exception as err:
        click.secho(f"Failed to download databases: {err}", err=True, fg="red")
        sys.exit(1)


@click.group(name="get-databases")
@click.version_option(__version__)
def get_databases():
    """Download RefSeq protein databases (and metadata) via NCBI datasets.

    Each group has its own subcommand: 'virus' and 'bacteria' produce a protein
    FASTA + metadata CSV (the screening -db/-mt inputs); 'host' produces the -bt
    baits FASTA (no CSV).
    """


@get_databases.command(name="virus")
@_common_download_options
@click.option(
    "-tx",
    "--taxon",
    help="NCBI taxon name or tax id to download (e.g. Flaviviridae, 10239). "
    "default = 10239 (Viruses).",
    default=10239,
    type=str,
)
@click.option(
    "--exclude-uninformative/--keep-uninformative",
    help="Drop 'hypothetical protein' and 'uncharacterized protein' records "
    "from the downloaded database. default = exclude",
    default=True,
)
@click.option(
    "--standardize-proteins/--raw-proteins",
    help="Rewrite the metadata CSV 'Protein' column to canonical names using "
    "the bundled viral protein map (also removes special characters and "
    "capitalises the first letter). default = standardize",
    default=True,
)
def get_databases_virus(
    outdir,
    prefix,
    cluster,
    refseq,
    attempts,
    stall_timeout,
    keep_download,
    released_before,
    debug,
    taxon,
    exclude_uninformative,
    standardize_proteins,
    exclude_taxa,
):
    """Download the RefSeq viral protein DB + metadata CSV (screening -db/-mt)."""
    _run_get_databases(
        dataset="virus",
        taxon=taxon or DEFAULT_TAXA["virus"],
        outdir=outdir,
        prefix=prefix,
        refseq=refseq,
        exclude_uninformative=exclude_uninformative,
        standardize_proteins=standardize_proteins,
        cluster=cluster,
        debug=debug,
        attempts=attempts,
        stall_timeout=stall_timeout,
        keep_download=keep_download,
        released_before=released_before,
        exclude_taxa=_resolve_exclusions(exclude_taxa, DEFAULT_VIRUS_EXCLUSIONS),
    )


@get_databases.command(name="bacteria")
@_common_download_options
@click.option(
    "-tx",
    "--taxon",
    help="NCBI taxon name or tax id to download (e.g. Rickettsiales, 2). "
    "default = 2 (Bacteria).",
    default=2,
    type=str,
)
@click.option(
    "--exclude-uninformative/--keep-uninformative",
    help="Drop 'hypothetical protein' and 'uncharacterized protein' records "
    "from the downloaded database. default = exclude",
    default=True,
)
@click.option(
    "--standardize-proteins/--raw-proteins",
    help="Clean the metadata CSV 'Protein' column (remove NCBI '[key=value]' "
    "tags and special characters, fix common misspellings, capitalise the first "
    "letter, drop bare CDS/ORF records). No bacterial name map is applied yet. "
    "default = standardize",
    default=True,
)
def get_databases_bacteria(
    outdir,
    prefix,
    cluster,
    refseq,
    attempts,
    stall_timeout,
    keep_download,
    released_before,
    debug,
    taxon,
    exclude_uninformative,
    standardize_proteins,
    exclude_taxa,
):
    """Download the RefSeq bacterial protein DB + metadata CSV (screening -db/-mt)."""
    _run_get_databases(
        dataset="bacteria",
        taxon=taxon or DEFAULT_TAXA["bacteria"],
        outdir=outdir,
        prefix=prefix,
        refseq=refseq,
        exclude_uninformative=exclude_uninformative,
        standardize_proteins=standardize_proteins,
        cluster=cluster,
        debug=debug,
        attempts=attempts,
        stall_timeout=stall_timeout,
        keep_download=keep_download,
        released_before=released_before,
        exclude_taxa=_resolve_exclusions(exclude_taxa),
    )


@get_databases.command(name="host")
@_common_download_options
@click.option(
    "-tx",
    "--taxon",
    help="NCBI taxon name or tax id of the host to download (e.g. 'Aedes "
    "aegypti', 7159). Required.",
    required=True,
    type=str,
)
@click.option(
    "--exclude-uninformative/--keep-uninformative",
    help="Drop 'hypothetical protein' and 'uncharacterized protein' records "
    "from the downloaded database. default = exclude",
    default=True,
)
def get_databases_host(
    outdir,
    prefix,
    cluster,
    refseq,
    attempts,
    stall_timeout,
    keep_download,
    released_before,
    debug,
    taxon,
    exclude_uninformative,
    exclude_taxa,
):
    """Download the host protein baits FASTA (screening -bt); no metadata CSV."""
    _run_get_databases(
        dataset="host",
        taxon=taxon,
        outdir=outdir,
        prefix=prefix,
        refseq=refseq,
        exclude_uninformative=exclude_uninformative,
        standardize_proteins=False,
        cluster=cluster,
        debug=debug,
        attempts=attempts,
        stall_timeout=stall_timeout,
        keep_download=keep_download,
        released_before=released_before,
        exclude_taxa=_resolve_exclusions(exclude_taxa),
    )
