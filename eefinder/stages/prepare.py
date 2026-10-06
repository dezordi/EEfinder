"""Input preparation: validate the reference inputs and build the indexes.

Runs before any search, so a malformed metadata table or an unindexed database
fails up front.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from eefinder.log import logger
from eefinder.make_database import MakeDB
from eefinder.stages.base import Stage, StageOutputs
from eefinder.utils import check_metadata_file

_REPORTED_MISSING = 5


@dataclass
class PrepareInputsOutputs(StageOutputs):
    """Paths validated by :class:`PrepareInputs`.

    A search index lives beside the FASTA it indexes, so the stage echoes back
    the inputs it checked rather than naming new files.
    """

    database: str
    dbmetadata: str
    hostgenesbaits: str
    indexed: bool = False


def _check_protein_fasta(path: str, role: str) -> None:
    """Fail unless ``path`` looks like a non-empty FASTA.

    Parameters
    ----------
    path : str
        File to check.
    role : str
        How to refer to it in an error message (e.g. ``"-db database"``).
    """
    file = Path(path)
    if not file.is_file():
        raise FileNotFoundError(f"the {role} does not exist: {path}")
    if file.stat().st_size == 0:
        raise ValueError(f"the {role} is empty: {path}")
    with open(path) as handle:
        for line in handle:
            if line.strip():
                if not line.startswith(">"):
                    raise ValueError(
                        f"the {role} does not look like a FASTA -- its first "
                        f"content line does not start with '>': {path}"
                    )
                return
    raise ValueError(f"the {role} has no sequence records: {path}")


def _accessions(fasta: str) -> "list[str]":
    """Return the first token of every FASTA header, in file order."""
    accessions = []
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].split()
                accessions.append(fields[0] if fields else "")
    return accessions


def check_database_pairing(database: str, dbmetadata: str) -> int:
    """Fail unless every protein in ``-db`` has a row in ``-mt``.

    The taxonomy join is a left join, so a protein with no metadata row does
    not fail the join -- it reaches the output with empty ``Family``/``Genus``/
    ``Species``. Every accession is therefore checked up front. Extra metadata
    rows are allowed: a database may be a subset of the table that describes it.

    Parameters
    ----------
    database : str
        Protein database FASTA (``-db``).
    dbmetadata : str
        Metadata CSV (``-mt``).

    Returns
    -------
    int
        How many accessions were checked.

    Raises
    ------
    ValueError
        If any accession in ``database`` is absent from ``dbmetadata``.
    """
    accessions = _accessions(database)
    known = set(pd.read_csv(dbmetadata, usecols=["Accession"])["Accession"].astype(str))
    missing = [accession for accession in accessions if accession not in known]
    if missing:
        quoted = ", ".join(missing[:_REPORTED_MISSING])
        if len(missing) > _REPORTED_MISSING:
            quoted += f", ... ({len(missing) - _REPORTED_MISSING} more)"
        raise ValueError(
            f"{len(missing)} of the {len(accessions)} protein(s) in {database} "
            f"have no row in {dbmetadata}: {quoted}. Every protein needs one, "
            "or it reaches the output with no taxonomy. Check that the two "
            "files are the matching pair -- 'eefinder get-databases' writes "
            "them together."
        )
    logger.debug(
        f"all {len(accessions)} accession(s) in {database} are described in "
        f"{dbmetadata}"
    )
    return len(accessions)


class PrepareInputs(Stage):
    """Validate the reference inputs and build the search indexes.

    Checks the ``-mt`` header, that ``-db`` and ``-bt`` are usable protein
    FASTAs, and that the database and its metadata describe the same protein
    set; then builds the BLAST or DIAMOND index for both FASTAs when asked.

    Parameters
    ----------
    database : str
        Protein database FASTA (``-db``).
    dbmetadata : str
        Protein metadata CSV (``-mt``).
    hostgenesbaits : str
        Host-gene bait FASTA (``-bt``).
    outdir, prefix : str
        Run output directory and prefix.
    mode : str
        ``"blastx"`` or a DIAMOND sensitivity, selecting which index to build.
    threads : int
        Threads for the index build.
    index_databases : bool
        Build the indexes. When ``False`` only the checks run.

    Example
    -------
    >>> PrepareInputs(                                        # doctest: +SKIP
    ...     database="db/virus.fa",
    ...     dbmetadata="db/virus.csv",
    ...     hostgenesbaits="db/host.fa",
    ...     outdir="out",
    ...     prefix="run",
    ...     index_databases=True,
    ... ).run().indexed
    True
    """

    name = "prepare"
    title = "Prepare inputs"
    stage_id = "01"

    def __init__(
        self,
        database: str,
        dbmetadata: str,
        hostgenesbaits: str,
        outdir: str,
        prefix: str,
        mode: str = "blastx",
        threads: int = 1,
        index_databases: bool = False,
    ) -> None:
        super().__init__(outdir, prefix)
        self.database = database
        self.dbmetadata = dbmetadata
        self.hostgenesbaits = hostgenesbaits
        self.mode = mode
        self.threads = threads
        self.index_databases = index_databases

    def _execute(self):
        logger.info("Checking metadata file")
        check_metadata_file(self.dbmetadata)

        logger.info("Checking database structure")
        _check_protein_fasta(self.database, "-db database")
        _check_protein_fasta(self.hostgenesbaits, "-bt host-gene baits")
        check_database_pairing(self.database, self.dbmetadata)

        if self.index_databases:
            logger.info("Indexing databases")
            logger.debug(f"MakeDB ({self.mode}) protein DB: {self.database}")
            MakeDB(self.mode, self.database, "prot", self.threads)
            logger.debug(f"MakeDB ({self.mode}) baits DB: {self.hostgenesbaits}")
            MakeDB(self.mode, self.hostgenesbaits, "prot", self.threads)
            message = (
                f"Validated {self.dbmetadata}, {self.database} and "
                f"{self.hostgenesbaits}; built the {self.mode} indexes."
            )
        else:
            logger.warning("index_databases step will not be performed")
            message = (
                f"Validated {self.dbmetadata}, {self.database} and "
                f"{self.hostgenesbaits}; indexing skipped."
            )

        outputs = PrepareInputsOutputs(
            step_info=None,
            database=self.database,
            dbmetadata=self.dbmetadata,
            hostgenesbaits=self.hostgenesbaits,
            indexed=self.index_databases,
        )
        return outputs, message
