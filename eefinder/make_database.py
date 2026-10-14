"""Build BLAST or DIAMOND protein databases from a FASTA file."""

from __future__ import annotations
import os
import shlex
import subprocess
from Bio.Blast.Applications import NcbimakeblastdbCommandline

BLAST_MODES = ("blastx", "tblastn")


def makeblastdb(data: str, db_type: str) -> None:
    """Create an NCBI BLAST database from ``data`` (``nucl`` or ``prot``)."""
    clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file=data)
    clinedb()


def has_index(mode: str, data: str) -> bool:
    """Whether the ``mode`` search index for ``data`` already exists.

    BLAST writes ``.pin`` (or a ``.pal`` alias for multi-volume databases);
    DIAMOND writes ``.dmnd``.
    """
    suffixes = (".pin", ".pal") if mode in BLAST_MODES else (".dmnd",)
    return any(os.path.exists(data + suffix) for suffix in suffixes)


def makediamonddb(data: str, threads: int) -> None:
    """Create a DIAMOND database (``{data}.dmnd``) with the BLOSUM45 matrix."""
    clinedb = (
        f"diamond makedb "
        f"--db {data} "
        f"--in {data} "
        f"--threads {int(threads)} "
        f"--matrix BLOSUM45"
    )
    subprocess.run(
        shlex.split(clinedb),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )


class MakeDB:
    """Create the similarity-search database for the chosen ``mode``.

    BLAST modes use ``makeblastdb``; DIAMOND modes use ``diamond makedb``.
    Runs on instantiation.

    Parameters
    ----------
    mode : str
        Search mode selected with ``--mode``.
    data : str
        Path to the database FASTA file, parsed from ``--database``.
    db_type : str
        ``"nucl"`` or ``"prot"`` (only used by BLAST).
    threads : int
        Number of threads (only used by DIAMOND).
    """

    def __init__(self, mode: str, data: str, db_type: str, threads: int) -> None:
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads

        self.make_db()

    def make_db(self) -> None:
        """Dispatch to the BLAST or DIAMOND database builder."""
        if self.mode in BLAST_MODES:
            makeblastdb(self.data, self.db_type)
        else:
            makediamonddb(self.data, self.threads)
