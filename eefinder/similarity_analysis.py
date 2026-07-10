"""Run the translated similarity search (blastx / DIAMOND blastx)."""

from __future__ import annotations

import shlex
import subprocess

from Bio.Blast.Applications import NcbiblastxCommandline

#: E-value cutoff shared by both search backends.
EVALUE_CUTOFF = 0.00001


def runblastx(query_file: str, database_file: str, threads: int) -> None:
    """Run NCBI ``blastx`` writing tabular (outfmt 6) to ``{query_file}.blastx``."""
    cline = NcbiblastxCommandline(
        query=query_file,
        db=database_file,
        out=f"{query_file}.blastx",
        outfmt=6,
        word_size=3,
        evalue=EVALUE_CUTOFF,
        num_threads=threads,
        matrix="BLOSUM45",
        max_intron_length=100,
        soft_masking="true",
    )
    cline()


def rundiamond(query_file: str, database_file: str, threads: int, mode: str) -> None:
    """Run ``diamond blastx`` in the given sensitivity ``mode``."""
    clinedmd = (
        f"diamond blastx "
        f"-p {int(threads)} "
        f"-d {database_file}.dmnd "
        f"-f 6 "
        f"-q {query_file} "
        f"-o {query_file}.blastx "
        f"-e {EVALUE_CUTOFF} "
        f"--matrix BLOSUM45 "
        f"-k 500 "
        f"--max-hsps 0 "
        f"--{mode}"
    )
    subprocess.run(
        shlex.split(clinedmd),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


class SimilaritySearch:
    """Run a BLAST or DIAMOND translated search of ``query`` vs ``database``.

    Runs on instantiation and writes tabular output to ``{query_file}.blastx``.

    Parameters
    ----------
    query_file : str
        Nucleotide FASTA used as the query.
    database_file : str
        Protein database created by :class:`~eefinder.make_database.MakeDB`.
    threads : int
        Number of threads for the search.
    mode : str
        ``"blastx"`` for NCBI BLAST, otherwise a DIAMOND sensitivity mode.
    """

    def __init__(
        self, query_file: str, database_file: str, threads: int, mode: str
    ) -> None:
        self.query_file = query_file
        self.database_file = database_file
        self.threads = threads
        self.mode = mode

        self.similarity_search()

    def similarity_search(self) -> None:
        """Dispatch to the BLAST or DIAMOND backend based on ``mode``."""
        if self.mode == "blastx":
            runblastx(self.query_file, self.database_file, self.threads)
        else:
            rundiamond(self.query_file, self.database_file, self.threads, self.mode)
