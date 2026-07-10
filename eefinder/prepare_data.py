"""Input preparation: tag genome FASTA headers with a run-specific prefix."""

from __future__ import annotations

from pathlib import Path


class InsertPrefix:
    """Prefix every FASTA header of a genome so element IDs are traceable.

    Each header ``>contig`` becomes ``>{prefix}/contig``. The class runs on
    instantiation and writes ``{outdir}/{prefix}.rn``.

    Parameters
    ----------
    input_file : str
        Path to the input genome FASTA file (nucleotides).
    prefix : str
        Prefix to prepend to each sequence header.
    outdir : str
        Output directory for the prefixed file.

    Example
    -------
    >>> InsertPrefix("genome.fasta", "Aaeg", "results")  # doctest: +SKIP
    """

    def __init__(self, input_file: str, prefix: str, outdir: str) -> None:
        self.input_file = input_file
        self.prefix = prefix
        self.outdir = outdir

        self.insert_prefix()

    def insert_prefix(self) -> None:
        """Stream the FASTA, prepending ``{prefix}/`` to each header line."""
        output_file = Path(self.outdir) / f"{self.prefix}.rn"
        with open(self.input_file) as reader, open(output_file, "w") as writer:
            for line in reader:
                if line.startswith(">"):
                    line = f">{self.prefix}/{line[1:]}"
                writer.write(line)
