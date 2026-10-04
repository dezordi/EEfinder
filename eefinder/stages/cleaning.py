"""Data cleaning: drop short contigs and prefix the FASTA headers."""

from __future__ import annotations
from dataclasses import dataclass
from eefinder.models import GENOME_FASTA
from eefinder.log import logger
from eefinder.prepare_data import PrepareGenome
from eefinder.stages.base import Stage, StageOutputs


@dataclass
class DataCleaningOutputs(StageOutputs):
    """Output of :class:`DataCleaning`.

    Attributes
    ----------
    cleaned_genome : str
        Prefixed, length-filtered genome. Every later stage searches or slices
        this file, never the original input.
    kept, total : int
        Contigs retained and seen, for the run log.
    """

    cleaned_genome: str
    kept: int = 0
    total: int = 0


class DataCleaning(Stage):
    """Prefix the FASTA headers and drop contigs below a length cutoff.

    Parameters
    ----------
    genome_file : str
        Input genome FASTA (``-in``).
    outdir, prefix : str
        Run output directory and prefix.
    length : int
        Minimum contig length to keep (``--length``).

    Example
    -------
    >>> result = DataCleaning(                       # doctest: +SKIP
    ...     genome_file="genome.fa", outdir="out", prefix="run", length=1000
    ... ).run()
    >>> result.cleaned_genome                        # doctest: +SKIP
    'out/run.cleaned_genome.fa'
    """

    name = "clean"
    title = "Data cleaning"
    stage_id = "02"

    models = {"cleaned_genome": GENOME_FASTA}

    def __init__(
        self,
        genome_file: str,
        outdir: str,
        prefix: str,
        length: int = 10000,
    ) -> None:
        super().__init__(outdir, prefix)
        self.genome_file = genome_file
        self.length = length

    def _execute(self):
        self._require(self.genome_file)
        logger.debug(
            f"PrepareGenome: {self.genome_file} -> {self.paths.cleaned_genome} "
            f"(prefixing headers, dropping contigs < {self.length} nt)"
        )
        prepared = PrepareGenome(
            self.genome_file, self.prefix, self.outdir, self.length
        )
        self._rename(f"{self.outdir}/{self.prefix}.rn.fmt", self.paths.cleaned_genome)
        logger.info(
            f"Prepared {prepared.kept} of {prepared.total} contig(s) "
            f"(>= {self.length} nt)"
        )

        outputs = DataCleaningOutputs(
            step_info=None,
            cleaned_genome=self.paths.cleaned_genome,
            kept=prepared.kept,
            total=prepared.total,
        )
        message = (
            f"{self.prefix} prefix included in {self.genome_file} sequences "
            f"header and sequences bellow than {self.length} nt are removed "
            f"from {self.genome_file}."
        )
        return outputs, message
