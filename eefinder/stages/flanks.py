"""Flanks extraction: widen each element by ``--flank`` nt and cut it out."""

from __future__ import annotations
from dataclasses import dataclass
from eefinder.bed import BedFlank, GetBed, GetFasta
from eefinder.models import CANDIDATE_FASTA
from eefinder.get_length import GetLength
from eefinder.log import logger
from eefinder.stages.base import Stage, StageOutputs


@dataclass
class FlanksExtractionOutputs(StageOutputs):
    """Outputs of :class:`FlanksExtraction`.

    Attributes
    ----------
    genome_lengths : str
        ``<contig>\\t<length>`` index that ``bedtools slop`` needs in order not
        to run off the end of a contig.
    flanks_bed : str
        Element intervals.
    flanks_slop_bed : str
        The same intervals widened by ``--flank`` nt on each side.
    flanks : str
        Element sequences plus flanks -- published as ``PREFIX.EEs.flanks.fa``.
    """

    genome_lengths: str
    flanks_bed: str
    flanks_slop_bed: str
    flanks: str


class FlanksExtraction(Stage):
    """Extract each element plus its flanking regions.

    Parameters
    ----------
    ee_elements : str
        Merged element sequences, whose headers carry the coordinates.
    genome : str
        Cleaned genome the flanks are cut from.
    outdir, prefix : str
        Run output directory and prefix.
    flank : int
        Flank length in nt on each side (``--flank``).

    Example
    -------
    >>> result = FlanksExtraction(                       # doctest: +SKIP
    ...     ee_elements="out/run.elements.fa",
    ...     genome="out/run.cleaned_genome.fa",
    ...     outdir="out",
    ...     prefix="run",
    ...     flank=10000,
    ... ).run()
    >>> result.flanks                                    # doctest: +SKIP
    'out/run.flanks.fa'
    """

    name = "flanks"
    title = "Flanks extraction"
    stage_id = "08"

    models = {"flanks": CANDIDATE_FASTA}

    def __init__(
        self,
        ee_elements: str,
        genome: str,
        outdir: str,
        prefix: str,
        flank: int = 10000,
    ) -> None:
        super().__init__(outdir, prefix)
        self.ee_elements = ee_elements
        self.genome = genome
        self.flank = flank

    def _execute(self):
        self._require(self.ee_elements, self.genome)
        paths = self.paths

        logger.debug(f"Extracting {self.flank} nt flanks around each EE")
        GetLength(self.genome)
        self._rename(f"{self.genome}.rn.fmt.lenght", paths.genome_lengths)

        GetBed(self.ee_elements)
        self._rename(f"{self.ee_elements}.bed", paths.flanks_bed)

        BedFlank(paths.flanks_bed, paths.genome_lengths, self.flank)
        self._rename(f"{paths.flanks_bed}.flank", paths.flanks_slop_bed)

        GetFasta(self.genome, paths.flanks_slop_bed, paths.flanks)

        outputs = FlanksExtractionOutputs(
            step_info=None,
            genome_lengths=paths.genome_lengths,
            flanks_bed=paths.flanks_bed,
            flanks_slop_bed=paths.flanks_slop_bed,
            flanks=paths.flanks,
        )
        return outputs, f"Extracted {self.flank}nt of each flanking region of EEs."
