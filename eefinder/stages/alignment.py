"""Sequence alignment: search the genome and extract the candidate elements."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
from eefinder.bed import GetFasta
from eefinder.models import (
    CANDIDATE_BED,
    CANDIDATE_FASTA,
    FILTERED_HIT_TABLE,
    HIT_TABLE,
)
from eefinder.filter_table import FilterTable
from eefinder.log import logger
from eefinder.search_methods import SearchMethod, resolve_search_method
from eefinder.stages.base import Stage, StageOutputs


@dataclass
class SequenceAlignmentOutputs(StageOutputs):
    """Outputs of :class:`SequenceAlignment`.

    Attributes
    ----------
    ee_hits : str
        Raw tabular hits against the reference database.
    ee_hits_filtered : str
        The same hits with redundant alignments over one locus collapsed.
    ee_candidates_bed : str
        Candidate intervals on the cleaned genome.
    ee_candidates : str
        Candidate sequences, the query for Stage III.
    """

    ee_hits: str
    ee_hits_filtered: str
    ee_candidates_bed: str
    ee_candidates: str


class SequenceAlignment(Stage):
    """Search the genome against the reference proteins and cut out candidates.

    Parameters
    ----------
    genome : str
        Cleaned genome from :class:`~eefinder.stages.cleaning.DataCleaning`.
    database : str
        Protein database FASTA (``-db``), already indexed.
    outdir, prefix : str
        Run output directory and prefix.
    mode : str
        ``"blastx"`` or a DIAMOND sensitivity (``--mode``).
    threads : int
        Threads for the search (``--threads``).
    translation_method : str
        ``"default"`` six-frame search, or ``gv``/``rv``/``gv-rv`` protein
        prediction (``--translation_method``).
    range_junction : int
        Window for collapsing redundant hits (``--range_junction``).
    search_method : SearchMethod, optional
        Produce the hits with this method instead of the one implied by
        ``mode``/``translation_method``. This is the extension point: any object
        satisfying :data:`~eefinder.models.HIT_TABLE` can be passed in, and
        the rest of the pipeline is unaffected.

    Example
    -------
    >>> result = SequenceAlignment(                   # doctest: +SKIP
    ...     genome="out/run.cleaned_genome.fa",
    ...     database="db/virus.fa",
    ...     outdir="out",
    ...     prefix="run",
    ...     threads=8,
    ... ).run()
    >>> result.ee_candidates                          # doctest: +SKIP
    'out/run.ee_candidates.fa'
    """

    name = "align"
    title = "Sequence alignment"
    stage_id = "03"

    models = {
        "ee_hits": HIT_TABLE,
        "ee_hits_filtered": FILTERED_HIT_TABLE,
        "ee_candidates_bed": CANDIDATE_BED,
        "ee_candidates": CANDIDATE_FASTA,
    }

    def __init__(
        self,
        genome: str,
        database: str,
        outdir: str,
        prefix: str,
        mode: str = "blastx",
        threads: int = 1,
        translation_method: str = "default",
        range_junction: int = 100,
        search_method: "Optional[SearchMethod]" = None,
    ) -> None:
        super().__init__(outdir, prefix)
        self.genome = genome
        self.database = database
        self.mode = mode
        self.threads = threads
        self.translation_method = translation_method
        self.range_junction = range_junction
        self.search_method = search_method or resolve_search_method(
            mode, translation_method
        )

    def _execute(self):
        self._require(self.genome, self.database)
        paths = self.paths

        logger.debug(f"Search method: {self.search_method.describe()}")
        self.search_method.run(self.genome, self.database, self.threads, paths.ee_hits)

        logger.debug(
            f"FilterTable EE: {paths.ee_hits} "
            f"(range_junction={self.range_junction})"
        )
        FilterTable(paths.ee_hits, self.range_junction, "EE", self.outdir)
        self._rename(f"{paths.ee_hits}.filtred", paths.ee_hits_filtered)
        self._rename(f"{paths.ee_hits}.filtred.bed", paths.ee_candidates_bed)

        logger.debug(f"GetFasta putative EEs from {paths.ee_candidates_bed}")
        GetFasta(self.genome, paths.ee_candidates_bed, paths.ee_candidates)

        outputs = SequenceAlignmentOutputs(
            step_info=None,
            ee_hits=paths.ee_hits,
            ee_hits_filtered=paths.ee_hits_filtered,
            ee_candidates_bed=paths.ee_candidates_bed,
            ee_candidates=paths.ee_candidates,
        )
        message = (
            f"Similarity analysis with {self.search_method.describe()} was "
            f"performed using {self.genome} against {self.database}."
            f"Matches against same subject sequence in a {self.range_junction}nt "
            "range junction are filtered, mantaining the one with the greatest "
            "bitscore. The surviving intervals were extracted as putative EEs."
        )
        return outputs, message
