"""Putative elements filter: drop candidates that match the host better.

Reference viral and bacterial proteins share conserved domains with host
proteins, so a locus can match the database without being donor-derived.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
from eefinder.compare_results import CompareResults
from eefinder.models import FILTERED_HIT_TABLE, HIT_TABLE
from eefinder.filter_table import FilterTable
from eefinder.log import logger
from eefinder.search_methods import SearchMethod, resolve_search_method
from eefinder.stages.base import Stage, StageOutputs


@dataclass
class PutativeElementsFilterOutputs(StageOutputs):
    """Outputs of :class:`PutativeElementsFilter`.

    Attributes
    ----------
    host_hits : str
        Raw tabular hits of the candidates against the host-gene baits.
    host_hits_filtered : str
        The same hits with redundant alignments collapsed.
    ee_hits_validated : str
        The candidate hits that survived the bitscore comparison -- the input
        to taxonomy assignment.
    """

    host_hits: str
    host_hits_filtered: str
    ee_hits_validated: str


class PutativeElementsFilter(Stage):
    """Drop candidates that match the host better than the reference database.

    Parameters
    ----------
    candidates : str
        Candidate sequences from
        :class:`~eefinder.stages.alignment.SequenceAlignment`.
    ee_hits_filtered : str
        The filtered database hits for those candidates, also from Stage II.
    hostgenesbaits : str
        Host-gene bait FASTA (``-bt``), already indexed.
    outdir, prefix : str
        Run output directory and prefix.
    mode, threads, translation_method, range_junction, search_method
        As in :class:`~eefinder.stages.alignment.SequenceAlignment`; the search
        settings are deliberately the same ones, so the two searches of a run
        can never diverge.

    Example
    -------
    >>> result = PutativeElementsFilter(                  # doctest: +SKIP
    ...     candidates="out/run.ee_candidates.fa",
    ...     ee_hits_filtered="out/run.ee_hits.filtered.tsv",
    ...     hostgenesbaits="db/host.fa",
    ...     outdir="out",
    ...     prefix="run",
    ... ).run()
    >>> result.ee_hits_validated                          # doctest: +SKIP
    'out/run.ee_hits.validated.tsv'
    """

    name = "filter"
    title = "Putative elements filter"
    stage_id = "04"

    models = {
        "host_hits": HIT_TABLE,
        "host_hits_filtered": FILTERED_HIT_TABLE,
        "ee_hits_validated": FILTERED_HIT_TABLE,
    }

    def __init__(
        self,
        candidates: str,
        ee_hits_filtered: str,
        hostgenesbaits: str,
        outdir: str,
        prefix: str,
        mode: str = "blastx",
        threads: int = 1,
        translation_method: str = "default",
        range_junction: int = 100,
        search_method: "Optional[SearchMethod]" = None,
    ) -> None:
        super().__init__(outdir, prefix)
        self.search_method = search_method or resolve_search_method(
            mode, translation_method
        )
        self.candidates = candidates
        self.ee_hits_filtered = ee_hits_filtered
        self.hostgenesbaits = hostgenesbaits
        self.mode = mode
        self.threads = threads
        self.translation_method = translation_method
        self.range_junction = range_junction

    def _execute(self):
        self._require(self.candidates, self.ee_hits_filtered, self.hostgenesbaits)
        paths = self.paths

        logger.debug(f"Host-bait search method: {self.search_method.describe()}")
        self.search_method.run(
            self.candidates, self.hostgenesbaits, self.threads, paths.host_hits
        )

        logger.debug(
            f"FilterTable HOST: {paths.host_hits} "
            f"(range_junction={self.range_junction})"
        )
        FilterTable(paths.host_hits, self.range_junction, "HOST", self.outdir)
        self._rename(f"{paths.host_hits}.filtred", paths.host_hits_filtered)

        logger.debug("CompareResults: dropping EEs that hit host baits harder")
        CompareResults(self.ee_hits_filtered, paths.host_hits_filtered)
        self._rename(f"{paths.host_hits_filtered}.concat.nr", paths.ee_hits_validated)

        outputs = PutativeElementsFilterOutputs(
            step_info=None,
            host_hits=paths.host_hits,
            host_hits_filtered=paths.host_hits_filtered,
            ee_hits_validated=paths.ee_hits_validated,
        )
        message = (
            f"Filter step based on similarity analysis with {self.mode} was "
            f"performed using {self.candidates} against {self.hostgenesbaits}. "
            f"Matches against same subject sequence in a {self.range_junction}nt "
            "range junction are filtered, mantaining the one with the greatest "
            "bitscore. The results are compared, and the putative EEs with the "
            "greatest bitscore on host genes baits database are removed"
        )
        return outputs, message
