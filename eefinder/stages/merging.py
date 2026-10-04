"""Merging fragmented elements, and the per-element taxonomy table.

One integration often survives as several fragments. Neighbouring fragments on
the same strand sharing a taxon at ``--merge_level`` within ``--limit`` nt are
merged into one element.
"""

from __future__ import annotations
from dataclasses import dataclass
from eefinder.bed import GetAnnotBed, GetFasta, MergeBed, RemoveAnnotation
from eefinder.models import ELEMENT_FASTA, ELEMENT_TABLE
from eefinder.get_taxonomy import GetFinalTaxonomy
from eefinder.log import logger
from eefinder.stages.base import Stage, StageOutputs
from eefinder.tag_elements import TagElements


@dataclass
class MergeFragmentedElementsOutputs(StageOutputs):
    """Outputs of :class:`MergeFragmentedElements`.

    Attributes
    ----------
    merge_bed : str
        Per-hit intervals annotated with the merge key (contig, taxon, sense).
    merge_bed_merged : str
        Result of ``bedtools merge`` over those intervals.
    merge_elements_bed : str
        Merged intervals with the annotation stripped back off.
    ee_elements : str
        Merged element sequences -- published as ``PREFIX.EEs.fa``.
    ee_elements_tax : str
        Per-element taxonomy table, with ``Average_pident`` and the overlap
        tags -- published as ``PREFIX.EEs.tax.tsv``.
    """

    merge_bed: str
    merge_bed_merged: str
    merge_elements_bed: str
    ee_elements: str
    ee_elements_tax: str


class MergeFragmentedElements(Stage):
    """Merge neighbouring same-taxon fragments and build the element table.

    Parameters
    ----------
    taxonomy_signature : str
        Per-hit signature from
        :class:`~eefinder.stages.taxonomy.TaxonomyAssignment`.
    genome : str
        Cleaned genome, re-sliced to produce the merged element sequences.
    outdir, prefix : str
        Run output directory and prefix.
    limit : int
        Maximum gap between fragments to merge, in nt (``--limit``).
    merge_level : str
        ``"family"`` or ``"genus"`` -- the rank at which two fragments count as
        the same taxon (``--merge_level``).

    Example
    -------
    >>> result = MergeFragmentedElements(                  # doctest: +SKIP
    ...     taxonomy_signature="out/run.taxonomy_signature.csv",
    ...     genome="out/run.cleaned_genome.fa",
    ...     outdir="out",
    ...     prefix="run",
    ...     limit=100,
    ...     merge_level="family",
    ... ).run()
    >>> result.ee_elements_tax                             # doctest: +SKIP
    'out/run.elements.tax.tsv'
    """

    name = "merge"
    title = "Merging fragmented elements"
    stage_id = "06"

    models = {
        "ee_elements": ELEMENT_FASTA,
        "ee_elements_tax": ELEMENT_TABLE,
    }

    def __init__(
        self,
        taxonomy_signature: str,
        genome: str,
        outdir: str,
        prefix: str,
        limit: int = 1,
        merge_level: str = "family",
    ) -> None:
        super().__init__(outdir, prefix)
        self.taxonomy_signature = taxonomy_signature
        self.genome = genome
        self.limit = limit
        self.merge_level = merge_level

    def _execute(self):
        self._require(self.taxonomy_signature, self.genome)
        paths = self.paths

        logger.debug(
            f"Merging truncated elements by {self.merge_level} "
            f"within {self.limit} nt"
        )
        GetAnnotBed(self.taxonomy_signature, self.merge_level)
        self._rename(f"{self.taxonomy_signature}.bed", paths.merge_bed)

        MergeBed(paths.merge_bed, str(self.limit))
        self._rename(f"{paths.merge_bed}.merge", paths.merge_bed_merged)

        RemoveAnnotation(paths.merge_bed_merged)
        self._rename(f"{paths.merge_bed_merged}.fmt", paths.merge_elements_bed)

        GetFasta(self.genome, paths.merge_elements_bed, paths.ee_elements)

        logger.debug("GetFinalTaxonomy + TagElements (Average_pident, overlap tags)")
        GetFinalTaxonomy(paths.merge_elements_bed, self.taxonomy_signature)
        self._rename(f"{paths.merge_elements_bed}.fa.tax", paths.ee_elements_tax)
        TagElements(paths.ee_elements_tax)

        outputs = MergeFragmentedElementsOutputs(
            step_info=None,
            merge_bed=paths.merge_bed,
            merge_bed_merged=paths.merge_bed_merged,
            merge_elements_bed=paths.merge_elements_bed,
            ee_elements=paths.ee_elements,
            ee_elements_tax=paths.ee_elements_tax,
        )
        message = (
            f"Merge EEs near of {self.limit}nt based on {self.merge_level} "
            "taxonomy information, and assembled the per-element taxonomy table."
        )
        return outputs, message
