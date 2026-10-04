"""Taxonomy assignment: join the surviving hits to the ``-mt`` metadata.

Produces the per-hit signature. The per-element table is produced by
:mod:`eefinder.stages.merging`, since an element only exists once neighbouring
fragments have been merged.
"""

from __future__ import annotations
from dataclasses import dataclass
from eefinder.models import TAXONOMY_SIGNATURE
from eefinder.get_taxonomy import GetTaxonomy
from eefinder.log import logger
from eefinder.stages.base import Stage, StageOutputs


@dataclass
class TaxonomyAssignmentOutputs(StageOutputs):
    """Output of :class:`TaxonomyAssignment`.

    Attributes
    ----------
    taxonomy_signature : str
        Validated hits joined to the metadata, one row per hit. Carries the
        accession, the coordinates, the sense and the taxonomic columns that
        Stage V merges on.
    """

    taxonomy_signature: str


class TaxonomyAssignment(Stage):
    """Join the surviving hits to the reference metadata table.

    Parameters
    ----------
    ee_hits_validated : str
        Surviving hits from
        :class:`~eefinder.stages.filtering.PutativeElementsFilter`.
    dbmetadata : str
        Protein metadata CSV (``-mt``). Its header is validated by Stage 0 and
        again here, since the stage can be run on its own.
    outdir, prefix : str
        Run output directory and prefix.

    Example
    -------
    >>> result = TaxonomyAssignment(                      # doctest: +SKIP
    ...     ee_hits_validated="out/run.ee_hits.validated.tsv",
    ...     dbmetadata="db/virus.csv",
    ...     outdir="out",
    ...     prefix="run",
    ... ).run()
    >>> result.taxonomy_signature                         # doctest: +SKIP
    'out/run.taxonomy_signature.csv'
    """

    name = "taxonomy"
    title = "Taxonomy assignment"
    stage_id = "05"

    models = {"taxonomy_signature": TAXONOMY_SIGNATURE}

    def __init__(
        self,
        ee_hits_validated: str,
        dbmetadata: str,
        outdir: str,
        prefix: str,
    ) -> None:
        super().__init__(outdir, prefix)
        self.ee_hits_validated = ee_hits_validated
        self.dbmetadata = dbmetadata

    def _execute(self):
        self._require(self.ee_hits_validated, self.dbmetadata)

        logger.debug(f"GetTaxonomy: joining hits to metadata {self.dbmetadata}")
        GetTaxonomy(self.ee_hits_validated, self.dbmetadata)
        self._rename(f"{self.ee_hits_validated}.tax", self.paths.taxonomy_signature)

        outputs = TaxonomyAssignmentOutputs(
            step_info=None,
            taxonomy_signature=self.paths.taxonomy_signature,
        )
        return outputs, "Performed initial taxonomy"
