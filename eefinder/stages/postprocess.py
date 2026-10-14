"""Post-processing: repeat filtering, overlap resolution and GFF3/protein output.

Each narrows or re-renders the merged elements and their taxonomy table, so
they share one stage.
"""

from __future__ import annotations
import os
import shutil
from dataclasses import dataclass
from typing import Optional, Sequence
from eefinder.cds import element_cds, write_faa
from eefinder.clean_data import MaskClean
from eefinder.models import ELEMENT_FASTA, ELEMENT_PROTEINS, ELEMENT_TABLE
from eefinder.get_taxonomy import GetCleanedTaxonomy
from eefinder.gff import WriteGFF3
from eefinder.log import logger
from eefinder.overlap import FilterOverlap
from eefinder.stages.base import Stage, StageOutputs
from eefinder.tag_elements import TagElements


@dataclass
class PostProcessingOutputs(StageOutputs):
    """Outputs of :class:`PostProcessing`.

    Attributes
    ----------
    ee_elements_resolved : str
        The elements left after overlap resolution -- a copy of the input
        under ``--overlap keep``. Published as ``PREFIX.EEs.fa``.
    ee_elements_resolved_tax : str
        Their taxonomy table -- published as ``PREFIX.EEs.tax.tsv``.
    ee_elements_gff3 : str
        GFF3 annotation of the resolved elements, with their CDSs.
    ee_elements_faa : str
        The proteins of those CDSs. Published as ``PREFIX.EEs.faa``.
    ee_elements_cleaned : str or None
        Mask-cleaned element sequences, when ``--clean_masked`` was given.
    ee_elements_cleaned_tax : str or None
        Taxonomy table restricted to those elements.
    ee_elements_cleaned_gff3 : str or None
        GFF3 annotation of those elements.
    removed_fasta, removed_tax : str or None
        Elements dropped by overlap resolution, preserved for inspection.
    """

    ee_elements_resolved: str
    ee_elements_resolved_tax: str
    ee_elements_gff3: str
    ee_elements_faa: str
    ee_elements_cleaned: Optional[str] = None
    ee_elements_cleaned_tax: Optional[str] = None
    ee_elements_cleaned_gff3: Optional[str] = None
    removed_fasta: Optional[str] = None
    removed_tax: Optional[str] = None


class PostProcessing(Stage):
    """Apply the repeat filter, resolve overlaps and write the GFF3 and proteins.

    Parameters
    ----------
    ee_elements : str
        Merged element sequences from
        :class:`~eefinder.stages.merging.MergeFragmentedElements`.
    ee_elements_tax : str
        The matching per-element taxonomy table.
    taxonomy_signature : str
        The per-hit table from
        :class:`~eefinder.stages.taxonomy.TaxonomyAssignment`; its hits become
        the elements' CDSs (:mod:`eefinder.cds`).
    outdir, prefix : str
        Run output directory and prefix.
    clean_masked : bool
        Emit the mask-cleaned second set of outputs (``--clean_masked``).
    mask_per : int
        Soft-masked percentage above which an element is dropped from that set
        (``--mask_per``).
    overlap : str
        ``"keep"``, ``"longest"`` or ``"targets"`` (``--overlap``).
    target_families, non_target_families : sequence of str
        Keep-list / drop-list for ``overlap="targets"``; exactly one must be
        non-empty.
    analysis : str
        ``"virus"`` or ``"bacteria"``, selecting the GFF3 feature type
        (``--analysis``).

    Example
    -------
    >>> result = PostProcessing(                           # doctest: +SKIP
    ...     ee_elements="out/run.elements.fa",
    ...     ee_elements_tax="out/run.elements.tax.tsv",
    ...     taxonomy_signature="out/run.taxonomy_signature.csv",
    ...     outdir="out",
    ...     prefix="run",
    ...     overlap="longest",
    ... ).run()
    >>> result.ee_elements_gff3                            # doctest: +SKIP
    'out/run.elements.gff3'
    """

    name = "postprocess"
    title = "Post-processing"
    stage_id = "07"

    models = {
        "ee_elements_resolved": ELEMENT_FASTA,
        "ee_elements_resolved_tax": ELEMENT_TABLE,
        "ee_elements_faa": ELEMENT_PROTEINS,
        "ee_elements_cleaned": ELEMENT_FASTA,
        "ee_elements_cleaned_tax": ELEMENT_TABLE,
    }

    def __init__(
        self,
        ee_elements: str,
        ee_elements_tax: str,
        taxonomy_signature: str,
        outdir: str,
        prefix: str,
        clean_masked: bool = False,
        mask_per: int = 50,
        overlap: str = "keep",
        target_families: Sequence[str] = (),
        non_target_families: Sequence[str] = (),
        analysis: str = "virus",
    ) -> None:
        super().__init__(outdir, prefix)
        self.ee_elements = ee_elements
        self.ee_elements_tax = ee_elements_tax
        self.taxonomy_signature = taxonomy_signature
        self.clean_masked = clean_masked
        self.mask_per = mask_per
        self.overlap = overlap
        self.target_families = list(target_families)
        self.non_target_families = list(non_target_families)
        self.analysis = analysis

    def _execute(self):
        self._require(self.ee_elements, self.ee_elements_tax, self.taxonomy_signature)
        paths = self.paths
        actions = []

        cleaned = cleaned_tax = cleaned_gff3 = None
        if self.clean_masked:
            logger.info("Cleaning elements")
            logger.debug(
                f"MaskClean: removing EEs with > {self.mask_per}% soft-masking"
            )
            MaskClean(self.ee_elements, self.mask_per)
            cleaned = self._rename(f"{self.ee_elements}.cl", paths.ee_elements_cleaned)
            GetCleanedTaxonomy(cleaned, self.ee_elements_tax)
            cleaned_tax = self._rename(f"{cleaned}.tax", paths.ee_elements_cleaned_tax)
            TagElements(cleaned_tax)
            actions.append(
                f"EEs with {self.mask_per} percent of lower-case letters are removed"
            )

        # Overlap resolution rewrites its files in place, so it works on copies:
        # the inputs may belong to another run (e.g. a symlinked Nextflow input).
        resolved = paths.ee_elements_resolved
        resolved_tax = paths.ee_elements_resolved_tax
        shutil.copyfile(self.ee_elements, resolved)
        shutil.copyfile(self.ee_elements_tax, resolved_tax)

        removed_fasta = removed_tax = None
        if self.overlap != "keep":
            logger.info(f"Filtering overlaping elements ({self.overlap})")
            os.makedirs(paths.removed_dir, exist_ok=True)
            removed_fasta = f"{paths.removed_dir}/{self.prefix}.EEs.removed.fa"
            removed_tax = f"{paths.removed_dir}/{self.prefix}.EEs.removed.tax.tsv"
            logger.debug(
                f"FilterOverlap strategy={self.overlap} "
                f"target_families={self.target_families} "
                f"non_target_families={self.non_target_families}; "
                f"removed elements -> {paths.removed_dir}"
            )
            FilterOverlap(
                resolved,
                resolved_tax,
                self.overlap,
                self.target_families,
                removed_fasta,
                removed_tax,
                non_target_families=self.non_target_families,
            )
            if self.clean_masked:
                FilterOverlap(
                    cleaned,
                    cleaned_tax,
                    self.overlap,
                    self.target_families,
                    f"{paths.removed_dir}/{self.prefix}.EEs.cleaned.removed.fa",
                    f"{paths.removed_dir}/{self.prefix}.EEs.cleaned.removed.tax.tsv",
                    non_target_families=self.non_target_families,
                )
            actions.append(
                f"resolved overlaping elements with the '{self.overlap}' strategy, "
                f"saving the filtered-out elements to {paths.removed_dir}"
            )

        logger.info("Generating GFF3 annotation and element proteins")
        logger.debug(f"WriteGFF3 (analysis={self.analysis}) for prefix {self.prefix!r}")
        cds = element_cds(resolved, resolved_tax, self.taxonomy_signature, self.prefix)
        WriteGFF3(
            resolved_tax,
            paths.ee_elements_gff3,
            prefix=self.prefix,
            analysis=self.analysis,
            cds=cds,
        )
        write_faa(cds, paths.ee_elements_faa)
        if self.clean_masked:
            cleaned_gff3 = paths.ee_elements_cleaned_gff3
            WriteGFF3(
                cleaned_tax,
                cleaned_gff3,
                prefix=self.prefix,
                analysis=self.analysis,
                cds=element_cds(
                    cleaned, cleaned_tax, self.taxonomy_signature, self.prefix
                ),
            )
        actions.append(
            "generated GFF3 annotation of endogenous elements and their CDS proteins"
        )

        outputs = PostProcessingOutputs(
            step_info=None,
            ee_elements_resolved=resolved,
            ee_elements_resolved_tax=resolved_tax,
            ee_elements_gff3=paths.ee_elements_gff3,
            ee_elements_faa=paths.ee_elements_faa,
            ee_elements_cleaned=cleaned,
            ee_elements_cleaned_tax=cleaned_tax,
            ee_elements_cleaned_gff3=cleaned_gff3,
            removed_fasta=removed_fasta,
            removed_tax=removed_tax,
        )
        return outputs, "; ".join(actions).capitalize() + "."
