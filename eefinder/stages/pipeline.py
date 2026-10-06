"""The whole screening pipeline.

:class:`ScreeningPipeline` threads each stage's outputs into the next one's
inputs, publishes the ``PREFIX.EEs.*`` outputs and archives the intermediates.
"""

from __future__ import annotations
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence
from eefinder.log import logger
from eefinder.stages.alignment import SequenceAlignment
from eefinder.stages.base import Stage, StageOutputs
from eefinder.stages.cleaning import DataCleaning
from eefinder.stages.filtering import PutativeElementsFilter
from eefinder.stages.flanks import FlanksExtraction
from eefinder.stages.merging import MergeFragmentedElements
from eefinder.stages.postprocess import PostProcessing
from eefinder.stages.prepare import PrepareInputs
from eefinder.stages.taxonomy import TaxonomyAssignment
from eefinder.utils import StepInfo


@dataclass
class ScreeningPipelineOutputs(StageOutputs):
    """Published outputs of a full run.

    Attributes
    ----------
    final_fasta, final_tax, final_gff3, final_flanks : str
        The four main outputs.
    final_cleaned_fasta, final_cleaned_tax, final_cleaned_gff3 : str or None
        Their ``--clean_masked`` counterparts.
    step_infos : list of StepInfo
        One record per stage, in execution order, for the run log.
    """

    final_fasta: str
    final_tax: str
    final_gff3: str
    final_flanks: str
    final_cleaned_fasta: Optional[str] = None
    final_cleaned_tax: Optional[str] = None
    final_cleaned_gff3: Optional[str] = None
    step_infos: "list[StepInfo]" = field(default_factory=list)


class ScreeningPipeline(Stage):
    """Run every stage of the screening pipeline in order.

    Parameters
    ----------
    genome_file, database, dbmetadata, hostgenesbaits : str
        The four required inputs (``-in``, ``-db``, ``-mt``, ``-bt``).
    outdir, prefix : str
        Run output directory and prefix.
    mode, threads, translation_method, range_junction, length
        Search and cleaning settings.
    limit, merge_level, flank
        Element assembly and flank settings.
    clean_masked, mask_per, overlap, target_families, non_target_families, analysis
        Post-processing settings.
    index_databases : bool
        Build the search indexes in stage 0.
    removetmp : bool
        Delete the intermediates instead of archiving them under
        ``tmp_files/``.

    Example
    -------
    >>> result = ScreeningPipeline(                        # doctest: +SKIP
    ...     genome_file="genome.fa",
    ...     database="db/virus.fa",
    ...     dbmetadata="db/virus.csv",
    ...     hostgenesbaits="db/host.fa",
    ...     outdir="out",
    ...     prefix="run",
    ... ).run()
    >>> result.final_tax                                   # doctest: +SKIP
    'out/run.EEs.tax.tsv'
    """

    name = "all"
    title = "Full screening pipeline"
    stage_id = "01-08"

    def __init__(
        self,
        genome_file: str,
        database: str,
        dbmetadata: str,
        hostgenesbaits: str,
        outdir: str,
        prefix: str,
        mode: str = "blastx",
        threads: int = 1,
        translation_method: str = "default",
        range_junction: int = 100,
        length: int = 10000,
        limit: int = 1,
        merge_level: str = "family",
        flank: int = 10000,
        clean_masked: bool = False,
        mask_per: int = 50,
        overlap: str = "keep",
        target_families: Sequence[str] = (),
        non_target_families: Sequence[str] = (),
        analysis: str = "virus",
        index_databases: bool = False,
        removetmp: bool = False,
    ) -> None:
        super().__init__(outdir, prefix)
        self.genome_file = genome_file
        self.database = database
        self.dbmetadata = dbmetadata
        self.hostgenesbaits = hostgenesbaits
        self.mode = mode
        self.threads = threads
        self.translation_method = translation_method
        self.range_junction = range_junction
        self.length = length
        self.limit = limit
        self.merge_level = merge_level
        self.flank = flank
        self.clean_masked = clean_masked
        self.mask_per = mask_per
        self.overlap = overlap
        self.target_families = list(target_families)
        self.non_target_families = list(non_target_families)
        self.analysis = analysis
        self.index_databases = index_databases
        self.removetmp = removetmp

    def _execute(self):
        paths = self.paths
        step_infos: "list[StepInfo]" = []

        def record(outputs: StageOutputs) -> StageOutputs:
            step_infos.append(outputs.step_info)
            return outputs

        logger.info("Preparing inputs")
        record(
            PrepareInputs(
                database=self.database,
                dbmetadata=self.dbmetadata,
                hostgenesbaits=self.hostgenesbaits,
                outdir=self.outdir,
                prefix=self.prefix,
                mode=self.mode,
                threads=self.threads,
                index_databases=self.index_databases,
            ).run()
        )

        logger.info("Cleaning input data")
        cleaned = record(
            DataCleaning(
                genome_file=self.genome_file,
                outdir=self.outdir,
                prefix=self.prefix,
                length=self.length,
            ).run()
        )

        logger.info("Running similarity search")
        aligned = record(
            SequenceAlignment(
                genome=cleaned.cleaned_genome,
                database=self.database,
                outdir=self.outdir,
                prefix=self.prefix,
                mode=self.mode,
                threads=self.threads,
                translation_method=self.translation_method,
                range_junction=self.range_junction,
            ).run()
        )

        logger.info("Running filter steps")
        filtered = record(
            PutativeElementsFilter(
                candidates=aligned.ee_candidates,
                ee_hits_filtered=aligned.ee_hits_filtered,
                hostgenesbaits=self.hostgenesbaits,
                outdir=self.outdir,
                prefix=self.prefix,
                mode=self.mode,
                threads=self.threads,
                translation_method=self.translation_method,
                range_junction=self.range_junction,
            ).run()
        )

        logger.info("Assigning taxonomy")
        assigned = record(
            TaxonomyAssignment(
                ee_hits_validated=filtered.ee_hits_validated,
                dbmetadata=self.dbmetadata,
                outdir=self.outdir,
                prefix=self.prefix,
            ).run()
        )

        logger.info("Merging truncated elements")
        merged = record(
            MergeFragmentedElements(
                taxonomy_signature=assigned.taxonomy_signature,
                genome=cleaned.cleaned_genome,
                outdir=self.outdir,
                prefix=self.prefix,
                limit=self.limit,
                merge_level=self.merge_level,
            ).run()
        )

        logger.info("Post-processing elements")
        processed = record(
            PostProcessing(
                ee_elements=merged.ee_elements,
                ee_elements_tax=merged.ee_elements_tax,
                outdir=self.outdir,
                prefix=self.prefix,
                clean_masked=self.clean_masked,
                mask_per=self.mask_per,
                overlap=self.overlap,
                target_families=self.target_families,
                non_target_families=self.non_target_families,
                analysis=self.analysis,
            ).run()
        )

        logger.info("Extracting flanking regions")
        flanked = record(
            FlanksExtraction(
                ee_elements=merged.ee_elements,
                genome=cleaned.cleaned_genome,
                outdir=self.outdir,
                prefix=self.prefix,
                flank=self.flank,
            ).run()
        )

        logger.info("Organizing final outputs")
        self._publish(merged, processed, flanked)
        self._archive()

        outputs = ScreeningPipelineOutputs(
            step_info=None,
            final_fasta=paths.final_fasta,
            final_tax=paths.final_tax,
            final_gff3=paths.final_gff3,
            final_flanks=paths.final_flanks,
            final_cleaned_fasta=(
                paths.final_cleaned_fasta if self.clean_masked else None
            ),
            final_cleaned_tax=(paths.final_cleaned_tax if self.clean_masked else None),
            final_cleaned_gff3=(
                paths.final_cleaned_gff3 if self.clean_masked else None
            ),
            step_infos=step_infos,
        )
        return outputs, "Ran every screening stage and published the outputs."

    def _publish(self, merged, processed, flanked) -> None:
        """Rename the stage outputs to their published ``PREFIX.EEs.*`` names."""
        paths = self.paths
        renames = [
            (merged.ee_elements, paths.final_fasta),
            (merged.ee_elements_tax, paths.final_tax),
            (processed.ee_elements_gff3, paths.final_gff3),
            (flanked.flanks, paths.final_flanks),
        ]
        if self.clean_masked:
            renames += [
                (processed.ee_elements_cleaned, paths.final_cleaned_fasta),
                (processed.ee_elements_cleaned_tax, paths.final_cleaned_tax),
                (processed.ee_elements_cleaned_gff3, paths.final_cleaned_gff3),
            ]
        for produced, published in renames:
            self._rename(produced, published)

    def _archive(self) -> None:
        """Archive or delete every file in the output directory that is not published.

        Defined as the complement of the published outputs, so the
        translation methods' side files are covered without naming them.
        """
        paths = self.paths
        published = set(paths.final_outputs())
        leftovers = [
            entry
            for entry in sorted(Path(self.outdir).iterdir())
            if entry.is_file() and str(entry) not in published
        ]

        if self.removetmp:
            logger.warning("Removing temporary files.\n")
            for entry in leftovers:
                logger.debug(f"Removing temporary file {entry}")
                entry.unlink()
            return

        os.makedirs(paths.tmp_dir, exist_ok=True)
        for entry in leftovers:
            shutil.move(str(entry), f"{paths.tmp_dir}/{entry.name}")
        logger.info(
            f"Temporary files were moved to {paths.tmp_dir}. Check the tool "
            "documentation to access the description of each temporary file.\n"
        )
