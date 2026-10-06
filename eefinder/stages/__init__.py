"""The screening logic as composable stages.

===============  ===========================================================
Subcommand       What it does
===============  ===========================================================
``prepare``      Validate ``-mt``/``-db``/``-bt``, build the search indexes
``clean``        Drop short contigs, prefix the headers
``align``        Translated search, collapse redundant hits
``filter``       Reverse search against the host-gene baits
``taxonomy``     Join the surviving hits to the metadata
``merge``        Merge fragments, build the element table
``postprocess``  Repeat filter, overlap resolution, GFF3
``flanks``       Extract the flanking regions
``all``          Every stage, in order
===============  ===========================================================

Each is a class with the same shape: configure it, call ``run()``, read the
returned outputs.
"""

from eefinder.stages.prepare import PrepareInputs, PrepareInputsOutputs
from eefinder.stages.cleaning import DataCleaning, DataCleaningOutputs
from eefinder.stages.alignment import SequenceAlignment, SequenceAlignmentOutputs
from eefinder.stages.filtering import (
    PutativeElementsFilter,
    PutativeElementsFilterOutputs,
)
from eefinder.stages.taxonomy import TaxonomyAssignment, TaxonomyAssignmentOutputs
from eefinder.stages.merging import (
    MergeFragmentedElements,
    MergeFragmentedElementsOutputs,
)
from eefinder.stages.postprocess import PostProcessing, PostProcessingOutputs
from eefinder.stages.flanks import FlanksExtraction, FlanksExtractionOutputs
from eefinder.stages.pipeline import ScreeningPipeline, ScreeningPipelineOutputs
from eefinder.stages.base import Stage, StageOutputs
from eefinder.stages.paths import ScreeningPaths

STAGES = (
    PrepareInputs,
    DataCleaning,
    SequenceAlignment,
    PutativeElementsFilter,
    TaxonomyAssignment,
    MergeFragmentedElements,
    PostProcessing,
    FlanksExtraction,
    ScreeningPipeline,
)

STAGES_BY_NAME = {stage.name: stage for stage in STAGES}

__all__ = [
    "Stage",
    "StageOutputs",
    "ScreeningPaths",
    "STAGES",
    "STAGES_BY_NAME",
    "PrepareInputs",
    "PrepareInputsOutputs",
    "DataCleaning",
    "DataCleaningOutputs",
    "SequenceAlignment",
    "SequenceAlignmentOutputs",
    "PutativeElementsFilter",
    "PutativeElementsFilterOutputs",
    "TaxonomyAssignment",
    "TaxonomyAssignmentOutputs",
    "MergeFragmentedElements",
    "MergeFragmentedElementsOutputs",
    "PostProcessing",
    "PostProcessingOutputs",
    "FlanksExtraction",
    "FlanksExtractionOutputs",
    "ScreeningPipeline",
    "ScreeningPipelineOutputs",
]
