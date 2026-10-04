"""Shared scaffolding for the ``screening`` stages.

A stage takes explicit input paths, runs on :meth:`Stage.run`, and returns a
dataclass naming every file it produced. The step classes in the package root
run on ``__init__`` and return nothing; a stage returns its outputs so it can be
composed.
"""

from __future__ import annotations
import os
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
from pathlib import Path
from eefinder.models import validate_all
from eefinder.stages.paths import ScreeningPaths
from eefinder.utils import StepInfo, check_outdir


@dataclass
class StageOutputs:
    """Base class for a stage's result.

    Subclasses add one ``str`` field per produced file.

    Attributes
    ----------
    step_info : StepInfo
        Timing and description record, appended to the run log by the caller.
    """

    step_info: StepInfo

    def files(self) -> "dict[str, str]":
        """Return the produced files, keyed by output name."""
        return {
            field.name: getattr(self, field.name)
            for field in fields(self)
            if field.name != "step_info" and isinstance(getattr(self, field.name), str)
        }

    def existing(self) -> "dict[str, str]":
        """Return the produced files that are actually on disk."""
        return {
            name: path for name, path in self.files().items() if Path(path).exists()
        }


class Stage(ABC):
    """One stage of the screening pipeline.

    Subclasses set the class attributes below and implement :meth:`_execute`.

    Attributes
    ----------
    name : str
        CLI subcommand name (``eefinder screening <name>``).
    title : str
        Human-readable stage name, as it appears in the run log.
    stage_id : str
        Zero-padded position in execution order, for ordering and display.
    """

    name: str = ""
    title: str = ""
    stage_id: str = ""

    #: Output models, as ``{outputs attribute: model}``, checked by
    #: :meth:`run`. A ``None`` attribute is skipped.
    models: dict = {}

    def __init__(self, outdir: str, prefix: str) -> None:
        # Created here, not in the CLI, so a stage can run on its own.
        self.outdir = check_outdir(outdir)
        self.prefix = prefix
        self.paths = ScreeningPaths(outdir=self.outdir, prefix=prefix)

    @abstractmethod
    def _execute(self) -> "tuple[StageOutputs, str]":
        """Do the work; return the outputs and the run-log message.

        The outputs' ``step_info`` is overwritten by :meth:`run`, so
        implementations may pass any placeholder.
        """

    def run(self) -> StageOutputs:
        """Execute the stage and return its outputs.

        Returns
        -------
        StageOutputs
            A stage-specific subclass naming every file produced, with
            :attr:`StageOutputs.step_info` filled in.
        """
        start = time.time()
        outputs, message = self._execute()
        validate_all(
            [
                (model, getattr(outputs, attribute, None))
                for attribute, model in self.models.items()
            ]
        )
        outputs.step_info = StepInfo.from_times(
            step=self.title,
            start_time=start,
            end_time=time.time(),
            message=message,
        )
        return outputs

    # -- helpers for subclasses -------------------------------------------
    @staticmethod
    def _rename(produced: str, canonical: str) -> str:
        """Move a step class's output to its canonical name.

        A no-op when the two paths are the same.
        """
        if produced == canonical:
            return canonical
        if not Path(produced).exists():
            raise FileNotFoundError(
                f"expected {produced} to have been produced, but it is missing"
            )
        os.replace(produced, canonical)
        return canonical

    def _require(self, *paths: str) -> None:
        """Fail with a usable message on a missing stage input."""
        missing = [path for path in paths if not Path(path).exists()]
        if missing:
            raise FileNotFoundError(
                f"stage '{self.name}' is missing required input(s): "
                + ", ".join(missing)
            )
