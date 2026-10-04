"""Output models for the screening stages.

A model is a tabular or FASTA schema plus a check. Each stage declares the
models its outputs satisfy, and :meth:`eefinder.stages.base.Stage.run` enforces
them, so a stage states its guarantees instead of documenting them.

:data:`HIT_TABLE` is the model every search method must satisfy, whatever
produces the hits.
"""

from __future__ import annotations
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from eefinder.filter_table import FILTERED_COLUMNS, OUTFMT6_COLUMNS
from eefinder.get_taxonomy import TAXONOMY_COLUMNS
from eefinder.log import logger


class ModelError(ValueError):
    """A stage output does not satisfy the model it declared."""


@dataclass(frozen=True)
class TableModel:
    """The shape of a tabular file passed between stages.

    Parameters
    ----------
    name : str
        Human-readable name, used in error messages.
    columns : tuple of str
        Required columns, in order.
    sep : str
        Field separator.
    header : bool
        Whether the file carries a header row. Headerless files are checked by
        field count instead.
    allow_empty : bool
        Whether zero data rows is acceptable.
    extra_columns : bool
        Whether columns beyond :attr:`columns` are tolerated.
    """

    name: str
    columns: "tuple[str, ...]"
    sep: str = "\t"
    header: bool = True
    allow_empty: bool = True
    extra_columns: bool = True

    def validate(self, path: str) -> None:
        """Check ``path`` against this model.

        Parameters
        ----------
        path : str
            File to check.

        Raises
        ------
        ModelError
            If the file is missing, empty when that is not allowed, or does not
            carry the required columns.
        """
        file = Path(path)
        if not file.is_file():
            raise ModelError(f"{self.name}: {path} was not produced")

        with open(path, newline="") as handle:
            rows = csv.reader(handle, delimiter=self.sep)
            try:
                first = next(rows)
            except StopIteration:
                if self.allow_empty:
                    logger.debug(f"{self.name}: {path} is empty (allowed)")
                    return
                raise ModelError(f"{self.name}: {path} is empty")

            if self.header:
                missing = [name for name in self.columns if name not in first]
                if missing:
                    raise ModelError(
                        f"{self.name}: {path} is missing column(s) "
                        f"{', '.join(missing)}. Expected: "
                        f"{', '.join(self.columns)}"
                    )
                if not self.extra_columns and len(first) != len(self.columns):
                    raise ModelError(
                        f"{self.name}: {path} has {len(first)} column(s), "
                        f"expected exactly {len(self.columns)}"
                    )
            else:
                expected = len(self.columns)
                if len(first) < expected:
                    raise ModelError(
                        f"{self.name}: {path} has {len(first)} field(s) in its "
                        f"first row, expected at least {expected} "
                        f"({', '.join(self.columns)})"
                    )
                if not self.extra_columns and len(first) != expected:
                    raise ModelError(
                        f"{self.name}: {path} has {len(first)} field(s), "
                        f"expected exactly {expected}"
                    )
        logger.debug(f"{self.name}: {path} satisfies the model")


@dataclass(frozen=True)
class FastaModel:
    """The shape of a FASTA passed between stages.

    Parameters
    ----------
    name : str
        Human-readable name, used in error messages.
    allow_empty : bool
        Whether zero records is acceptable.
    """

    name: str
    allow_empty: bool = True

    def validate(self, path: str) -> None:
        """Check that ``path`` exists and is a FASTA.

        Raises
        ------
        ModelError
            If the file is missing, or its first content line is not a header.
        """
        file = Path(path)
        if not file.is_file():
            raise ModelError(f"{self.name}: {path} was not produced")
        with open(path) as handle:
            for line in handle:
                if line.strip():
                    if not line.startswith(">"):
                        raise ModelError(
                            f"{self.name}: {path} does not start with a FASTA " "header"
                        )
                    return
        if not self.allow_empty:
            raise ModelError(f"{self.name}: {path} has no records")
        logger.debug(f"{self.name}: {path} is empty (allowed)")


#: Hit table every search method must emit: BLAST ``outfmt 6``, tab-separated,
#: no header, query coordinates in nucleotides on the cleaned genome.
HIT_TABLE = TableModel(
    name="hit table (outfmt6)",
    columns=tuple(OUTFMT6_COLUMNS),
    header=False,
)

#: Hits after redundant-alignment collapsing: outfmt6 plus the EEfinder
#: annotations (``sense``, ``bed_name``, ``tag``).
FILTERED_HIT_TABLE = TableModel(
    name="filtered hit table",
    columns=tuple(FILTERED_COLUMNS),
    header=True,
)

#: Candidate element intervals: headerless BED3 on the cleaned genome.
CANDIDATE_BED = TableModel(
    name="candidate BED",
    columns=("contig", "start", "end"),
    header=False,
)

#: Per-hit taxonomy signature: the filtered hits joined to the ``-mt`` metadata.
#: Comma-separated, because it is written by pandas as a CSV.
TAXONOMY_SIGNATURE = TableModel(
    name="taxonomy signature",
    columns=("sseqid", "Species", "Genus", "Family", "Molecule_type", "Host"),
    sep=",",
    header=True,
)

#: Per-element taxonomy table -- what the user reads as ``PREFIX.EEs.tax.tsv``.
ELEMENT_TABLE = TableModel(
    name="element taxonomy table",
    columns=tuple(TAXONOMY_COLUMNS),
    header=True,
)

#: Candidate element sequences cut out of the cleaned genome.
CANDIDATE_FASTA = FastaModel(name="candidate sequences")

#: Merged element sequences.
ELEMENT_FASTA = FastaModel(name="element sequences")

#: The cleaned genome: the query of every search in a run.
GENOME_FASTA = FastaModel(name="cleaned genome", allow_empty=False)


def validate_all(pairs: "list[tuple[object, Optional[str]]]") -> None:
    """Validate several (model, path) pairs, skipping ``None`` paths.

    Parameters
    ----------
    pairs : list of (model, path or None)
        What to check.

    Raises
    ------
    ModelError
        On the first violation.
    """
    for model, path in pairs:
        if path is not None:
            model.validate(path)
