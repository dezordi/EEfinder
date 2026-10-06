"""Pluggable ways of producing the hit table.

A search method takes a nucleotide query and a protein database and writes a
hit table satisfying :data:`~eefinder.models.HIT_TABLE`.
:meth:`SearchMethod.run` enforces that model, so the stages downstream read
the model rather than knowing which method ran.

To add one: subclass :class:`SearchMethod`, set ``name`` and ``requires``,
implement :meth:`SearchMethod.search`, and decorate with
:func:`register_search_method`.
"""

from __future__ import annotations
import os
import shutil
from abc import ABC, abstractmethod
from typing import Dict, Type
from eefinder.models import HIT_TABLE, ModelError
from eefinder.log import logger
from eefinder.similarity_analysis import SimilaritySearch
from eefinder.translation import TRANSLATION_METHODS

SEARCH_METHODS: "Dict[str, Type[SearchMethod]]" = {}


def register_search_method(cls: "Type[SearchMethod]") -> "Type[SearchMethod]":
    """Class decorator that adds a search method to :data:`SEARCH_METHODS`."""
    if not cls.name:
        raise ValueError(f"{cls.__name__} must set a non-empty 'name'")
    SEARCH_METHODS[cls.name] = cls
    return cls


class SearchMethod(ABC):
    """One way of turning (genome, protein database) into a hit table.

    Attributes
    ----------
    name : str
        Registry key.
    requires : tuple of str
        External binaries the method needs on ``PATH``.
    produces : TableModel
        Model the output must satisfy.
    """

    name: str = ""
    requires: "tuple[str, ...]" = ()
    produces = HIT_TABLE

    def check_available(self) -> None:
        """Raise unless every required binary is on ``PATH``."""
        missing = [tool for tool in self.requires if shutil.which(tool) is None]
        if missing:
            raise FileNotFoundError(
                f"search method '{self.name}' needs {', '.join(missing)} on PATH"
            )

    @abstractmethod
    def search(self, query: str, database: str, threads: int, out_table: str) -> None:
        """Write the hit table for ``query`` vs ``database`` to ``out_table``.

        Implementations must emit nucleotide query coordinates on ``query``.
        """

    def run(self, query: str, database: str, threads: int, out_table: str) -> str:
        """Search, then enforce the output model.

        Returns
        -------
        str
            ``out_table``.

        Raises
        ------
        ModelError
            If the method's output does not satisfy :attr:`produces`.
        """
        self.check_available()
        logger.debug(
            f"search method '{self.name}': {query} vs {database} -> {out_table}"
        )
        self.search(query, database, threads, out_table)
        try:
            self.produces.validate(out_table)
        except ModelError as err:
            raise ModelError(
                f"search method '{self.name}' violated its output model: {err}"
            ) from err
        return out_table

    def describe(self) -> str:
        """One-line description for logs and ``--help``."""
        return f"{self.name} (needs {', '.join(self.requires) or 'nothing'})"


class _SimilaritySearchMethod(SearchMethod):
    """Base for the methods backed by :class:`SimilaritySearch`.

    :class:`SimilaritySearch` always writes ``{query}.blastx``; this moves it to
    the requested path.
    """

    mode: str = "blastx"
    translation_method: str = "default"

    def search(self, query: str, database: str, threads: int, out_table: str) -> None:
        SimilaritySearch(query, database, threads, self.mode, self.translation_method)
        produced = f"{query}.blastx"
        if produced != out_table:
            os.replace(produced, out_table)


@register_search_method
class BlastxSearch(_SimilaritySearchMethod):
    """Six-frame translated search with NCBI ``blastx`` -- the default."""

    name = "blastx"
    requires = ("blastx",)
    mode = "blastx"
    translation_method = "default"


@register_search_method
class DiamondSearch(_SimilaritySearchMethod):
    """Six-frame translated search with ``diamond blastx``.

    Parameters
    ----------
    sensitivity : str
        A DIAMOND sensitivity mode (``fast`` ... ``ultra-sensitive``).
    """

    name = "diamond"
    requires = ("diamond",)
    translation_method = "default"

    def __init__(self, sensitivity: str = "very-sensitive") -> None:
        self.mode = sensitivity

    def describe(self) -> str:
        return f"{self.name} --{self.mode} (needs diamond)"


@register_search_method
class PredictedProteinSearch(_SimilaritySearchMethod):
    """Align predicted proteins, then trace the coordinates back to the contig.

    Parameters
    ----------
    translation_method : str
        ``"gv"``, ``"rv"`` or ``"gv-rv"``.
    aligner : str
        ``"blastx"`` for ``blastp``, or a DIAMOND sensitivity for
        ``diamond blastp``.
    """

    name = "predicted"
    requires = ()

    def __init__(self, translation_method: str = "gv", aligner: str = "blastx") -> None:
        if translation_method not in TRANSLATION_METHODS:
            raise ValueError(
                f"unknown translation method {translation_method!r}; "
                f"expected one of {', '.join(TRANSLATION_METHODS)}"
            )
        self.translation_method = translation_method
        self.mode = aligner
        self.requires = ("blastp",) if aligner == "blastx" else ("diamond",)
        if translation_method == "gv-rv":
            self.requires = self.requires + ("cd-hit",)

    def describe(self) -> str:
        return (
            f"{self.name}:{self.translation_method} via {self.mode} "
            f"(needs {', '.join(self.requires)})"
        )


def resolve_search_method(mode: str, translation_method: str) -> SearchMethod:
    """Pick the search method for a ``--mode``/``--translation_method`` pair.

    Parameters
    ----------
    mode : str
        ``"blastx"`` or a DIAMOND sensitivity.
    translation_method : str
        ``"default"``, ``"gv"``, ``"rv"`` or ``"gv-rv"``.

    Returns
    -------
    SearchMethod
        A configured instance.

    Example
    -------
    >>> resolve_search_method("blastx", "default").name
    'blastx'
    >>> resolve_search_method("very-sensitive", "default").name
    'diamond'
    >>> resolve_search_method("blastx", "gv-rv").name
    'predicted'
    """
    if translation_method != "default":
        return PredictedProteinSearch(
            translation_method=translation_method, aligner=mode
        )
    if mode == "blastx":
        return BlastxSearch()
    return DiamondSearch(sensitivity=mode)
