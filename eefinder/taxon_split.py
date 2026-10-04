"""Split a download into one request per taxon at a given rank.

A single request for a broad taxon (every virus, every bacterium) is the
transfer that fails and hangs. Asking for one family -- or one genus -- at a
time keeps each request small, and the packages are concatenated before the
database is built.

Records whose lineage has no taxon at the chosen rank cannot be reached this
way. They are reported as :class:`SkippedTaxon` rather than downloaded: without
a family there is no ``Family``/``Genus`` to assign and no ``Molecule_type`` to
look up, so the elements they would produce carry no taxonomy.
"""

from __future__ import annotations
import json
import shlex
import subprocess
from typing import NamedTuple
from eefinder.log import logger
from eefinder.taxon_exclusion import (
    DATASETS_BINARY,
    TaxonNode,
    _Taxonomy,
    resolve_taxon,
    summarize_taxa,
)

#: Ranks a download may be split at.
SPLIT_LEVELS = ("family", "genus")

#: Value of ``--split-level`` that downloads the requested taxon in one request.
NO_SPLIT = "none"


class SkippedTaxon(NamedTuple):
    """A subtree left out because it has no taxon at the split rank."""

    tax_id: int
    name: str
    rank: str
    assembly_count: int
    reason: str


class SplitPlan(NamedTuple):
    """The taxa to download, and what was left out."""

    #: Rank the split was made at, or :data:`NO_SPLIT`.
    level: str
    root: TaxonNode
    #: One request per entry.
    taxa: "tuple[TaxonNode, ...]"
    #: Subtrees with no taxon at ``level``, reported instead of downloaded.
    skipped: "tuple[SkippedTaxon, ...]"

    @property
    def split(self) -> bool:
        """Whether the plan is more than a single request for the root."""
        return len(self.taxa) > 1

    @property
    def skipped_assemblies(self) -> int:
        """Assemblies the skipped subtrees hold, as NCBI counts them."""
        return sum(entry.assembly_count for entry in self.skipped)


def taxa_at_rank(
    root: str, rank: str, datasets_bin: str = DATASETS_BINARY
) -> "dict[int, TaxonNode]":
    """Return every taxon of ``rank`` below ``root``, keyed by tax id.

    One ``datasets`` call: ``--rank`` makes the CLI walk the subtree itself.

    Parameters
    ----------
    root : str
        Tax id or scientific name of the taxon being downloaded.
    rank : str
        ``"family"`` or ``"genus"``.
    datasets_bin : str
        Path/name of the ``datasets`` executable.

    Raises
    ------
    RuntimeError
        If the ``datasets`` call fails.
    """
    command = (
        f"{datasets_bin} summary taxonomy taxon {shlex.quote(str(root))} "
        f"--rank {shlex.quote(rank)} --as-json-lines"
    )
    logger.debug(f"taxonomy command: {command}")
    result = subprocess.run(
        shlex.split(command),
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"could not list {rank} taxa below {root!r}: "
            f"{(result.stderr or '').strip() or 'datasets exited ' + str(result.returncode)}"
        )
    nodes: "dict[int, TaxonNode]" = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            continue
        record = payload.get("taxonomy", payload)
        tax_id = record.get("tax_id")
        if tax_id is None or str(record.get("rank") or "").lower() != rank:
            continue
        counts = {
            entry.get("type"): entry.get("count", 0)
            for entry in (record.get("counts") or [])
        }
        nodes[int(tax_id)] = TaxonNode(
            tax_id=int(tax_id),
            name=(record.get("current_scientific_name") or {}).get("name", ""),
            parents=tuple(record.get("parents") or ()),
            children=tuple(record.get("children") or ()),
            rank=rank,
            assembly_count=int(counts.get("COUNT_TYPE_ASSEMBLY") or 0),
        )
    return nodes


def _uncovered_roots(
    taxonomy: "_Taxonomy",
    root_id: int,
    covered: "dict[int, TaxonNode]",
) -> "list[int]":
    """Highest subtrees below ``root_id`` that contain none of ``covered``.

    Every ancestor of a covered taxon has at least one, so walking down from the
    root and stopping at the first node that is neither covered nor such an
    ancestor yields the minimal set of uncovered subtrees. Children are looked
    up one level at a time, in batches, to keep the number of calls small.
    """
    bearing = {parent for node in covered.values() for parent in node.parents}
    uncovered: "list[int]" = []
    level = [root_id]
    while level:
        taxonomy.prime(summarize_taxa([str(t) for t in level], taxonomy.datasets_bin))
        next_level: "list[int]" = []
        for tax_id in level:
            if tax_id in covered:
                continue
            if tax_id in bearing:
                next_level.extend(taxonomy.node(tax_id).children)
            else:
                uncovered.append(tax_id)
        level = next_level
    return uncovered


def plan_split(
    root: str,
    level: str,
    datasets_bin: str = DATASETS_BINARY,
) -> SplitPlan:
    """Work out one download per taxon of ``level`` below ``root``.

    When ``root`` is itself at or below ``level`` -- a single family, say --
    there is nothing to split and the plan is one request for it.

    Parameters
    ----------
    root : str
        Tax id or scientific name of the taxon being downloaded.
    level : str
        ``"family"``, ``"genus"`` or :data:`NO_SPLIT`.
    datasets_bin : str
        Path/name of the ``datasets`` executable.

    Returns
    -------
    SplitPlan

    Raises
    ------
    ValueError
        If ``level`` is not a known split level.
    RuntimeError
        If ``root`` cannot be resolved.
    """
    if level != NO_SPLIT and level not in SPLIT_LEVELS:
        raise ValueError(
            f"unknown split level {level!r}; expected one of "
            f"{', '.join(SPLIT_LEVELS + (NO_SPLIT,))}"
        )

    node = resolve_taxon(root, datasets_bin)
    if node is None:
        raise RuntimeError(f"NCBI taxonomy has no taxon {root!r}")

    if level == NO_SPLIT:
        return SplitPlan(level=NO_SPLIT, root=node, taxa=(node,), skipped=())

    if node.rank == level:
        logger.info(f"{node.name} is already a {level}; downloading it in one request")
        return SplitPlan(level=level, root=node, taxa=(node,), skipped=())

    covered = taxa_at_rank(root, level, datasets_bin)
    if not covered:
        logger.warning(
            f"no {level} taxa were found below {node.name} ({node.tax_id}); "
            "downloading it in one request"
        )
        return SplitPlan(level=level, root=node, taxa=(node,), skipped=())

    logger.info(
        f"{node.name} covers {len(covered)} {level} taxa; downloading one at a time"
    )
    taxonomy = _Taxonomy(datasets_bin)
    taxonomy.prime({node.tax_id: node})
    taxonomy.prime(covered)
    uncovered = _uncovered_roots(taxonomy, node.tax_id, covered)

    skipped = tuple(
        SkippedTaxon(
            tax_id=entry.tax_id,
            name=entry.name,
            rank=entry.rank,
            assembly_count=entry.assembly_count,
            reason=f"no {level} in its lineage",
        )
        for entry in (taxonomy.node(tax_id) for tax_id in uncovered)
    )
    taxa = tuple(covered[tax_id] for tax_id in sorted(covered))
    return SplitPlan(level=level, root=node, taxa=taxa, skipped=skipped)
