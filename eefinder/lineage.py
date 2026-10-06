"""T
EEfinder stores the whole lineage, realm to species, as
``r__{realm};k__{kingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species}``
-- one field per rank, always all eight, with :data:`UNKNOWN` where NCBI has no
taxon at that rank.

Family- and genus-based logic reads the rank it needs out of the string. When
that rank is :data:`UNKNOWN` the whole lineage is used instead, so two
unclassified viruses from different orders are not treated as the same taxon.
"""

from __future__ import annotations

RANKS = (
    ("realm", "r"),
    ("kingdom", "k"),
    ("phylum", "p"),
    ("class", "c"),
    ("order", "o"),
    ("family", "f"),
    ("genus", "g"),
    ("species", "s"),
)
RANK_BY_PREFIX = {prefix: rank for rank, prefix in RANKS}
SEPARATOR = ";"
UNKNOWN = "Unk"


def format_lineage(ranks: "dict[str, str]") -> str:
    """Build a ``Taxonomy`` string from a rank -> name mapping.

    Parameters
    ----------
    ranks : dict[str, str]
        Rank names (lower-cased, as NCBI reports them) to taxon names. A rank
        that is absent, empty or literally ``"Unknown"`` becomes
        :data:`UNKNOWN`.

    Returns
    -------
    str
        All eight fields, in :data:`RANKS` order.

    Example
    -------
    >>> format_lineage({"realm": "Riboviria", "family": "Famviridae"})
    'r__Riboviria;k__Unk;p__Unk;c__Unk;o__Unk;f__Famviridae;g__Unk;s__Unk'
    """
    fields = []
    for rank, prefix in RANKS:
        name = (ranks.get(rank) or "").strip()
        if not name or name.lower() in ("unknown", UNKNOWN.lower()):
            name = UNKNOWN
        fields.append(f"{prefix}__{name}")
    return SEPARATOR.join(fields)


def parse_lineage(taxonomy: str) -> "dict[str, str]":
    """Read a ``Taxonomy`` string back into a rank -> name mapping.

    Unknown prefixes and malformed fields are ignored, so a hand-edited or
    future-extended string does not raise.

    Parameters
    ----------
    taxonomy : str
        A ``Taxonomy`` value.

    Returns
    -------
    dict[str, str]
        Only the ranks that carry a real name; :data:`UNKNOWN` counts as
        absent, so the fallbacks below treat it as no taxon at all.

    Example
    -------
    >>> parse_lineage("r__Riboviria;f__Famviridae;g__Unk;s__Unk")
    {'realm': 'Riboviria', 'family': 'Famviridae'}
    """
    found: "dict[str, str]" = {}
    for field in str(taxonomy or "").split(SEPARATOR):
        prefix, sep, name = field.strip().partition("__")
        if not sep:
            continue
        rank = RANK_BY_PREFIX.get(prefix)
        name = name.strip()
        if rank and name and name.lower() not in ("unknown", UNKNOWN.lower()):
            found[rank] = name
    return found


def rank_name(taxonomy: str, rank: str) -> str:
    """Return the taxon at ``rank``, or an empty string when absent.

    Example
    -------
    >>> rank_name("r__Riboviria;f__Famviridae", "family")
    'Famviridae'
    >>> rank_name("r__Riboviria;f__Unk", "family")
    ''
    """
    return parse_lineage(taxonomy).get(rank, "")


def merge_key(taxonomy: str, level: str) -> str:
    """Return what identifies a taxon for merging at ``level``.

    The taxon name at ``level`` when it has one, otherwise the lineage itself,
    so unranked records only share a key with others of the same lineage.

    Parameters
    ----------
    taxonomy : str
        A ``Taxonomy`` value.
    level : str
        ``"family"`` or ``"genus"``.

    Returns
    -------
    str

    Example
    -------
    >>> merge_key("r__Riboviria;f__Famviridae;g__Unk", "family")
    'Famviridae'
    >>> merge_key("r__Duplodnaviria;c__Caudoviricetes;f__Unk", "family")
    'r__Duplodnaviria;c__Caudoviricetes;f__Unk'
    """
    return parse_lineage(taxonomy).get(level) or str(taxonomy or "")


def trim_lineage(taxonomy: str) -> str:
    """Drop the trailing :data:`UNKNOWN` ranks of a lineage.

    A lineage resolved only down to the family is written
    ``...;f__Famviridae`` rather than padded with ``g__Unk;s__Unk``.

    Example
    -------
    >>> trim_lineage("r__Riboviria;f__Famviridae;g__Unk;s__Unk")
    'r__Riboviria;f__Famviridae'
    """
    fields = str(taxonomy or "").split(SEPARATOR)
    while fields and fields[-1].partition("__")[2].strip() in ("", UNKNOWN):
        fields.pop()
    return SEPARATOR.join(fields)


def lca(taxonomies: "list[str]") -> str:
    """Return the lowest common lineage of several ``Taxonomy`` values.

    A rank is kept while every lineage names the same taxon there; from the
    first rank they disagree on down, every rank is :data:`UNKNOWN`. Merging
    elements of two species of one genus therefore yields the genus, and of two
    genera the family. A rank *none* of them names is not a disagreement, so a
    shared family survives a lineage with no realm.

    Parameters
    ----------
    taxonomies : list[str]
        ``Taxonomy`` values; empty ones are ignored.

    Returns
    -------
    str
        The common lineage, trimmed of its unresolved tail
        (:func:`trim_lineage`), or an empty string when nothing was given.

    Example
    -------
    >>> lca(["r__Riboviria;f__Fam;g__Gen;s__A", "r__Riboviria;f__Fam;g__Gen;s__B"])
    'r__Riboviria;k__Unk;p__Unk;c__Unk;o__Unk;f__Fam;g__Gen'
    """
    parsed = [parse_lineage(t) for t in taxonomies if str(t or "").strip()]
    if not parsed:
        return ""
    common: "dict[str, str]" = {}
    for rank, _ in RANKS:
        names = {entry.get(rank, "") for entry in parsed}
        if len(names) != 1:
            break
        name = names.pop()
        if name:
            common[rank] = name
    return trim_lineage(format_lineage(common))


def is_ranked(taxonomy: str, level: str) -> bool:
    """Whether ``taxonomy`` carries a taxon at ``level``."""
    return bool(parse_lineage(taxonomy).get(level))


def lineage_from_nodes(
    lineage: "list[dict]",
    ranks: "dict[int, str]",
    species: str = "",
) -> str:
    """Build a ``Taxonomy`` string from a report lineage and a rank lookup.

    The ``data_report.jsonl`` lineage is a list of ``{"name", "taxId"}`` with no
    ranks, so the ranks come from a separate lookup keyed by tax id.

    Parameters
    ----------
    lineage : list[dict]
        Lineage entries, each with ``name`` and ``taxId``.
    ranks : dict[int, str]
        Tax id -> rank name.
    species : str
        Organism name, used for ``s__`` when the lineage has no species entry.

    Returns
    -------
    str
    """
    found: "dict[str, str]" = {}
    for entry in lineage or ():
        tax_id = entry.get("taxId")
        rank = ranks.get(int(tax_id)) if tax_id is not None else None
        if rank in RANK_BY_PREFIX.values() or rank in dict(RANKS):
            found.setdefault(rank, entry.get("name", ""))
    if species:
        found["species"] = species
    return format_lineage(found)
