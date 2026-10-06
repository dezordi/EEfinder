"""Protein-name standardisation used by ``get-databases``.

:func:`standardize_protein` is the entry point. It runs a shared cleaning
pipeline and then, for the ``virus`` target, the rules in
``data/protein_rules.yaml``; ``bacteria`` and ``host`` are cleaned only.
``data/README.md`` documents the rule file and its scopes.
"""

from __future__ import annotations
import re
from pathlib import Path
from typing import Callable, NamedTuple, Optional
import yaml
from eefinder.lineage import parse_lineage
from eefinder.log import logger

_RULE_FILE = Path(__file__).resolve().parent / "data" / "protein_rules.yaml"


def _load_rules() -> dict:
    """Read the bundled rule file, empty when it is missing or unreadable."""
    if not _RULE_FILE.is_file():
        logger.warning(f"no protein rule file at {_RULE_FILE}; names are only cleaned")
        return {}
    try:
        with open(_RULE_FILE) as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as err:
        logger.warning(
            f"could not read {_RULE_FILE.name} ({err}); names are only cleaned"
        )
        return {}
    return loaded if isinstance(loaded, dict) else {}


_RULES = _load_rules()


_MW_RE = re.compile(
    r"\b(\d+)[\s-]*(?:kda|kd|k|l)(?:[- ]*(?:putative\s+)?(?:nonstructural\s+)?(?:protein|like\s+protein))?\b",
    re.IGNORECASE,
)


def normalize_molecular_weight(name: str) -> str:
    """Rewrite a molecular-weight token to ``"X kDa protein"``.

    Examples
    --------
    >>> normalize_molecular_weight("100 kDa")
    '100 kDa protein'
    >>> normalize_molecular_weight("33-kDa")
    '33 kDa protein'
    >>> normalize_molecular_weight("33K-like protein")
    '33 kDa protein'
    >>> normalize_molecular_weight("33KD putative nonstructural protein")
    '33 kDa protein'
    """
    return _MW_RE.sub(r"\1 kDa protein", name)


_BRACKET_TAG_RE = re.compile(r"\[[A-Za-z_]+=(?:[^\[\]]|\[[^\[\]]*\])*\]")


def strip_bracket_tags(name: str) -> str:
    """Remove NCBI ``[key=value]`` metadata tags, leaving whitespace uncollapsed.

    Examples
    --------
    >>> strip_bracket_tags("nucleoprotein [organism=Some virus]").strip()
    'nucleoprotein'
    """
    return _BRACKET_TAG_RE.sub(" ", name)


_SPECIAL_CHARS = ":,/\\?!\"'"

_QUALIFIER_RE = re.compile(r"^(putative|predicted|probable|hypothetical)\s+")
_TRAILING_RE = re.compile(r"(,\s*partial|\s+precursor)$")

_OUTPUT_QUALIFIER_RE = re.compile(
    r"^\s*(?:putative|putatively|predicted|probable|possible|presumed|presumptive)\s+",
    re.IGNORECASE,
)


def _strip_qualifiers(text: str) -> str:
    """Drop leading hedging qualifiers (e.g. "putative ") from a protein name."""
    prev = None
    while prev != text:
        prev = text
        text = _OUTPUT_QUALIFIER_RE.sub("", text)
    return text


_NS_DESIGNATION_RE = re.compile(
    r"^(NS\d+[A-Za-z]?)(?:[\s-]+(?:like[\s-]+)?protein|[\s-]+peptide)$",
    re.IGNORECASE,
)


def _strip_designation_suffix(text: str) -> str:
    """Drop a redundant "protein"/"peptide" suffix from an ``NSxx`` designation."""
    match = _NS_DESIGNATION_RE.match(text.strip())
    return match.group(1).upper() if match else text


_LEADING_DIRECTIVE_RE = re.compile(r"^\s*(cds|orf)\s*:\s*", re.IGNORECASE)

_UNKNOWN_TOKENS = {"", "cds", "orf"}

_TYPO_CORRECTIONS: dict[str, str] = {
    str(wrong).strip().lower(): str(right).strip()
    for right, wrongs in (_RULES.get("typos") or {}).items()
    for wrong in (wrongs or ())
    if str(wrong).strip()
}

_TYPO_RE = re.compile(
    r"\b(?:%s)\b"
    % "|".join(re.escape(t) for t in sorted(_TYPO_CORRECTIONS, key=len, reverse=True)),
    re.IGNORECASE,
)

_REWRITES: "list[tuple[re.Pattern, str]]" = [
    (re.compile(entry["pattern"], re.IGNORECASE), entry.get("replacement") or "")
    for entry in (_RULES.get("rewrites") or ())
    if isinstance(entry, dict) and (entry.get("pattern") or "").strip()
]


def _apply_rewrites(text: str) -> str:
    """Apply the bundled shape-fix rewrites, in file order."""
    for pattern, replacement in _REWRITES:
        text = pattern.sub(replacement, text)
    return text


def _apply_typos(text: str) -> str:
    """Correct known misspellings (whole-word, case-insensitive)."""
    return _TYPO_RE.sub(lambda m: _TYPO_CORRECTIONS[m.group(0).lower()], text)


def _normalize_product(name: str) -> str:
    """Normalise a raw protein product into a match key for a protein map."""
    text = re.sub(r"[-_/\\();,]", " ", name.lower())
    text = re.sub(r"\s+", " ", text).strip()
    text = _QUALIFIER_RE.sub("", text)
    text = _TRAILING_RE.sub("", text)
    text = text.strip()
    text = _apply_typos(text)
    text = re.sub(r"\bpolymerases\b", "polymerase", text)
    return text


def _clean_and_capitalize(name: str) -> str:
    """Strip special characters and capitalise a lower-case leading letter."""
    for char in _SPECIAL_CHARS:
        name = name.replace(char, "")
    name = re.sub(r"\s+", " ", name).strip()
    if name and name[0].islower():
        name = name[0].upper() + name[1:]
    return name


ProteinMapper = Callable[[str, str, "frozenset[str]"], Optional[str]]


def _standardize(
    name: str,
    mol_type: str,
    mapper: Optional[ProteinMapper],
    taxonomy: str = "",
) -> str:
    """Run the shared cleaning pipeline, optionally applying a target ``mapper``.

    A ``mapper`` hit yields its canonical name, otherwise the cleaned name is
    kept. A name that is only a ``CDS``/``ORF`` directive, or that begins with
    ``hypothetical``, becomes ``"Unknown"``.

    Parameters
    ----------
    name : str
        Raw protein product (e.g. a FASTA-header product).
    mol_type : str
        The record's ``Molecule_type``, passed through to ``mapper`` for scoping.
    mapper : ProteinMapper, optional
        Target-specific canonicalisation; ``None`` for generic cleaning only.
    taxonomy : str
        The record's lineage, passed through to ``mapper`` for scoping.

    Returns
    -------
    str
        The standardised, cleaned and capitalised protein name.
    """
    stripped = strip_bracket_tags(name)
    stripped = _LEADING_DIRECTIVE_RE.sub("", stripped)
    stripped = normalize_molecular_weight(stripped)
    stripped = _apply_typos(stripped)
    stripped = _apply_rewrites(stripped)
    stripped = _strip_qualifiers(stripped)
    if stripped.strip().lower().startswith("hypothetical"):
        return "Unknown"
    stripped = _strip_designation_suffix(stripped)

    if mapper is not None:
        suggested = mapper(
            _normalize_product(stripped), mol_type, lineage_taxa(taxonomy)
        )
        if suggested is not None:
            return _clean_and_capitalize(suggested)

    result = _clean_and_capitalize(stripped)
    if result.lower() in _UNKNOWN_TOKENS:
        return "Unknown"
    return result


class _Rule(NamedTuple):
    """One rule of the ``proteins`` section, flattened out of the nesting."""

    suggested: str
    current: str
    match_type: str
    mol_scope: str
    taxon_scope: str
    pattern: Optional["re.Pattern"]

    def match_position(self, normalized: str) -> Optional[int]:
        """Where this rule's key matches, or ``None`` when it does not."""
        if self.match_type == "exact":
            return 0 if normalized == self.current else None
        if not self.pattern:
            return None
        found = self.pattern.search(normalized)
        return found.start() if found else None


def _compile_rule(
    suggested: str, current: str, match_type: str, mol_scope: str, taxon_scope: str
) -> Optional[_Rule]:
    """Build a :class:`_Rule`, or ``None`` for an unusable match type/pattern."""
    if match_type == "exact":
        pattern = None
    elif match_type == "contains":
        pattern = re.compile(rf"\b{re.escape(current)}\b")
    elif match_type == "regex":
        try:
            pattern = re.compile(current)
        except re.error:
            return None
    else:
        return None
    return _Rule(
        suggested=suggested,
        current=current,
        match_type=match_type,
        mol_scope=mol_scope,
        taxon_scope=taxon_scope,
        pattern=pattern,
    )


def _load_protein_map() -> "list[_Rule]":
    """Flatten the nested ``proteins`` section into an ordered list of rules."""
    rules: "list[_Rule]" = []
    for taxon_scope, by_mol in (_RULES.get("proteins") or {}).items():
        for mol_scope, by_name in (by_mol or {}).items():
            for suggested, by_type in (by_name or {}).items():
                for match_type, keys in (by_type or {}).items():
                    for key in keys or ():
                        rule = _compile_rule(
                            str(suggested).strip(),
                            str(key).strip(),
                            str(match_type).strip(),
                            str(mol_scope).strip(),
                            str(taxon_scope).strip(),
                        )
                        if rule is None:
                            logger.warning(
                                f"ignoring protein rule {suggested!r} <- {key!r}: "
                                f"unusable match type {match_type!r}"
                            )
                            continue
                        rules.append(rule)
    return rules


_PROTEIN_RULES = _load_protein_map()

_ANY_TAXON = "any"


def _in_scope(mol_type: str, scope: str) -> bool:
    """Whether ``mol_type`` satisfies a ``molecule_type_scope`` expression."""
    mol = mol_type or ""
    for token in (part.strip() for part in scope.split(";")):
        if token == "any":
            return True
        if token == "RT" and "RT" in mol:
            return True
        if token == "RNA" and "RNA" in mol and "RT" not in mol:
            return True
        if token == "+ssRNA" and mol.startswith("ssRNA(+)"):
            return True
        if token == "-ssRNA" and mol.startswith("ssRNA(-)"):
            return True
        if token == "dsRNA" and mol.startswith("dsRNA"):
            return True
        if token == "dsDNA" and mol.startswith("dsDNA"):
            return True
        if token == "ssDNA" and mol.startswith("ssDNA"):
            return True
    return False


def lineage_taxa(taxonomy: str) -> "frozenset[str]":
    """Return the lower-cased taxon names of a lineage, for ``taxon_scope``.

    Parameters
    ----------
    taxonomy : str
        A ``Taxonomy`` value, or a bare taxon name.

    Returns
    -------
    frozenset[str]
        Every named rank, or the bare name on its own.
    """
    text = str(taxonomy or "").strip()
    if not text:
        return frozenset()
    parsed = parse_lineage(text)
    if parsed:
        return frozenset(name.lower() for name in parsed.values())
    return frozenset({text.lower()})


def _taxon_in_scope(taxa: "frozenset[str]", scope: str) -> bool:
    """Whether any taxon of the record's lineage satisfies ``taxon_scope``."""
    return any(
        token.strip().lower() in taxa
        for token in scope.split(";")
        if token.strip() and token.strip() != _ANY_TAXON
    )


def _viral_mapper(
    normalized: str, mol_type: str, taxa: "frozenset[str]" = frozenset()
) -> Optional[str]:
    """Look ``normalized`` up in the viral map, honouring both scopes.

    Rules are ranked by: taxon-scoped before unscoped, ``exact`` before
    ``contains``/``regex``, earliest match in the name, then file order.
    """
    for taxon_scoped in (True, False):
        best: "Optional[tuple[int, int, str]]" = None
        for order, rule in enumerate(_PROTEIN_RULES):
            if (rule.taxon_scope != _ANY_TAXON) is not taxon_scoped:
                continue
            if taxon_scoped and not _taxon_in_scope(taxa, rule.taxon_scope):
                continue
            if not _in_scope(mol_type, rule.mol_scope):
                continue
            position = rule.match_position(normalized)
            if position is None:
                continue
            if rule.match_type == "exact":
                return rule.suggested
            if best is None or (position, order) < best[:2]:
                best = (position, order, rule.suggested)
        if best is not None:
            return best[2]
    return None


def _standardize_virus(name: str, mol_type: str = "", taxonomy: str = "") -> str:
    """Standardise a viral protein name via the bundled viral protein map."""
    return _standardize(name, mol_type, _viral_mapper, taxonomy)


def _standardize_bacteria(name: str, mol_type: str = "", taxonomy: str = "") -> str:
    """Standardise a bacterial protein name (generic cleaning only, for now).

    There is no bacterial protein map yet; this is the extension point for
    bacteria-specific canonicalisation.
    """
    return _standardize(name, mol_type, None)


def _standardize_host(name: str, mol_type: str = "", taxonomy: str = "") -> str:
    """Standardise a host protein name (generic cleaning only).

    Host baits are gene/protein names kept as-is aside from generic cleaning;
    this is the extension point for host-specific rules.
    """
    return _standardize(name, mol_type, None)


_STANDARDIZERS: dict[str, Callable[[str, str, str], str]] = {
    "virus": _standardize_virus,
    "bacteria": _standardize_bacteria,
    "host": _standardize_host,
}


def standardize_protein(
    name: str, mol_type: str = "", target: str = "virus", taxonomy: str = ""
) -> str:
    """Standardise a raw protein name for the given ``target`` database.

    Parameters
    ----------
    name : str
        Raw protein product (e.g. a FASTA-header product).
    mol_type : str
        The record's ``Molecule_type``, used to scope the map rules.
    target : str
        The database target — one of ``"virus"``, ``"bacteria"`` or ``"host"``;
        selects the target-specific standardisation logic.
    taxonomy : str
        The record's lineage, used to scope the taxon-specific rules. Without
        it only the unscoped rules apply.

    Returns
    -------
    str
        The standardised, cleaned and capitalised protein name.

    Raises
    ------
    ValueError
        If ``target`` is not a known database target.
    """
    try:
        standardizer = _STANDARDIZERS[target]
    except KeyError:
        raise ValueError(f"Unknown target type: {target!r}")
    return standardizer(name, mol_type, taxonomy)
