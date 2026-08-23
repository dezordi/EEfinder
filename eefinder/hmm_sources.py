"""Prebuilt profile-HMM sources, fetched and mapped onto EEfinder metadata.

The ``ncbi-refseq`` source *builds* profiles from a RefSeq download; a prebuilt
source instead *downloads* somebody else's profiles and has to reconstruct the
metadata EEfinder needs from the tables that ship with them. Two fields are
mandatory (see ``proposal.md`` §1.1): viral taxonomy at least to genus, and a
protein name. ``Host`` is optional and filled ``Unknown``.

`VOGDB <https://vogdb.org>`_ is the one prebuilt set that supplies both, plus a
flag for the profiles that also occur in cellular genomes:

======================== ==================================================
file                     what it supplies
======================== ==================================================
``vog.annotations.tsv``  ``Protein`` (consensus description + category)
``vog.lca.tsv``          ``Genus``/``Family`` from the LCA lineage
``vog.members.tsv``      member proteins, i.e. the species of a group
``vog.virusonly.tsv``    "only in viruses" flags (host-contamination filter)
``vogdb.species.txt``    species name per taxid
``vogdb.host.txt``       ``phage``/``nonphage`` + host per taxid
``vog.hmm.tar.gz``       the profiles themselves (554 MB)
======================== ==================================================

Only the last one is large. Every filter below runs on the small tables, so
``--metadata-only`` can report exactly how many profiles a configuration would
keep for a few MB of download -- which is how a storage-constrained machine sizes
a database before committing the disk to it.

Two more sources are supported, each with a different metadata problem:

* **RVDB-prot** ships a SQLite annotation database (6 MB) with one row per family
  (``size``, ``nbseq``, ``LCAtaxid``) plus keyword tables, but **no taxon names**
  -- only NCBI taxids. Those are resolved through the ``datasets`` CLI EEfinder
  already depends on (1,457 distinct taxids for release 32.0, cached after the
  first run), which also supplies ``genomic_moltype``. Its "protein name" is a
  bag of scored keywords rather than a phrase, so
  :func:`protein_from_keywords` composes one.
* **NeoRdRp** ships a 179-column InterProScan/Palmscan annotation of its seed
  sequences that contains **no taxonomy field at all**. What can be recovered
  comes from the *names* of the CDD and RdRp-scan signatures that matched
  (``ps-ssRNAv_Solemoviridae_RdRp``), which is indirect and incomplete: of the
  19,394 profiles in release 2.1, 17,669 carry any annotation at all, 9,980
  resolve to a family and only 1,346 to a genus. It is
  therefore an RdRp-sensitivity booster to merge with another source, not a
  standalone database -- with the default ``--min-lca-rank genus`` most of it is
  dropped, which is the honest outcome rather than a bug.

Two caveats this module encodes rather than hides:

* a VOG's taxonomy is its **LCA**, which for broad groups sits above family; the
  ``min_lca_rank`` filter is what keeps the mandatory genus/family field honest;
* ``vogdb.host.txt`` fills ``host``/``superkingdom`` sparsely, so the eukaryotic
  filter keys on the ``phage``/``nonphage`` flag (cross-checked against the
  lineage), not on the host columns.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import lzma
import os
import re
import sqlite3
import subprocess
import tarfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field

from eefinder.get_databases import _genus_family_from_lineage
from eefinder.log import logger

#: Mirrors serving the same VOGDB tree; tried in order.
VOGDB_MIRRORS = (
    "https://fileshare.lisc.univie.ac.at/vog",
    "https://fileshare.csb.univie.ac.at/vog",
)

#: Small tables every filter and the metadata CSV are built from (~6 MB total).
VOGDB_METADATA_FILES = (
    "vog.annotations.tsv.gz",
    "vog.lca.tsv.gz",
    "vog.members.tsv.gz",
    "vog.virusonly.tsv.gz",
    "vogdb.species.txt",
    "vogdb.host.txt",
)

#: The profile archive -- the only large download (554 MB as of release 236).
VOGDB_PROFILE_ARCHIVE = "vog.hmm.tar.gz"

#: File recording the release number of a VOGDB tree.
VOGDB_RELEASE_FILE = "release.txt"

#: Virus-only stringency levels, mapped to their column in ``vog.virusonly.tsv``.
VIRUSONLY_COLUMNS = {"high": 1, "medium": 2, "low": 3}

#: Taxonomic ranks a profile's LCA can be required to resolve to.
LCA_RANKS = ("family", "genus", "none")

#: Lineage names that mark a prokaryotic-virus (phage) clade. Used to
#: cross-check the per-species phage flag, and as the *only* eukaryotic filter
#: for the sources that publish no host information at all (RVDB, NeoRdRp).
#: The RNA-phage families matter for NeoRdRp in particular: ``Leviviridae`` and
#: its ICTV successors are the single largest group in its profile set.
PHAGE_CLADES = (
    # dsDNA phages
    "caudoviricetes",
    "caudovirales",
    "myoviridae",
    "siphoviridae",
    "podoviridae",
    "microviridae",
    "inoviridae",
    "tectiviridae",
    "corticoviridae",
    "plasmaviridae",
    "tubulavirales",
    # RNA phages (Leviviridae was split into these by ICTV)
    "leviviricetes",
    "leviviridae",
    "fiersviridae",
    "steitzviridae",
    "blumeviridae",
    "duinviridae",
    "solspiviridae",
    "cystoviridae",
    # archaeal viruses
    "fuselloviridae",
    "lipothrixviridae",
    "rudiviridae",
    "bicaudaviridae",
)

#: UniProt-style prefix VOGDB puts in front of many consensus descriptions
#: (``sp|P13316|GRCA_BPT4 Autonomous glycyl radical cofactor``): 4,666 of the
#: 49,116 groups of release 236 carry one, and a far higher share of the
#: eukaryotic-virus groups do. The accession is provenance, not a protein name,
#: and its ``|`` characters would travel into the taxonomy table.
_UNIPROT_PREFIX_RE = re.compile(r"^\s*(?:sp|tr)\|[^|]*\|\S+\s+", re.IGNORECASE)


def clean_vogdb_description(description: str) -> str:
    """Reduce a VOGDB consensus description to a protein name.

    Strips the UniProt ``sp|ACC|ENTRY`` provenance prefix and any remaining
    ``|`` (which :mod:`eefinder.bed` uses as a field separator), leaving the text
    :func:`~eefinder.normalization.standardize_protein` can canonicalise.

    Examples
    --------
    >>> clean_vogdb_description("sp|P84400|MB43_EHV1V Membrane protein UL43 homolog")
    'Membrane protein UL43 homolog'
    >>> clean_vogdb_description("terminase large subunit")
    'terminase large subunit'
    """
    cleaned = _UNIPROT_PREFIX_RE.sub("", description or "")
    return cleaned.replace("|", " ").strip()


#: Functional categories / descriptions that carry no protein information.
VOGDB_UNINFORMATIVE = ("xu",)
VOGDB_UNINFORMATIVE_DESCRIPTIONS = ("hypothetical", "uncharacterized", "unknown")


@dataclass
class ProfileRecord:
    """Everything known about one prebuilt profile before it becomes a row.

    Shared by every prebuilt source so one set of filters
    (:func:`filter_vogs`) applies to all of them; a source simply leaves the
    fields it cannot supply empty.
    """

    group: str
    protein_count: int = 0
    species_count: int = 0
    category: str = ""
    description: str = ""
    lineage: "list[str]" = field(default_factory=list)
    lca_taxid: str = ""
    member_taxids: "list[str]" = field(default_factory=list)
    virus_only: "dict[str, bool]" = field(default_factory=dict)
    #: Filled by the sources that can supply them directly (RVDB via the
    #: ``datasets`` taxonomy lookup); otherwise derived from ``lineage``.
    species: str = ""
    mol_type: str = ""

    @property
    def genus_family(self) -> "tuple[str, str]":
        """``(genus, family)`` inferred from the LCA lineage by ICTV suffix."""
        return _genus_family_from_lineage([{"name": name} for name in self.lineage])

    @property
    def lca_rank(self) -> str:
        """Finest ICTV rank the LCA resolves to: ``genus``/``family``/``none``."""
        genus, family = self.genus_family
        if genus:
            return "genus"
        if family:
            return "family"
        return "none"

    def is_phage_clade(self) -> bool:
        """Whether the LCA lineage names a known prokaryotic-virus clade."""
        lowered = [name.lower() for name in self.lineage]
        return any(clade in lowered for clade in PHAGE_CLADES)

    def is_uninformative(self) -> bool:
        """Whether the group's functional annotation says nothing."""
        if self.category.lower() in VOGDB_UNINFORMATIVE:
            return True
        lowered = self.description.lower()
        return any(term in lowered for term in VOGDB_UNINFORMATIVE_DESCRIPTIONS)


#: Backwards-compatible alias -- the record type used to be VOGDB-specific.
VogRecord = ProfileRecord


# --------------------------------------------------------------------------- #
# Download
# --------------------------------------------------------------------------- #
def _fetch(url: str, dest: str) -> None:
    """Download ``url`` to ``dest`` (no-op when the file is already there)."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        logger.debug(f"reusing already downloaded {dest}")
        return
    logger.debug(f"downloading {url}")
    with urllib.request.urlopen(url) as response, open(dest, "wb") as out:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def _md5(path: str) -> str:
    """Return the MD5 checksum of a file."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify(path: str, md5_path: str) -> bool:
    """Compare a download against its ``.md5`` sidecar, if one was fetched."""
    if not os.path.exists(md5_path):
        return True
    expected = open(md5_path).read().split()[0].strip()
    actual = _md5(path)
    if expected != actual:
        logger.warning(f"{path} failed its MD5 check ({actual} != {expected})")
        return False
    return True


def download_vogdb(
    workdir: str,
    release: str = "latest",
    metadata_only: bool = True,
    mirrors: "tuple[str, ...]" = VOGDB_MIRRORS,
) -> "dict[str, str]":
    """Fetch the VOGDB tables (and, unless ``metadata_only``, the profiles).

    Parameters
    ----------
    workdir : str
        Directory the files are downloaded into (existing files are reused, so a
        re-run after a failure does not start over).
    release : str
        ``"latest"`` or a release number (e.g. ``"236"``) to pin.
    metadata_only : bool
        Skip :data:`VOGDB_PROFILE_ARCHIVE`, the only large file.
    mirrors : tuple[str, ...]
        Base URLs to try in order.

    Returns
    -------
    dict[str, str]
        ``filename -> local path`` for everything fetched.

    Raises
    ------
    RuntimeError
        If no mirror served a required file.
    """
    os.makedirs(workdir, exist_ok=True)
    wanted = list(VOGDB_METADATA_FILES)
    if not metadata_only:
        wanted.append(VOGDB_PROFILE_ARCHIVE)

    paths: dict[str, str] = {}
    errors: list[str] = []
    for base in mirrors:
        root = f"{base}/{'latest' if release == 'latest' else f'vog{release}'}"
        try:
            for name in [VOGDB_RELEASE_FILE] + wanted:
                dest = os.path.join(workdir, name)
                _fetch(f"{root}/{name}", dest)
                try:
                    _fetch(f"{root}/{name}.md5", f"{dest}.md5")
                except urllib.error.URLError:
                    pass  # not every mirror serves the sidecars
                if not _verify(dest, f"{dest}.md5"):
                    raise RuntimeError(f"{name} is corrupt on {base}")
                paths[name] = dest
            logger.info(f"Downloaded VOGDB {release} from {base}")
            return paths
        except (urllib.error.URLError, RuntimeError, OSError) as err:
            errors.append(f"{base}: {err}")
            logger.debug(f"VOGDB mirror {base} failed: {err}")
    raise RuntimeError("could not download VOGDB from any mirror: " + "; ".join(errors))


def vogdb_release(workdir: str) -> str:
    """Return the release number recorded in a downloaded VOGDB tree."""
    path = os.path.join(workdir, VOGDB_RELEASE_FILE)
    if not os.path.exists(path):
        return "unknown"
    return open(path).read().strip().splitlines()[0]


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def _open_table(path: str):
    """Open a VOGDB table, transparently handling the gzipped ones."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _rows(path: str) -> "list[list[str]]":
    """Read a headed, tab-separated VOGDB table into rows (header dropped)."""
    with _open_table(path) as handle:
        return [
            line.rstrip("\n").split("\t")
            for line in handle
            if line.strip() and not line.startswith("#")
        ]


def parse_annotations(path: str) -> "dict[str, ProfileRecord]":
    """Parse ``vog.annotations.tsv``: counts, functional category, description."""
    records: dict[str, VogRecord] = {}
    for row in _rows(path):
        if len(row) < 5:
            continue
        records[row[0]] = ProfileRecord(
            group=row[0],
            protein_count=int(row[1]) if row[1].isdigit() else 0,
            species_count=int(row[2]) if row[2].isdigit() else 0,
            category=row[3],
            description=row[4],
        )
    return records


def parse_lca(path: str) -> "dict[str, tuple[list[str], str]]":
    """Parse ``vog.lca.tsv`` into ``group -> (lineage names, LCA taxid)``."""
    lineages: dict[str, tuple[list[str], str]] = {}
    for row in _rows(path):
        if len(row) < 5:
            continue
        names = [name.strip() for name in row[3].split(";") if name.strip()]
        lineages[row[0]] = (names, row[4])
    return lineages


def parse_members(path: str) -> "dict[str, list[str]]":
    """Parse ``vog.members.tsv`` into ``group -> [member taxids]``.

    Member ids look like ``2178927.YP_009889867.1``; only the leading taxid is
    needed, since it is the join key to the species and host tables.
    """
    members: dict[str, list[str]] = {}
    for row in _rows(path):
        if len(row) < 5:
            continue
        members[row[0]] = [
            protein.split(".", 1)[0] for protein in row[4].split(",") if protein
        ]
    return members


def parse_virusonly(path: str) -> "dict[str, dict[str, bool]]":
    """Parse ``vog.virusonly.tsv`` into ``group -> {stringency: only in viruses}``."""
    flags: dict[str, dict[str, bool]] = {}
    for row in _rows(path):
        if len(row) < 4:
            continue
        flags[row[0]] = {
            level: row[column] == "1" for level, column in VIRUSONLY_COLUMNS.items()
        }
    return flags


def parse_species(path: str) -> "dict[str, str]":
    """Parse ``vogdb.species.txt`` into ``taxid -> species name``."""
    species: dict[str, str] = {}
    for row in _rows(path):
        if len(row) < 2:
            continue
        species[row[1]] = row[0]
    return species


def parse_hosts(path: str) -> "dict[str, tuple[bool | None, str]]":
    """Parse ``vogdb.host.txt`` into ``taxid -> (is_phage, host name)``.

    ``is_phage`` is ``None`` when the table does not say; the ``host`` column is
    sparsely filled, so an empty host is normal and becomes ``Unknown`` later.
    """
    hosts: dict[str, tuple[bool | None, str]] = {}
    for row in _rows(path):
        if len(row) < 2:
            continue
        flag = row[1].strip().lower()
        is_phage = True if flag == "phage" else False if flag == "nonphage" else None
        hosts[row[0]] = (is_phage, row[2].strip() if len(row) > 2 else "")
    return hosts


def load_vogdb(paths: "dict[str, str]") -> "dict[str, ProfileRecord]":
    """Merge the downloaded VOGDB tables into one record per group."""
    records = parse_annotations(paths["vog.annotations.tsv.gz"])
    lineages = parse_lca(paths["vog.lca.tsv.gz"])
    members = parse_members(paths["vog.members.tsv.gz"])
    flags = parse_virusonly(paths["vog.virusonly.tsv.gz"])
    for group, record in records.items():
        record.lineage, record.lca_taxid = lineages.get(group, ([], ""))
        record.member_taxids = members.get(group, [])
        record.virus_only = flags.get(group, {})
    logger.debug(f"VOGDB: loaded {len(records)} group(s)")
    return records


# --------------------------------------------------------------------------- #
# Filtering
# --------------------------------------------------------------------------- #
@dataclass
class VogdbFilters:
    """The filters applied to a VOGDB release before profiles are extracted."""

    taxon: str = "10239"
    eukaryotic_only: bool = True
    euk_fraction: float = 0.9
    virus_only_stringency: str = "medium"
    min_lca_rank: str = "genus"
    min_profile_seqs: int = 3
    exclude_uninformative: bool = True


@dataclass
class FilterCounts:
    """How many groups each filter removed (reported in the build log)."""

    total: int = 0
    dropped_taxon: int = 0
    dropped_not_eukaryotic: int = 0
    dropped_not_virus_only: int = 0
    dropped_lca_rank: int = 0
    dropped_too_small: int = 0
    dropped_uninformative: int = 0
    kept: int = 0


def _eukaryotic_fraction(
    record: ProfileRecord, hosts: "dict[str, tuple[bool | None, str]]"
) -> "float | None":
    """Fraction of member species flagged non-phage (``None`` when unflagged)."""
    flags = [
        hosts[taxid][0]
        for taxid in record.member_taxids
        if taxid in hosts and hosts[taxid][0] is not None
    ]
    if not flags:
        return None
    return sum(1 for flag in flags if not flag) / len(flags)


def filter_vogs(
    records: "dict[str, ProfileRecord]",
    hosts: "dict[str, tuple[bool | None, str]]",
    filters: VogdbFilters,
) -> "tuple[list[ProfileRecord], FilterCounts]":
    """Apply the configured filters, returning the survivors and the tallies.

    The order matters only for the counts, which attribute each dropped group to
    the first filter that rejected it.
    """
    counts = FilterCounts(total=len(records))
    kept: list[ProfileRecord] = []
    taxon = filters.taxon.strip().lower()
    wants_all = taxon in {"", "10239", "viruses"}

    for record in records.values():
        lineage = [name.lower() for name in record.lineage]
        if not wants_all and taxon not in lineage:
            counts.dropped_taxon += 1
            continue

        if filters.eukaryotic_only:
            fraction = _eukaryotic_fraction(record, hosts)
            if record.is_phage_clade() or (
                fraction is not None and fraction < filters.euk_fraction
            ):
                counts.dropped_not_eukaryotic += 1
                continue

        stringency = filters.virus_only_stringency
        # Only VOGDB publishes these flags; for a source that does not, the test
        # is not applicable rather than failed.
        if (
            stringency != "none"
            and record.virus_only
            and not record.virus_only.get(stringency, False)
        ):
            counts.dropped_not_virus_only += 1
            continue

        if filters.min_lca_rank != "none":
            rank = record.lca_rank
            if rank == "none" or (filters.min_lca_rank == "genus" and rank != "genus"):
                counts.dropped_lca_rank += 1
                continue

        if record.protein_count < filters.min_profile_seqs:
            counts.dropped_too_small += 1
            continue

        if filters.exclude_uninformative and record.is_uninformative():
            counts.dropped_uninformative += 1
            continue

        kept.append(record)

    counts.kept = len(kept)
    logger.debug(f"VOGDB filters: {counts}")
    return kept, counts


# --------------------------------------------------------------------------- #
# Profile extraction
# --------------------------------------------------------------------------- #
def extract_profiles(
    archive_path: str, groups: "set[str]", outdir: str, namespace: str = "VOGDB"
) -> "dict[str, str]":
    """Extract the kept profiles from ``vog.hmm.tar.gz``, renaming each one.

    The ``NAME`` line is rewritten to ``{namespace}__{group}`` so profile ids stay
    unique (and say where they came from) when several databases are searched
    together.

    Parameters
    ----------
    archive_path : str
        The downloaded ``vog.hmm.tar.gz``.
    groups : set[str]
        Group names (``VOG00001``) to extract.
    outdir : str
        Directory the individual profiles are written to.
    namespace : str
        Id namespace prefix.

    Returns
    -------
    dict[str, str]
        ``group -> path`` of the extracted profiles.
    """
    os.makedirs(outdir, exist_ok=True)
    extracted: dict[str, str] = {}
    with tarfile.open(archive_path, "r:*") as archive:
        for member in archive:
            if not member.isfile():
                continue
            group = os.path.basename(member.name)
            if not group.endswith(".hmm"):
                continue
            group = group[: -len(".hmm")]
            if group not in groups:
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            text = handle.read().decode()
            profile_id = f"{namespace}__{group}"
            text = _rename_profile(text, profile_id)
            path = os.path.join(outdir, f"{profile_id}.hmm")
            with open(path, "w") as out:
                out.write(text)
            extracted[group] = path
    logger.debug(f"VOGDB: extracted {len(extracted)} of {len(groups)} profile(s)")
    return extracted


def _rename_profile(text: str, profile_id: str) -> str:
    """Rewrite the ``NAME`` line of a single HMMER profile."""
    lines = text.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.startswith("NAME "):
            lines[index] = f"NAME  {profile_id}\n"
            break
    return "".join(lines)


def profile_length(hmm_path: str) -> int:
    """Return the ``LENG`` of a single profile (0 when unreadable)."""
    with open(hmm_path) as handle:
        for line in handle:
            if line.startswith("LENG"):
                return int(line.split()[1])
            if line.startswith("HMM "):
                break
    return 0


# --------------------------------------------------------------------------- #
# Shared helpers for the sources that give taxids instead of names
# --------------------------------------------------------------------------- #
#: NCBI datasets CLI, already required by ``get-databases``.
DATASETS_BINARY = "datasets"

#: How many taxids to resolve per ``datasets`` invocation.
TAXONOMY_BATCH = 200

#: Ranks pulled out of a ``datasets`` taxonomy classification.
_WANTED_RANKS = ("species", "genus", "family")


def resolve_taxids(
    taxids: "list[str]", cache_path: "str | None" = None
) -> "dict[str, dict[str, str]]":
    """Resolve NCBI taxids to ``{species, genus, family, mol_type}`` names.

    RVDB-prot records the LCA of each family as a **taxid**, so the names the
    mandatory taxonomy field needs have to be looked up. This uses the
    ``datasets`` CLI EEfinder already depends on rather than a 100 MB taxonomy
    dump, and caches the result next to the download because the lookup is the
    slow part (roughly a second per taxid, and RVDB 32.0 has 1,457 distinct ones).

    Parameters
    ----------
    taxids : list[str]
        Taxids to resolve; duplicates and invalid entries are tolerated.
    cache_path : str, optional
        JSON file to read already-resolved taxids from and write new ones to.

    Returns
    -------
    dict[str, dict[str, str]]
        ``taxid -> {"species", "genus", "family", "mol_type"}``; taxids NCBI does
        not recognise are simply absent.
    """
    resolved: dict[str, dict[str, str]] = {}
    if cache_path and os.path.exists(cache_path):
        try:
            with open(cache_path) as handle:
                resolved = json.load(handle)
            logger.debug(f"taxonomy cache: {len(resolved)} taxid(s) from {cache_path}")
        except (json.JSONDecodeError, OSError):
            resolved = {}

    pending = sorted(
        {str(taxid) for taxid in taxids if str(taxid).isdigit() and str(taxid) != "0"}
        - set(resolved)
    )
    if pending:
        logger.info(f"Resolving {len(pending)} taxid(s) through the datasets CLI")
    for start in range(0, len(pending), TAXONOMY_BATCH):
        batch = pending[start : start + TAXONOMY_BATCH]
        result = subprocess.run(
            [
                DATASETS_BINARY,
                "summary",
                "taxonomy",
                "taxon",
                *batch,
                "--as-json-lines",
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logger.warning(
                f"datasets taxonomy lookup failed for a batch: "
                f"{result.stderr.strip()[:200]}"
            )
            continue
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            try:
                payload = json.loads(line)["taxonomy"]
            except (json.JSONDecodeError, KeyError):
                continue
            classification = payload.get("classification", {})
            entry = {
                rank: classification.get(rank, {}).get("name", "")
                for rank in _WANTED_RANKS
            }
            entry["mol_type"] = payload.get("genomic_moltype", "")
            taxid = str(payload.get("tax_id", "")) or str(
                payload.get("current_scientific_name", {}).get("tax_id", "")
            )
            for queried in payload.get("query", []) or [taxid]:
                if queried:
                    resolved[str(queried)] = entry
        logger.debug(f"resolved {start + len(batch)}/{len(pending)} taxid(s)")

    if cache_path and pending:
        with open(cache_path, "w") as handle:
            json.dump(resolved, handle)
    return resolved


def extract_profiles_from_hmm(
    hmm_path: str, groups: "dict[str, str]", outdir: str
) -> "dict[str, str]":
    """Split a concatenated ``.hmm`` (optionally ``.xz``) into selected profiles.

    RVDB-prot and NeoRdRp ship one compressed file holding every profile, rather
    than VOGDB's tar of individual ones. Each kept profile is written out on its
    own with its ``NAME`` rewritten to the namespaced id.

    Parameters
    ----------
    hmm_path : str
        The downloaded ``.hmm`` or ``.hmm.xz``.
    groups : dict[str, str]
        ``original profile name -> namespaced profile id`` for the profiles to
        keep.
    outdir : str
        Directory the individual profiles are written to.

    Returns
    -------
    dict[str, str]
        ``original name -> path`` of the extracted profiles.
    """
    os.makedirs(outdir, exist_ok=True)
    extracted: dict[str, str] = {}
    opener = lzma.open if hmm_path.endswith(".xz") else open
    block: list[str] = []
    name = ""
    with opener(hmm_path, "rt") as handle:
        for line in handle:
            block.append(line)
            if line.startswith("NAME "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("//"):
                profile_id = groups.get(name)
                if profile_id:
                    path = os.path.join(outdir, f"{profile_id}.hmm")
                    with open(path, "w") as out:
                        out.write(_rename_profile("".join(block), profile_id))
                    extracted[name] = path
                block, name = [], ""
    logger.debug(
        f"extracted {len(extracted)} of {len(groups)} profile(s) from {hmm_path}"
    )
    return extracted


# --------------------------------------------------------------------------- #
# RVDB-prot
# --------------------------------------------------------------------------- #
#: Where RVDB-prot publishes its releases.
RVDB_BASE_URL = "https://rvdb-prot.pasteur.fr/files"

#: Default RVDB-prot release (the newest at the time of writing).
RVDB_DEFAULT_RELEASE = "32.0"

#: Keywords that say nothing about which protein a family is.
RVDB_STOPWORDS = frozenset(
    {
        "viral",
        "virus",
        "viruses",
        "protein",
        "proteins",
        "putative",
        "hypothetical",
        "unnamed",
        "uncharacterized",
        "uncharacterised",
        "predicted",
        "product",
        "orf",
        "gene",
        "partial",
        "family",
        "like",
        "domain",
        "containing",
        "unknown",
        "fragment",
        "isolate",
        "strain",
    }
)

#: Keyword sets that identify a canonical protein name. RVDB's annotation is a
#: bag of scored single tokens rather than a phrase, so the canonical-name map
#: (which matches phrases) cannot be applied directly; these token sets bridge the
#: two for the most frequent viral functions.
RVDB_KEYWORD_RULES = (
    ({"rna", "dependent", "polymerase"}, "RdRp"),
    ({"rna", "directed", "polymerase"}, "RdRp"),
    ({"rdrp"}, "RdRp"),
    ({"reverse", "transcriptase"}, "Reverse Transcriptase"),
    ({"nucleocapsid"}, "Nucleocapsid Protein"),
    ({"nucleoprotein"}, "Nucleocapsid Protein"),
    ({"glycoprotein"}, "Glycoprotein"),
    ({"capsid"}, "Capsid Protein"),
    ({"coat"}, "Capsid Protein"),
    ({"polyprotein"}, "Polyprotein"),
    ({"integrase"}, "integrase (IN)"),
    ({"helicase"}, "Helicase"),
    ({"protease"}, "Main Protease"),
    ({"terminase"}, "Terminase Large Subunit"),
    ({"movement"}, "Movement protein"),
)

#: How many informative keywords are joined when no rule matches.
RVDB_KEYWORDS_IN_NAME = 3


def rvdb_files(release: str = RVDB_DEFAULT_RELEASE) -> "dict[str, str]":
    """Return the ``{local name: URL}`` map of an RVDB-prot release."""
    return {
        "rvdb.sqlite.xz": f"{RVDB_BASE_URL}/U-RVDBv{release}-prot-hmm.sqlite.xz",
        "rvdb.hmm.xz": f"{RVDB_BASE_URL}/U-RVDBv{release}-prot.hmm.xz",
    }


def protein_from_keywords(keywords: "list[tuple[str, int]]") -> str:
    """Compose a protein name from RVDB's scored keyword bag.

    A rule from :data:`RVDB_KEYWORD_RULES` wins when its whole token set is
    present, giving the same canonical name the phrase map would; otherwise the
    top :data:`RVDB_KEYWORDS_IN_NAME` informative keywords are joined in
    frequency order.

    Parameters
    ----------
    keywords : list[tuple[str, int]]
        ``(keyword, frequency)`` pairs, any order.

    Returns
    -------
    str
        A protein name, or ``"Unknown"`` when every keyword is uninformative.

    Examples
    --------
    >>> protein_from_keywords([("rna", 2311), ("polymerase", 1000), ("dependent", 863)])
    'RdRp'
    >>> protein_from_keywords([("tail", 40), ("fiber", 38), ("viral", 90)])
    'Tail fiber'
    """
    ranked = sorted(keywords, key=lambda item: (-item[1], item[0]))
    tokens = {word.lower() for word, _ in ranked}
    for required, canonical in RVDB_KEYWORD_RULES:
        if required <= tokens:
            return canonical

    informative = [
        word
        for word, _ in ranked
        if word.lower() not in RVDB_STOPWORDS and word.strip()
    ]
    if not informative:
        return "Unknown"
    name = " ".join(informative[:RVDB_KEYWORDS_IN_NAME])
    return name[0].upper() + name[1:]


def download_rvdb(
    workdir: str,
    release: str = RVDB_DEFAULT_RELEASE,
    metadata_only: bool = True,
) -> "dict[str, str]":
    """Download an RVDB-prot release (annotations always, profiles on demand)."""
    os.makedirs(workdir, exist_ok=True)
    paths: dict[str, str] = {}
    for name, url in rvdb_files(release).items():
        if metadata_only and name.endswith(".hmm.xz"):
            continue
        dest = os.path.join(workdir, name)
        _fetch(url, dest)
        paths[name] = dest
    logger.info(f"Downloaded RVDB-prot v{release} into {workdir}")
    return paths


def load_rvdb(
    sqlite_path: str, cache_path: "str | None" = None
) -> "dict[str, ProfileRecord]":
    """Read the RVDB-prot annotation database into profile records.

    ``Protein`` comes from the Pfam-derived keyword table when a family has one
    (``fam_kw_ref``) and from the sequence-name keywords otherwise; the lineage is
    resolved from the family's LCA taxid through :func:`resolve_taxids`.
    """
    path = _decompress(sqlite_path)
    connection = sqlite3.connect(path)
    try:
        families = {
            str(row[0]): ProfileRecord(
                group=f"FAM{int(row[0]):06d}",
                protein_count=int(row[2] or 0),
                species_count=0,
                lca_taxid=str(row[3] or ""),
            )
            for row in connection.execute(
                "select id, size, nbseq, LCAtaxid from family"
            )
        }
        keywords: dict[str, list[tuple[str, int]]] = {}
        for table in ("fam_kw_ref", "fam_kw_seqnames"):
            for fam_id, word, freq in connection.execute(
                f"select f.famId, k.str, f.freq from {table} f "
                "join keyword k on k.id = f.kwId"
            ):
                fam = str(fam_id)
                if table == "fam_kw_seqnames" and keywords.get(fam):
                    continue  # the Pfam-derived keywords already named it
                keywords.setdefault(fam, []).append((word, int(freq or 0)))
    finally:
        connection.close()

    taxonomy = resolve_taxids(
        [record.lca_taxid for record in families.values()], cache_path=cache_path
    )
    for fam_id, record in families.items():
        record.description = protein_from_keywords(keywords.get(fam_id, []))
        entry = taxonomy.get(record.lca_taxid, {})
        record.species = entry.get("species", "")
        record.mol_type = entry.get("mol_type", "")
        record.lineage = [
            name for name in (entry.get("family", ""), entry.get("genus", "")) if name
        ]
    logger.debug(f"RVDB: loaded {len(families)} family(ies) from {sqlite_path}")
    return {record.group: record for record in families.values()}


# --------------------------------------------------------------------------- #
# NeoRdRp
# --------------------------------------------------------------------------- #
#: Zenodo record holding NeoRdRp 2.1.
NEORDRP_BASE_URL = "https://zenodo.org/records/10851672/files"

#: Default NeoRdRp release.
NEORDRP_DEFAULT_RELEASE = "2.1"

#: Columns of the seed annotation the taxonomy tokens are mined from.
NEORDRP_TAXONOMY_COLUMNS = (
    "InterProScan_CDD_Signature Description",
    "hmmsearch_RdRp-scan_description_of_target",
)

#: Column holding the profile a seed sequence belongs to.
NEORDRP_HMM_COLUMN = "Seed_RdRp_HMM Name"

#: Every NeoRdRp profile models the same protein.
NEORDRP_PROTEIN = "RdRp"

#: ...of the same genome type.
NEORDRP_MOL_TYPE = "RNA"


def neordrp_files(release: str = NEORDRP_DEFAULT_RELEASE) -> "dict[str, str]":
    """Return the ``{local name: URL}`` map of a NeoRdRp release."""
    return {
        "neordrp.annotation.tsv.xz": (
            f"{NEORDRP_BASE_URL}/Annotation_of_Seed_RdRp_datasets.tsv.xz?download=1"
        ),
        "neordrp.hmm.xz": f"{NEORDRP_BASE_URL}/NeoRdRp.{release}.hmm.xz?download=1",
    }


def download_neordrp(
    workdir: str,
    release: str = NEORDRP_DEFAULT_RELEASE,
    metadata_only: bool = True,
) -> "dict[str, str]":
    """Download NeoRdRp (annotation always, profiles on demand)."""
    os.makedirs(workdir, exist_ok=True)
    paths: dict[str, str] = {}
    for name, url in neordrp_files(release).items():
        if metadata_only and name.endswith(".hmm.xz"):
            continue
        dest = os.path.join(workdir, name)
        _fetch(url, dest)
        paths[name] = dest
    logger.info(f"Downloaded NeoRdRp {release} into {workdir}")
    return paths


def load_neordrp(annotation_path: str) -> "dict[str, ProfileRecord]":
    """Read the NeoRdRp seed annotation into one record per profile.

    The annotation has no taxonomy field, so the only taxonomic signal available
    is in the **names** of the CDD and RdRp-scan signatures that matched a seed
    (``ps-ssRNAv_Solemoviridae_RdRp``). Those tokens are collected per profile and
    read with the same ICTV suffix heuristic used elsewhere. This is inference,
    not an authoritative assignment, and it leaves a large share of the profiles
    without a genus -- which the ``--min-lca-rank`` filter then drops.

    The file has one row per **seed sequence**, so a profile's sequence count is
    how many of its seeds are annotated; without it every profile would look
    empty to the ``--min-profile-seqs`` filter.
    """
    genera: dict[str, dict[str, int]] = {}
    families: dict[str, dict[str, int]] = {}
    seeds: dict[str, int] = {}
    informative: dict[str, int] = {}
    opener = lzma.open if annotation_path.endswith(".xz") else open
    with opener(annotation_path, "rt") as handle:
        header = next(handle).rstrip("\n").split("\t")
        try:
            hmm_column = header.index(NEORDRP_HMM_COLUMN)
        except ValueError:
            raise ValueError(
                f"{annotation_path} has no '{NEORDRP_HMM_COLUMN}' column; the "
                "NeoRdRp annotation format may have changed."
            )
        taxonomy_columns = [
            header.index(name) for name in NEORDRP_TAXONOMY_COLUMNS if name in header
        ]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= hmm_column:
                continue
            profile = fields[hmm_column]
            if _NEORDRP_PLACEHOLDER_RE.match(profile):
                continue
            # One row per seed sequence, so the row count per profile is how many
            # sequences it was built from (a lower bound: only annotated seeds).
            seeds[profile] = seeds.get(profile, 0) + 1
            found = set()
            for column in taxonomy_columns:
                if column < len(fields) and fields[column] not in ("-", ""):
                    found.update(_taxonomic_tokens(fields[column]))
            if not found:
                continue
            informative[profile] = informative.get(profile, 0) + 1
            # Vote per seed, not per token: a profile whose seeds matched several
            # signatures would otherwise collect a genus and a family from
            # different lineages.
            genus, family = _genus_family_from_lineage(
                [{"name": name} for name in sorted(found)]
            )
            if genus:
                genera.setdefault(profile, {})[genus] = (
                    genera.setdefault(profile, {}).get(genus, 0) + 1
                )
            if family:
                families.setdefault(profile, {})[family] = (
                    families.setdefault(profile, {}).get(family, 0) + 1
                )

    records = {}
    for profile, count in seeds.items():
        voters = informative.get(profile, 0)
        genus = _majority(genera.get(profile, {}), voters, NEORDRP_TAXON_AGREEMENT)
        family = _majority(families.get(profile, {}), voters, NEORDRP_TAXON_AGREEMENT)
        records[profile] = ProfileRecord(
            group=profile,
            protein_count=count,
            description=NEORDRP_PROTEIN,
            mol_type=NEORDRP_MOL_TYPE,
            lineage=[name for name in (family, genus) if name],
        )
    logger.debug(f"NeoRdRp: loaded {len(records)} profile(s) from {annotation_path}")
    return records


#: Suffixes that mark an ICTV taxon name inside a signature description.
_TAXON_SUFFIXES = ("viridae", "virus", "virales", "viricetes", "virinae")

#: Values NeoRdRp uses to mean "this seed is in no profile". The annotation
#: spells it both as ``-`` and as ``-|-`` (1,399 seeds of release 2.1).
_NEORDRP_PLACEHOLDER_RE = re.compile(r"^[-|\s]*$")

#: Minimum share of a profile's **taxonomically informative** seeds that must
#: agree on a taxon for it to be adopted. Pooling every seed's tokens instead
#: would let one profile end up with a genus and a family from different
#: lineages; using all seeds as the denominator would instead punish a profile
#: for the seeds whose signatures simply say nothing.
NEORDRP_TAXON_AGREEMENT = 0.5


def _taxonomic_tokens(description: str) -> "set[str]":
    """Pull ICTV-looking taxon names out of a signature description.

    Tokens must be **capitalised**: ICTV taxon names are, while the generic words
    that share the suffixes (``coronavirus``, ``bacteriophage``) are not, and
    those would otherwise be adopted as genus names.
    """
    found = set()
    for token in re.split(r"[_\s;,|:()]+", description):
        token = token.strip()
        if (
            len(token) > 5
            and token[:1].isupper()
            and token.lower().endswith(_TAXON_SUFFIXES)
        ):
            found.add(token)
    return found


def _majority(counts: "dict[str, int]", total: int, minimum: float) -> str:
    """Most frequent value when it holds ``minimum`` of ``total``, else empty."""
    if not counts or total <= 0:
        return ""
    name, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
    return name if count / total >= minimum else ""


def _decompress(path: str) -> str:
    """Decompress an ``.xz`` file next to itself and return the usable path."""
    if not path.endswith(".xz"):
        return path
    target = path[: -len(".xz")]
    if not os.path.exists(target):
        logger.debug(f"decompressing {path}")
        with lzma.open(path, "rb") as source, open(target, "wb") as out:
            while True:
                chunk = source.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
    return target
