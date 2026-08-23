"""Build profile HMMs from protein clusters (the ``ncbi-refseq``/host builder).

One cluster of homologous proteins becomes one profile HMM:

``cd-hit`` cluster -> ``mafft --auto`` -> **CIAlign curation** -> ``hmmbuild``

The curation step is not cosmetic. RefSeq entries include partial CDS,
polyprotein fragments and occasional long unique insertions; ``mafft`` places
those as ragged termini and sparse insertion columns, and ``hmmbuild`` would turn
them into match states with almost no information content -- profiles that score
noise as readily as signal. `CIAlign <https://cialign.readthedocs.io>`_
(Tumescheit *et al.* 2022) trims them and logs exactly what it removed.

Two deliberate choices, both reversible from the CLI:

* the cd-hit **representative is protected** from every row-wise removal
  (``--retain_str``), so a cluster can never lose the sequence its id and
  metadata derive from;
* ``--remove_divergent`` is **off by default**: it drops members sharing less
  than 65% of positions with the consensus, and in a family-level viral cluster
  that is precisely the remote homolog whose inclusion gives the profile its
  sensitivity advantage over ``blastx``.

Curation also never silently deletes a profile: if the curated alignment falls
below the sequence/column minimums, the builder falls back to the uncurated
alignment and says so in the debug log.

:class:`BuilderSettings` is recorded in the log of every database this module
builds, because the symmetric host filter compares a viral profile's score with a
host profile's score -- a comparison that is only meaningful when both databases
were produced with the same settings.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from dataclasses import asdict, dataclass, field

from Bio import SeqIO

from eefinder.log import logger

#: External binaries used by the builder (conda: ``mafft``, ``cialign``,
#: ``hmmer``).
MAFFT_BINARY = "mafft"
CIALIGN_BINARY = "CIAlign"
HMMBUILD_BINARY = "hmmbuild"

#: Minimum sequences for CIAlign to run: its column-wise functions
#: (``crop_ends``, ``remove_insertions``, ``remove_divergent``) compare each
#: sequence against the majority and refuse to run on fewer.
MIN_CURATION_SEQS = 3

#: Characters allowed in a profile id. Profile ids travel through the pipeline as
#: ``sseqid``, so ``|`` (used as a separator by :mod:`eefinder.bed`) and
#: whitespace (the ``--domtblout`` field separator) must never appear in one.
_ID_SAFE_RE = re.compile(r"[^A-Za-z0-9.-]+")


@dataclass
class BuilderSettings:
    """Every setting that affects the profiles a build produces.

    Serialised into each database's log so that ``screening`` can verify a viral
    and a host database are comparable before running the symmetric filter.
    """

    cluster_identity: float = 0.5
    min_profile_seqs: int = 3
    min_alignment_columns: int = 40
    keep_singletons: bool = False
    msa_curation: bool = True
    curation_crop_ends: bool = True
    curation_remove_insertions: bool = True
    curation_remove_short: bool = True
    curation_remove_min_length: int = 30
    curation_remove_divergent: bool = False
    curation_divergent_minperc: float = 0.25
    hmmer_version: str = ""

    def to_dict(self) -> dict:
        """Return the settings as a plain dict (for JSON logs)."""
        return asdict(self)

    def differences(self, other: "BuilderSettings") -> "list[str]":
        """Return the names of the settings that differ from ``other``.

        ``hmmer_version`` is compared too: a profile built by a different HMMER
        release is not guaranteed to score identically.
        """
        mine, theirs = self.to_dict(), other.to_dict()
        return [key for key in mine if mine[key] != theirs.get(key)]


@dataclass
class ClusterSpec:
    """A cluster of protein ids to turn into one profile."""

    profile_id: str
    representative: str
    members: "list[str]"


@dataclass
class BuiltProfile:
    """The outcome of building one profile."""

    profile_id: str
    hmm_path: str
    representative: str
    members: "list[str]" = field(default_factory=list)
    n_seqs: int = 0  #: sequences in the alignment ``hmmbuild`` consumed
    length: int = 0  #: profile length (``LENG``)
    curated: bool = False
    removed_sequences: int = 0
    removed_columns: int = 0


def slugify(name: str) -> str:
    """Turn a protein name into an id-safe slug.

    Spaces and separators collapse to ``_`` and every other unsafe character is
    dropped, so the slug is safe as part of a profile id.

    Examples
    --------
    >>> slugify("Glycoprotein (G)")
    'Glycoprotein_G'
    >>> slugify("RNA-dependent RNA polymerase")
    'RNA-dependent_RNA_polymerase'
    """
    slug = _ID_SAFE_RE.sub("_", name.strip())
    return slug.strip("_") or "Unknown"


def make_profile_id(source: str, taxon: str, protein: str, index: int) -> str:
    """Compose a namespaced, collision-safe profile id.

    Format ``{SOURCE}__{Taxon}__{ProteinSlug}__{NNN}``, e.g.
    ``NCBIREFSEQ__Chuviridae__Glycoprotein__001``.

    The zero-padded counter is not cosmetic:
    :meth:`eefinder.get_taxonomy.GetFinalTaxonomy._build_row` resolves accessions
    with a **substring** test, so an id must never be a substring of another
    (``…__1`` would otherwise also match ``…__10``).

    Parameters
    ----------
    source : str
        Database source namespace (e.g. ``"NCBIREFSEQ"``, ``"HOST"``).
    taxon : str
        Taxon the cluster was built inside (may be empty for host profiles).
    protein : str
        Consensus protein name.
    index : int
        1-based counter within the (source, taxon, protein) group.

    Returns
    -------
    str
    """
    parts = [slugify(source)]
    if taxon:
        parts.append(slugify(taxon))
    parts.append(slugify(protein))
    parts.append(f"{index:03d}")
    return "__".join(parts)


def hmmer_version() -> str:
    """Return the detected HMMER version (empty string when unavailable)."""
    try:
        result = subprocess.run([HMMBUILD_BINARY, "-h"], capture_output=True, text=True)
    except (FileNotFoundError, OSError):
        return ""
    match = re.search(r"HMMER\s+(\d+\.\d+(?:\.\d+)?)", result.stdout)
    return match.group(1) if match else ""


def alignment_dimensions(path: str) -> "tuple[int, int]":
    """Return ``(sequences, columns)`` of an aligned FASTA (0, 0 when empty)."""
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        return 0, 0
    return len(records), len(records[0].seq)


def write_member_fasta(
    sequences: "dict[str, str]", members: "list[str]", out_path: str
) -> int:
    """Write the cluster members to a FASTA, sorted by id for determinism.

    ``mafft`` output depends on input order, so sorting is what makes two builds
    of the same download byte-identical.

    Parameters
    ----------
    sequences : dict[str, str]
        ``protein_id -> sequence`` for (at least) every member.
    members : list[str]
        Member ids to write; ids missing from ``sequences`` are skipped.
    out_path : str
        Destination FASTA.

    Returns
    -------
    int
        Number of records written.
    """
    written = 0
    with open(out_path, "w") as out:
        for member in sorted(members):
            sequence = sequences.get(member)
            if sequence is None:
                logger.debug(f"cluster member {member} not found in the FASTA index")
                continue
            out.write(f">{member}\n{sequence}\n")
            written += 1
    return written


def run_mafft(faa_in: str, aln_out: str, threads: int = 1) -> None:
    """Align ``faa_in`` with ``mafft --auto``, writing aligned FASTA.

    Raises
    ------
    RuntimeError
        If ``mafft`` exits non-zero (stderr included).
    """
    command = f"{MAFFT_BINARY} --auto --quiet --thread {int(threads)} {faa_in}"
    logger.debug(f"mafft command: {command}")
    with open(aln_out, "w") as out:
        result = subprocess.run(
            shlex.split(command), stdout=out, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"mafft failed (exit {result.returncode}): {result.stderr.strip()}"
        )


def parse_cialign_removed(removed_path: str) -> "dict[str, list[str]]":
    """Parse a CIAlign ``{stem}_removed.txt`` into ``function -> [items]``.

    Items are sequence names for the row-wise functions and column indices for
    ``remove_insertions``. Used for the debug log only -- the counts reported by
    :func:`curate_alignment` are computed from the alignments themselves, so a
    change in this file's layout cannot corrupt them.
    """
    removals: dict[str, list[str]] = {}
    if not os.path.exists(removed_path):
        return removals
    with open(removed_path) as handle:
        for line in handle:
            function, _, items = line.rstrip("\n").partition("\t")
            if function:
                removals[function] = [i for i in items.split(",") if i]
    return removals


def curate_alignment(
    aln_in: str,
    stem: str,
    settings: BuilderSettings,
    retain_id: str = "",
) -> "tuple[str, int, int]":
    """Clean an alignment with CIAlign, falling back to the input when needed.

    Parameters
    ----------
    aln_in : str
        Aligned FASTA from :func:`run_mafft`.
    stem : str
        ``--outfile_stem`` for CIAlign (its outputs are ``{stem}_cleaned.fasta``
        and ``{stem}_removed.txt``).
    settings : BuilderSettings
        Which cleaning functions to run and with what parameters.
    retain_id : str
        Sequence protected from every row-wise removal (the cluster
        representative).

    Returns
    -------
    tuple[str, int, int]
        ``(alignment_path, removed_sequences, removed_columns)``. The path is the
        curated alignment, or ``aln_in`` when curation was skipped, failed, or
        would have taken the alignment below
        ``settings.min_profile_seqs``/``settings.min_alignment_columns``.
    """
    seqs_in, cols_in = alignment_dimensions(aln_in)
    if seqs_in < MIN_CURATION_SEQS:
        # CIAlign's column-wise functions need a majority to compare against and
        # refuse to run below three sequences; curating is a no-op here anyway.
        logger.debug(
            f"skipping curation of {aln_in}: {seqs_in} sequence(s) < "
            f"{MIN_CURATION_SEQS}"
        )
        return aln_in, 0, 0

    options = []
    if settings.curation_crop_ends:
        options.append("--crop_ends")
    if settings.curation_remove_insertions:
        options.append("--remove_insertions")
    if settings.curation_remove_short:
        options.append(
            f"--remove_short --remove_min_length {settings.curation_remove_min_length}"
        )
    if settings.curation_remove_divergent:
        options.append(
            "--remove_divergent --remove_divergent_minperc "
            f"{settings.curation_divergent_minperc}"
        )
    if not options:
        return aln_in, 0, 0
    if retain_id:
        options.append(f"--retain_str {shlex.quote(retain_id)}")

    command = (
        f"{CIALIGN_BINARY} --infile {aln_in} --outfile_stem {stem} "
        f"{' '.join(options)} --silent"
    )
    logger.debug(f"CIAlign command: {command}")
    result = subprocess.run(
        shlex.split(command),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    cleaned = f"{stem}_cleaned.fasta"
    if result.returncode != 0 or not os.path.exists(cleaned):
        logger.debug(
            f"CIAlign did not curate {aln_in} (exit {result.returncode}): "
            f"{result.stderr.strip()[:200]} -- keeping the uncurated alignment"
        )
        return aln_in, 0, 0

    seqs_out, cols_out = alignment_dimensions(cleaned)
    if (
        seqs_out < settings.min_profile_seqs
        or cols_out < settings.min_alignment_columns
    ):
        logger.debug(
            f"curation of {aln_in} left {seqs_out} sequence(s) / {cols_out} column(s), "
            f"below the minimums -- keeping the uncurated alignment"
        )
        return aln_in, 0, 0

    logger.debug(
        f"CIAlign {aln_in}: {seqs_in - seqs_out} sequence(s) and "
        f"{cols_in - cols_out} column(s) removed "
        f"({parse_cialign_removed(f'{stem}_removed.txt')})"
    )
    return cleaned, seqs_in - seqs_out, cols_in - cols_out


def run_hmmbuild(aln_path: str, hmm_out: str, name: str, threads: int = 1) -> None:
    """Build one profile from an alignment, naming it ``name`` (``hmmbuild -n``).

    Raises
    ------
    RuntimeError
        If ``hmmbuild`` exits non-zero (stderr included).
    """
    command = (
        f"{HMMBUILD_BINARY} -n {shlex.quote(name)} --amino --cpu {int(threads)} "
        f"-o /dev/null {hmm_out} {aln_path}"
    )
    logger.debug(f"hmmbuild command: {command}")
    result = subprocess.run(
        shlex.split(command),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"hmmbuild failed (exit {result.returncode}): {result.stderr.strip()}"
        )


def read_hmm_header(hmm_path: str) -> "dict[str, str]":
    """Return the header fields (``NAME``, ``LENG``, ``NSEQ``, ...) of a profile."""
    header: dict[str, str] = {}
    with open(hmm_path) as handle:
        for line in handle:
            if line.startswith("HMM "):
                break
            key, _, value = line.rstrip("\n").partition(" ")
            if key and value.strip():
                header.setdefault(key, value.strip())
    return header


def build_profile(
    sequences: "dict[str, str]",
    cluster: ClusterSpec,
    workdir: str,
    settings: BuilderSettings,
    threads: int = 1,
) -> "BuiltProfile | None":
    """Align, curate and build one profile for ``cluster``.

    Parameters
    ----------
    sequences : dict[str, str]
        ``protein_id -> sequence`` index covering the cluster members.
    cluster : ClusterSpec
        The cluster to build.
    workdir : str
        Directory for the intermediates (created by the caller).
    settings : BuilderSettings
        Clustering/curation/build settings.
    threads : int
        Threads for ``mafft``/``hmmbuild``.

    Returns
    -------
    BuiltProfile | None
        ``None`` when the cluster is too small, or when even the uncurated
        alignment is below the configured minimums.
    """
    # A singleton cluster is only buildable when the user asked for singletons.
    if settings.keep_singletons and len(cluster.members) == 1:
        threshold = 1
    else:
        threshold = settings.min_profile_seqs
    if len(cluster.members) < threshold:
        logger.debug(
            f"skipping {cluster.profile_id}: {len(cluster.members)} member(s) "
            f"< {threshold}"
        )
        return None

    member_faa = os.path.join(workdir, f"{cluster.profile_id}.faa")
    written = write_member_fasta(sequences, cluster.members, member_faa)
    if written < threshold:
        logger.debug(f"skipping {cluster.profile_id}: only {written} sequence(s) found")
        return None

    aln_path = os.path.join(workdir, f"{cluster.profile_id}.aln.fa")
    if written == 1:
        # A single sequence needs no alignment; hmmbuild accepts it as-is.
        shutil.copyfile(member_faa, aln_path)
        curated, removed_seqs, removed_cols = aln_path, 0, 0
    else:
        run_mafft(member_faa, aln_path, threads)
        if settings.msa_curation:
            curated, removed_seqs, removed_cols = curate_alignment(
                aln_path,
                os.path.join(workdir, f"{cluster.profile_id}.cialign"),
                settings,
                retain_id=cluster.representative,
            )
        else:
            curated, removed_seqs, removed_cols = aln_path, 0, 0

    n_seqs, n_cols = alignment_dimensions(curated)
    if n_seqs < threshold or n_cols < settings.min_alignment_columns:
        logger.debug(
            f"skipping {cluster.profile_id}: alignment has {n_seqs} sequence(s) / "
            f"{n_cols} column(s)"
        )
        return None

    hmm_path = os.path.join(workdir, f"{cluster.profile_id}.hmm")
    run_hmmbuild(curated, hmm_path, cluster.profile_id, threads)
    header = read_hmm_header(hmm_path)
    return BuiltProfile(
        profile_id=cluster.profile_id,
        hmm_path=hmm_path,
        representative=cluster.representative,
        members=list(cluster.members),
        n_seqs=int(header.get("NSEQ", n_seqs)),
        length=int(header.get("LENG", n_cols)),
        curated=curated != aln_path,
        removed_sequences=removed_seqs,
        removed_columns=removed_cols,
    )


def concat_profiles(hmm_paths: "list[str]", out_hmm: str) -> int:
    """Concatenate single profiles into one database file.

    Parameters
    ----------
    hmm_paths : list[str]
        Profile files, in the order they should appear.
    out_hmm : str
        Destination ``.hmm``.

    Returns
    -------
    int
        Number of profiles written.
    """
    with open(out_hmm, "w") as out:
        for path in hmm_paths:
            with open(path) as handle:
                out.write(handle.read())
    logger.debug(f"wrote {len(hmm_paths)} profile(s) to {out_hmm}")
    return len(hmm_paths)


def missing_binaries(settings: BuilderSettings) -> "list[str]":
    """Return the builder binaries that are not on ``PATH``.

    ``CIAlign`` is only required when curation is enabled.
    """
    required = [MAFFT_BINARY, HMMBUILD_BINARY]
    if settings.msa_curation:
        required.append(CIALIGN_BINARY)
    return [name for name in required if shutil.which(name) is None]
