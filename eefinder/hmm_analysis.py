"""Profile-HMM search (HMMER3) as an alternative to the translated BLAST search.

``hmmsearch`` scores protein queries against profile HMMs, which is measurably
more sensitive than pairwise search for the remote homologies endogenous
elements usually present. The engine is selected with ``screening -md hmmer``.

Two facts shape this module:

* ``hmmsearch`` needs **protein** queries, so the nucleotide contigs are
  translated first -- by default in all six frames
  (:func:`eefinder.translation.translate_sixframe`), since EEs are frequently
  pseudogenised and an ORF predictor would miss them.
* every downstream step of the pipeline consumes a ``blastx``-shaped ``outfmt 6``
  table, so :func:`domtbl_to_aa_table` rewrites ``--domtblout`` into those 12
  columns and :func:`eefinder.translation.traceback` maps the amino-acid
  coordinates back to contig nucleotides. Nothing after
  :class:`~eefinder.filter_table.FilterTable` changes.

Column mapping from ``--domtblout`` (``hmmsearch``: *target* = our translated
query, *query* = the profile):

===============  =======================================================
outfmt 6         source
===============  =======================================================
``qseqid``       ``target name`` -> contig via the coordinates TSV
``sseqid``       ``query name`` (the profile id; joins metadata ``Accession``)
``pident``       ``acc`` x 100 -- the mean posterior probability, **not**
                 percent identity (a profile hit has no subject sequence to
                 be identical to)
``length``       ``ali to - ali from + 1`` (aa, as in ``blastx``)
``mismatch``     0 (undefined for a profile hit)
``gapopen``      0 (undefined for a profile hit)
``qstart/qend``  ``ali from``/``ali to``, traced back to nucleotides
``sstart/send``  ``hmm from``/``hmm to`` (profile coordinates)
``evalue``       ``i-Evalue`` (the independent per-domain E-value)
``bitscore``     the per-domain ``score``
===============  =======================================================

E-values from ``hmmsearch`` scale with the size of the searched sequence set,
which varies with genome size, so ``-Z``/``--domZ`` are fixed to
:data:`DEFAULT_Z` -- significance is then comparable across runs and genomes.
"""

from __future__ import annotations

import math
import os
import shlex
import subprocess
from dataclasses import dataclass

from eefinder import translation
from eefinder.filter_table import OUTFMT6_COLUMNS
from eefinder.log import logger

#: HMMER binaries used by this module (from the ``hmmer`` conda package).
HMMSEARCH_BINARY = "hmmsearch"
HMMPRESS_BINARY = "hmmpress"
HMMBUILD_BINARY = "hmmbuild"

#: E-value cutoff, mirroring :data:`eefinder.similarity_analysis.EVALUE_CUTOFF`.
EVALUE_CUTOFF = 0.00001

#: Fixed database size for ``-Z``/``--domZ`` so E-values do not depend on how
#: many translated segments a particular genome happens to produce.
DEFAULT_Z = 1_000_000

#: Suffixes of the binary files ``hmmpress`` writes next to a ``.hmm``.
PRESS_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")


@dataclass
class DomainHit:
    """One domain hit parsed from a ``hmmsearch --domtblout`` row."""

    target: str  #: translated query segment id (``target name``)
    profile: str  #: profile id (``query name``)
    profile_length: int  #: ``qlen``
    evalue: float  #: ``i-Evalue``
    bitscore: float  #: per-domain ``score``
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    accuracy: float  #: ``acc``, the mean posterior probability

    @property
    def coverage(self) -> float:
        """Fraction of the profile covered by the alignment."""
        if self.profile_length <= 0:
            return 0.0
        return (self.hmm_to - self.hmm_from + 1) / self.profile_length


def parse_domtbl(domtbl_path: str):
    """Yield :class:`DomainHit` records from a ``--domtblout`` file.

    Comment lines (``#``) and malformed rows are skipped, so a truncated file
    (e.g. an interrupted search) degrades to fewer hits rather than an error.

    Parameters
    ----------
    domtbl_path : str
        Path to a ``hmmsearch --domtblout`` table.

    Yields
    ------
    DomainHit
    """
    with open(domtbl_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(None, 22)
            if len(fields) < 22:
                logger.debug(f"skipping malformed domtbl row: {line.rstrip()}")
                continue
            try:
                yield DomainHit(
                    target=fields[0],
                    profile=fields[3],
                    profile_length=int(fields[5]),
                    evalue=float(fields[12]),
                    bitscore=float(fields[13]),
                    hmm_from=int(fields[15]),
                    hmm_to=int(fields[16]),
                    ali_from=int(fields[17]),
                    ali_to=int(fields[18]),
                    accuracy=float(fields[21]),
                )
            except ValueError:
                logger.debug(f"skipping unparsable domtbl row: {line.rstrip()}")


def domtbl_to_aa_table(
    domtbl_path: str,
    out_path: str,
    min_coverage: float = 0.0,
    append: bool = False,
) -> int:
    """Convert a ``--domtblout`` table to ``outfmt 6`` with amino-acid coordinates.

    The output is the input :func:`eefinder.translation.traceback` expects: 12
    headerless :data:`~eefinder.filter_table.OUTFMT6_COLUMNS` columns whose
    ``qseqid`` is a translated-segment id and whose ``qstart``/``qend`` are
    amino-acid positions.

    Parameters
    ----------
    domtbl_path : str
        Table written by :func:`run_hmmsearch`.
    out_path : str
        Destination tabular file.
    min_coverage : float
        Drop hits covering less than this fraction of the profile (``0.0`` keeps
        everything). Short partial-domain hits are the dominant spurious-hit
        class in profile search.
    append : bool
        Append instead of truncating, so several databases can be merged into
        one table before the traceback.

    Returns
    -------
    int
        Number of hits written.
    """
    written = dropped = 0
    with open(out_path, "a" if append else "w") as out:
        for hit in parse_domtbl(domtbl_path):
            if min_coverage and hit.coverage < min_coverage:
                dropped += 1
                continue
            row = [
                hit.target,
                hit.profile,
                f"{hit.accuracy * 100:.1f}",
                str(hit.ali_to - hit.ali_from + 1),
                "0",
                "0",
                str(hit.ali_from),
                str(hit.ali_to),
                str(hit.hmm_from),
                str(hit.hmm_to),
                f"{hit.evalue:.3g}",
                f"{hit.bitscore:.1f}",
            ]
            out.write("\t".join(row) + "\n")
            written += 1
    if dropped:
        logger.debug(
            f"{dropped} hit(s) dropped below --hmm_min_coverage {min_coverage} "
            f"from {domtbl_path}"
        )
    logger.debug(f"{written} hit(s) written to {out_path} from {domtbl_path}")
    return written


def is_pressed(hmm_path: str) -> bool:
    """Whether ``hmmpress`` binary files already exist next to ``hmm_path``."""
    return all(os.path.exists(f"{hmm_path}{suffix}") for suffix in PRESS_SUFFIXES)


def press_database(hmm_path: str) -> None:
    """Index a profile database with ``hmmpress`` (overwriting any old index)."""
    command = f"{HMMPRESS_BINARY} -f {hmm_path}"
    logger.debug(f"hmmpress command: {command}")
    subprocess.run(
        shlex.split(command),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )


def has_ga_thresholds(hmm_path: str) -> bool:
    """Whether **every** profile in ``hmm_path`` carries a ``GA`` line.

    ``hmmsearch --cut_ga`` fails on a database where any model lacks its
    gathering threshold, so this is checked before the flag is used (the
    calibration of ``proposal.md`` §6.2 is what writes those lines).
    """
    models = thresholds = 0
    with open(hmm_path) as handle:
        for line in handle:
            if line.startswith("NAME "):
                models += 1
            elif line.startswith("GA "):
                thresholds += 1
    return models > 0 and models == thresholds


def run_hmmsearch(
    hmm_database: str,
    query_faa: str,
    domtbl_out: str,
    threads: int = 1,
    evalue: float = EVALUE_CUTOFF,
    z: int = DEFAULT_Z,
    use_ga: bool = False,
    sensitive: bool = False,
) -> None:
    """Search ``query_faa`` against ``hmm_database``, writing ``--domtblout``.

    Parameters
    ----------
    hmm_database : str
        Profile database (``.hmm``); its ``hmmpress`` index is used when present.
    query_faa : str
        Protein FASTA of translated queries.
    domtbl_out : str
        Destination ``--domtblout`` table.
    threads : int
        ``--cpu``.
    evalue : float
        ``-E``/``--domE`` significance cutoff (ignored when ``use_ga``).
    z : int
        Fixed ``-Z``/``--domZ`` database size (see :data:`DEFAULT_Z`).
    use_ga : bool
        Use each profile's ``GA`` threshold (``--cut_ga``) instead of ``evalue``.
    sensitive : bool
        Pass ``--max`` (disables the heuristic filters: much more sensitive on
        diverged elements, roughly an order of magnitude slower).

    Raises
    ------
    RuntimeError
        If ``hmmsearch`` exits non-zero (its stderr is included, unlike the
        DIAMOND path, so failures are never silent).
    ValueError
        If ``use_ga`` is requested for a database without ``GA`` lines.
    """
    if use_ga and not has_ga_thresholds(hmm_database):
        raise ValueError(
            f"{hmm_database} has profiles without GA thresholds: --cut_ga cannot "
            "be used. Calibrate the database or drop --hmm_use_ga."
        )
    cutoff = "--cut_ga" if use_ga else f"-E {evalue} --domE {evalue}"
    command = (
        f"{HMMSEARCH_BINARY} --domtblout {domtbl_out} --noali "
        f"--cpu {int(threads)} {cutoff} -Z {int(z)} --domZ {int(z)} "
        f"{'--max ' if sensitive else ''}"
        f"{hmm_database} {query_faa}"
    )
    logger.debug(f"hmmsearch command: {command}")
    result = subprocess.run(
        shlex.split(command),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"hmmsearch failed (exit {result.returncode}): {result.stderr.strip()}"
        )


class HmmSearch:
    """Search a translated query against one or more profile databases.

    The HMM counterpart of :class:`~eefinder.similarity_analysis.SimilaritySearch`
    and, like it, runs on instantiation and writes a ``blastx``-shaped table to
    ``{query_file}.blastx`` (the name is kept so every downstream step of
    ``main.py`` is untouched).

    With several databases, each is searched separately and the hit tables are
    concatenated before the coordinate traceback -- which is what lets each
    database keep its own thresholds and calibration state.

    Parameters
    ----------
    query_file : str
        Nucleotide FASTA (contigs, or extracted candidate regions).
    hmm_databases : str | list[str]
        One or more pressed profile databases.
    threads : int
        Threads for the translation and search.
    translation_method : str
        How query proteins are obtained: ``"sixframe"`` (default) or a
        prediction method (``"gv"``/``"rv"``/``"gv-rv"``).
    evalue, z, min_coverage, use_ga, sensitive
        Search parameters; see :func:`run_hmmsearch` and
        :func:`domtbl_to_aa_table`.
    out_table : str, optional
        Override the output path (defaults to ``{query_file}.blastx``).
    """

    def __init__(
        self,
        query_file: str,
        hmm_databases: "str | list[str]",
        threads: int = 1,
        translation_method: str = "sixframe",
        evalue: float = EVALUE_CUTOFF,
        z: int = DEFAULT_Z,
        min_coverage: float = 0.0,
        use_ga: bool = False,
        sensitive: bool = False,
        out_table: "str | None" = None,
    ) -> None:
        self.query_file = query_file
        self.hmm_databases = (
            [hmm_databases] if isinstance(hmm_databases, str) else list(hmm_databases)
        )
        self.threads = threads
        self.translation_method = (
            "sixframe" if translation_method == "default" else translation_method
        )
        self.evalue = evalue
        self.z = z
        self.min_coverage = min_coverage
        self.use_ga = use_ga
        self.sensitive = sensitive
        self.out_table = out_table or f"{query_file}.blastx"

        self.hmm_search()

    def hmm_search(self) -> None:
        """Translate, search every database, then trace hits back to nucleotides."""
        protein_faa, coords_tsv = translation.predict_and_cluster(
            self.query_file, self.translation_method, self.threads
        )
        aa_table = f"{self.query_file}.hmm.aa.tsv"
        # Truncate any table left by a previous run before appending per database.
        open(aa_table, "w").close()

        total = 0
        for index, database in enumerate(self.hmm_databases):
            domtbl = f"{self.query_file}.hmm.{index}.domtbl"
            logger.debug(f"HmmSearch: {protein_faa} vs {database}")
            run_hmmsearch(
                database,
                protein_faa,
                domtbl,
                threads=self.threads,
                evalue=self.evalue,
                z=self.z,
                use_ga=self.use_ga,
                sensitive=self.sensitive,
            )
            total += domtbl_to_aa_table(
                domtbl, aa_table, min_coverage=self.min_coverage, append=True
            )

        logger.debug(
            f"HmmSearch: {total} hit(s) from {len(self.hmm_databases)} database(s); "
            f"tracing coordinates back to {self.out_table}"
        )
        translation.traceback(aa_table, coords_tsv, self.out_table)


def empty_outfmt6(path: str) -> None:
    """Write an empty ``outfmt 6`` table (helper for no-hit edge cases)."""
    with open(path, "w"):
        pass


#: Re-exported for callers that build their own tables (keeps the column order
#: in one place).
__all__ = [
    "CALIBRATION_EVALUE",
    "DEFAULT_Z",
    "DomainHit",
    "EVALUE_CUTOFF",
    "HmmSearch",
    "OUTFMT6_COLUMNS",
    "calibrate_database",
    "collect_host_scores",
    "domtbl_to_aa_table",
    "empty_outfmt6",
    "has_ga_thresholds",
    "is_pressed",
    "parse_domtbl",
    "press_database",
    "run_hmmsearch",
    "score_for_evalue",
    "write_ga_thresholds",
]


# --------------------------------------------------------------------------- #
# Host calibration (writes per-profile GA thresholds)
# --------------------------------------------------------------------------- #
#: Permissive E-value used when scoring profiles against a host proteome: the
#: point is to observe how high each profile *can* score on host proteins, so a
#: strict cutoff would hide exactly the borderline hits the threshold protects
#: against.
CALIBRATION_EVALUE = 10.0


def collect_host_scores(
    hmm_database: str,
    host_proteome: str,
    workdir: str,
    threads: int = 1,
    evalue: float = CALIBRATION_EVALUE,
    quantile: float = 1.0,
) -> "dict[str, float]":
    """Score every profile against a host proteome and summarise per profile.

    Parameters
    ----------
    hmm_database : str
        Profile database to calibrate.
    host_proteome : str
        Host protein FASTA (the curated ``-bt`` bait set is a safer choice than a
        whole proteome, which may already contain endogenised elements).
    workdir : str
        Directory for the intermediate ``--domtblout``.
    threads : int
        ``--cpu``.
    evalue : float
        Permissive significance cutoff (see :data:`CALIBRATION_EVALUE`).
    quantile : float
        Which point of each profile's host-score distribution to take: ``1.0``
        (the default) is the maximum, ``0.99`` ignores a single outlier -- useful
        when the proteome may itself contain an endogenised element.

    Returns
    -------
    dict[str, float]
        ``profile id -> host score`` for the profiles that hit anything.
    """
    os.makedirs(workdir, exist_ok=True)
    domtbl = os.path.join(workdir, "calibration.domtbl")
    run_hmmsearch(
        hmm_database,
        host_proteome,
        domtbl,
        threads=threads,
        evalue=evalue,
    )
    scores: dict[str, list[float]] = {}
    for hit in parse_domtbl(domtbl):
        scores.setdefault(hit.profile, []).append(hit.bitscore)

    summarised: dict[str, float] = {}
    for profile, values in scores.items():
        values.sort()
        if quantile >= 1.0:
            summarised[profile] = values[-1]
        else:
            index = min(len(values) - 1, int(round(quantile * (len(values) - 1))))
            summarised[profile] = values[index]
    logger.debug(
        f"calibration: {len(summarised)} of the profiles in {hmm_database} scored "
        f"on {host_proteome}"
    )
    return summarised


def score_for_evalue(tau: float, lambda_: float, evalue: float, z: float) -> float:
    """Bit score at which a profile's Forward E-value equals ``evalue``.

    HMMER's Forward tail is exponential -- ``E = Z * exp(-lambda * (S - tau))``
    with ``tau``/``lambda`` from the profile's ``STATS LOCAL FORWARD`` line --
    so the equivalent score is ``tau + ln(Z / E) / lambda``.

    This is what a calibrated threshold falls back to for a profile that never
    scored on the host: ``--cut_ga`` *replaces* the E-value cutoff rather than
    adding to it, so a zero floor would make a calibrated database **more**
    permissive than an uncalibrated one.
    """
    if lambda_ <= 0:
        return 0.0
    return tau + math.log(z / evalue) / lambda_


def _forward_stats(lines: "list[str]") -> "tuple[float, float] | None":
    """Return ``(tau, lambda)`` from a profile's ``STATS LOCAL FORWARD`` line."""
    for line in lines:
        if line.startswith("STATS LOCAL FORWARD"):
            fields = line.split()
            try:
                return float(fields[3]), float(fields[4])
            except (IndexError, ValueError):
                return None
    return None


def write_ga_thresholds(
    hmm_database: str,
    scores: "dict[str, float]",
    out_path: str,
    margin: float = 0.1,
    evalue: float = EVALUE_CUTOFF,
    z: float = DEFAULT_Z,
) -> int:
    """Write a ``GA`` gathering threshold into every profile of a database.

    A profile's threshold is the score it achieved on the host proteome plus
    ``margin``, so at screening time (``hmmsearch --cut_ga``) a hit only counts
    when it beats what that profile can achieve on the host's own genes.

    Profiles that never hit the host fall back to the score equivalent to
    ``evalue`` at ``z`` (:func:`score_for_evalue`), so a calibrated database is
    never *more* permissive than an uncalibrated one -- ``--cut_ga`` replaces the
    E-value cutoff rather than adding to it.

    Both the sequence and domain thresholds are written (HMMER requires both).

    Parameters
    ----------
    hmm_database : str
        Input profile database.
    scores : dict[str, float]
        ``profile -> host score`` from :func:`collect_host_scores`.
    out_path : str
        Destination database (may be the same path as the input).
    margin : float
        Relative margin added to each observed host score.
    evalue, z : float
        The significance the floor corresponds to; pass the values the screening
        run will use so an uncalibrated-equivalent profile behaves identically.

    Returns
    -------
    int
        Number of profiles written.
    """
    with open(hmm_database) as handle:
        text = handle.read()

    written = 0
    out_blocks = []
    for block in text.split("//\n"):
        if not block.strip():
            continue
        lines = block.splitlines(keepends=True)
        name = ""
        insert_at = len(lines)
        for index, line in enumerate(lines):
            if line.startswith("NAME "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith(("GA ", "TC ", "NC ")):
                continue  # replaced below
            elif line.startswith(("STATS ", "HMM ")) and insert_at == len(lines):
                insert_at = index
        lines = [line for line in lines if not line.startswith(("GA ", "TC ", "NC "))]
        stats = _forward_stats(lines)
        floor = score_for_evalue(stats[0], stats[1], evalue, z) if stats else 0.0
        threshold = max(floor, scores.get(name, 0.0) * (1 + margin))
        lines.insert(insert_at, f"GA    {threshold:.2f} {threshold:.2f};\n")
        out_blocks.append("".join(lines))
        written += 1

    with open(out_path, "w") as out:
        out.write("//\n".join(out_blocks) + "//\n")
    logger.debug(f"calibration: wrote GA thresholds for {written} profile(s)")
    return written


def calibrate_database(
    hmm_database: str,
    host_proteome: str,
    workdir: str,
    margin: float = 0.1,
    quantile: float = 1.0,
    threads: int = 1,
    evalue: float = EVALUE_CUTOFF,
    z: float = DEFAULT_Z,
) -> "tuple[int, int]":
    """Calibrate a profile database against a host proteome, in place.

    Runs :func:`collect_host_scores`, writes the thresholds with
    :func:`write_ga_thresholds` and re-indexes the database with ``hmmpress``, so
    it can then be searched with ``--hmm_use_ga``.

    Returns
    -------
    tuple[int, int]
        ``(profiles written, profiles that scored on the host)``.
    """
    scores = collect_host_scores(
        hmm_database, host_proteome, workdir, threads=threads, quantile=quantile
    )
    calibrated = os.path.join(workdir, "calibrated.hmm")
    written = write_ga_thresholds(
        hmm_database, scores, calibrated, margin=margin, evalue=evalue, z=z
    )
    os.replace(calibrated, hmm_database)
    press_database(hmm_database)
    return written, len(scores)
