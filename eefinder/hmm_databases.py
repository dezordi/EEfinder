"""Assemble the profile-HMM databases the ``hmmer`` engine screens with.

Two databases are produced by two ``get-databases`` subcommands:

* ``virus-hmm`` -- the viral profile database (``screening -dbh``) plus its
  metadata CSV (``-mth``). The ``ncbi-refseq`` source builds it from the protein
  FASTA + metadata CSV ``get-databases virus`` already produces, which is what
  makes its taxonomy exact rather than inferred: **proteins are clustered inside
  one viral taxon at a time** (``--cluster-level family|genus``), so a profile
  cannot span families and its ``Family``/``Genus`` are known by construction.
* ``host-hmm`` -- the host profile database (``screening -bth``) the symmetric EE
  filter compares against. Hypothetical/uncharacterised products are excluded
  (an uninformative host profile is unauditable: when it wins a comparison the
  user cannot check the claim) and so are transposon/retroelement-like products,
  which would otherwise out-score viral profiles over genuine
  retroelement-derived elements and delete true positives.

Each cluster gets **one** protein name by majority vote
(:func:`eefinder.normalization.consensus_protein_name`) on the names as they
already appear in the metadata CSV -- i.e. after the bundled canonical-name map
has been applied, so only its blind spots reach the vote. Clusters whose vote was
close are listed in ``{prefix}.protein_votes.tsv``, which is the worklist for
extending ``eefinder/data/viral_proteins.tsv``.

Both databases record their :class:`~eefinder.hmm_builder.BuilderSettings` in
their log; the symmetric filter is only valid when the two agree.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from eefinder import __version__
from eefinder.get_databases import (
    GetDatabases,
    METADATA_COLUMNS,
    UNINFORMATIVE_PRODUCTS,
    cluster_identical_proteins,
    molecule_type_for_family,
    parse_protein_header,
)
from eefinder.hmm_analysis import press_database
from eefinder.hmm_builder import (
    BuilderSettings,
    BuiltProfile,
    ClusterSpec,
    build_profile,
    concat_profiles,
    hmmer_version,
    make_profile_id,
    slugify,
)
from eefinder.hmm_sources import (
    NEORDRP_DEFAULT_RELEASE,
    ProfileRecord,
    RVDB_DEFAULT_RELEASE,
    VOGDB_PROFILE_ARCHIVE,
    VogdbFilters,
    clean_vogdb_description,
    download_neordrp,
    download_rvdb,
    download_vogdb,
    extract_profiles,
    extract_profiles_from_hmm,
    filter_vogs,
    load_neordrp,
    load_rvdb,
    load_vogdb,
    parse_hosts,
    parse_species,
    profile_length,
    vogdb_release,
)
from eefinder.log import logger
from eefinder.normalization import (
    UNKNOWN_PROTEIN,
    consensus_protein_name,
    format_votes,
    is_transposon_like,
    standardize_protein,
)
from eefinder.translation import cluster_homologs
from eefinder.utils import (
    HmmBuildInfo,
    ProfileCounts,
    StepInfo,
    check_outdir,
)

#: Viral profile sources. ``ncbi-refseq`` builds profiles from a RefSeq download;
#: ``vogdb``, ``rvdb`` and ``neordrp`` download a prebuilt release and map it onto
#: the same metadata (see :mod:`eefinder.hmm_sources`).
VIRAL_SOURCES = ("ncbi-refseq", "vogdb", "rvdb", "neordrp")

#: Namespace prefix per source, so profile ids stay unique (and say where they
#: came from) when several databases are searched together.
SOURCE_NAMESPACES = {
    "ncbi-refseq": "NCBIREFSEQ",
    "vogdb": "VOGDB",
    "rvdb": "RVDB",
    "neordrp": "NEORDRP",
    "host": "HOST",
}

#: Taxonomic levels a viral build can be scoped to.
CLUSTER_LEVELS = ("family", "genus")

#: Metadata CSV columns of a viral profile database: the columns ``screening``
#: already consumes, plus provenance/quality columns **appended after** ``Host``
#: (``get_taxonomy.GetFinalTaxonomy`` reads the first columns positionally).
HMM_METADATA_COLUMNS = METADATA_COLUMNS + [
    "Protein_votes",
    "Protein_agreement",
    "Profile_seqs",
    "Profile_length",
    "Source",
    "LCA_rank",
]

#: Metadata CSV columns of a host profile database (no viral taxonomy).
HOST_METADATA_COLUMNS = [
    "Accession",
    "Protein",
    "Protein_votes",
    "Protein_agreement",
    "Profile_seqs",
    "Profile_length",
    "Source",
]

#: Taxon values that carry no taxonomic information (``nan``/``none`` cover the
#: string forms a missing pandas value takes).
_UNCLASSIFIED = {"", "unclassified", "unknown", "nan", "none", "undefined"}

#: Agreement below which a cluster is listed in ``{prefix}.protein_votes.tsv``.
LOW_AGREEMENT = 0.6


def is_unclassified(value: object) -> bool:
    """Whether a taxonomy value is missing/unclassified."""
    return str(value).strip().lower() in _UNCLASSIFIED


def read_sequences(fasta_path: str) -> "dict[str, str]":
    """Index a protein FASTA as ``id -> sequence`` (first header token as id)."""
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}


def plurality(values: "list[str]", minimum: float = 0.9, fallback: str = "") -> str:
    """Return the most frequent value when it holds ``minimum`` of the votes.

    Used for the metadata columns that are not the protein name: a cluster spans
    several species, so ``Species`` is only reported when one dominates.

    Parameters
    ----------
    values : list[str]
        Member values (missing/unclassified entries are ignored).
    minimum : float
        Minimum share the top value must hold.
    fallback : str
        Returned when no value reaches ``minimum``.

    Returns
    -------
    str
    """
    kept = [str(value).strip() for value in values if not is_unclassified(value)]
    if not kept:
        return fallback
    counts: dict[str, int] = {}
    for value in kept:
        counts[value] = counts.get(value, 0) + 1
    top, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
    return top if count / len(kept) >= minimum else fallback


def is_uninformative(product: str) -> bool:
    """Whether a product name is a hypothetical/uncharacterised placeholder."""
    lowered = str(product).lower()
    return any(term in lowered for term in UNINFORMATIVE_PRODUCTS) or (
        lowered.strip() == UNKNOWN_PROTEIN.lower()
    )


class HmmDatabaseBuilder:
    """Shared machinery for the viral and host profile-database builds.

    Not a step class: subclasses implement :meth:`build` and run it from their
    own ``__init__``, following the pipeline's side-effect convention.
    """

    #: Source name used for the id namespace and the ``Source`` column.
    source = "ncbi-refseq"

    def __init__(
        self,
        outdir: str,
        prefix: str,
        settings: BuilderSettings,
        threads: int = 1,
    ) -> None:
        self.outdir = check_outdir(outdir)
        self.prefix = prefix
        self.settings = settings
        self.threads = threads
        self.workdir = check_outdir(f"{self.outdir}/{self.prefix}_build")
        self.profiles: list[BuiltProfile] = []
        self.rows: list[dict] = []
        self.votes: list[dict] = []
        self.clusters_formed = 0
        self.clusters_skipped = 0

    # -- helpers shared by both builds -------------------------------------- #
    def _cluster(self, fasta_path: str, tag: str) -> "dict[str, list[str]]":
        """Cluster a protein FASTA into homologous groups with ``cd-hit``."""
        clustered = os.path.join(self.workdir, f"{tag}.nr.faa")
        clusters = cluster_homologs(
            fasta_path, clustered, self.settings.cluster_identity, self.threads
        )
        logger.debug(
            f"{tag}: {len(clusters)} cluster(s) at -c {self.settings.cluster_identity}"
        )
        return clusters

    def _build_cluster(
        self, sequences: "dict[str, str]", spec: ClusterSpec
    ) -> "BuiltProfile | None":
        """Build one profile, counting formed vs skipped clusters."""
        self.clusters_formed += 1
        profile = build_profile(
            sequences, spec, self.workdir, self.settings, self.threads
        )
        if profile is None:
            self.clusters_skipped += 1
            return None
        self.profiles.append(profile)
        return profile

    def _record_votes(self, profile_id: str, consensus, taxon: str = "") -> None:
        """Remember a low-agreement vote for ``{prefix}.protein_votes.tsv``."""
        if consensus.agreement < LOW_AGREEMENT and consensus.votes:
            self.votes.append(
                {
                    "Profile": profile_id,
                    "Taxon": taxon,
                    "Protein": consensus.name,
                    "Agreement": f"{consensus.agreement:.2f}",
                    "Votes": format_votes(consensus.votes),
                }
            )

    def _write_outputs(self, columns: "list[str]", metadata_only: bool = False) -> str:
        """Concatenate the profiles, press them, and write the sidecar tables.

        With ``metadata_only`` the profiles were never fetched, so only the CSV
        and the sidecar tables are written -- the cheap way to size a prebuilt
        database before committing the disk to it.
        """
        hmm_path = f"{self.outdir}/{self.prefix}.hmm"
        if metadata_only:
            logger.info(
                "--metadata-only: reporting what would be kept; no profile "
                "database was written."
            )
        else:
            concat_profiles([profile.hmm_path for profile in self.profiles], hmm_path)
            if self.profiles:
                press_database(hmm_path)
            else:
                logger.warning(
                    f"No profile survived the filters; {hmm_path} is empty and was "
                    "not indexed."
                )

        csv_path = f"{self.outdir}/{self.prefix}.csv"
        pd.DataFrame(self.rows, columns=columns).to_csv(csv_path, index=False)
        logger.info(f"Wrote {csv_path} ({len(self.rows)} profile(s))")

        clusters_path = f"{self.outdir}/{self.prefix}.clusters.tsv"
        with open(clusters_path, "w") as out:
            out.write("Profile\tRepresentative\tMembers\n")
            for profile in self.profiles:
                out.write(
                    f"{profile.profile_id}\t{profile.representative}\t"
                    f"{','.join(profile.members)}\n"
                )

        votes_path = f"{self.outdir}/{self.prefix}.protein_votes.tsv"
        vote_columns = ["Profile", "Taxon", "Protein", "Agreement", "Votes"]
        pd.DataFrame(self.votes, columns=vote_columns).to_csv(
            votes_path, sep="\t", index=False
        )
        if self.votes:
            logger.info(
                f"{len(self.votes)} cluster(s) with protein-name agreement below "
                f"{LOW_AGREEMENT} listed in {votes_path}"
            )
        return hmm_path

    def _write_log(
        self,
        arguments: dict,
        counts: ProfileCounts,
        start_time: float,
        steps: "list[StepInfo]",
    ) -> None:
        """Write the ``{prefix}.log`` JSON build summary."""
        info = HmmBuildInfo.from_run(
            eefinder_version=__version__,
            arguments=arguments,
            builder_settings=self.settings.to_dict(),
            profile_counts=counts,
            start_time=start_time,
            end_time=time.time(),
            steps_information=steps,
        )
        log_path = f"{self.outdir}/{self.prefix}.log"
        with open(log_path, "w") as json_out:
            json.dump(asdict(info), json_out, indent=4)
        logger.info(f"Wrote {log_path}")


class GetViralHmmDatabase(HmmDatabaseBuilder):
    """Build the viral profile database from one or more sources.

    Runs on instantiation. Writes ``{outdir}/{prefix}.hmm`` (+ ``hmmpress``
    index), ``{prefix}.csv``, ``{prefix}.clusters.tsv``,
    ``{prefix}.protein_votes.tsv`` and ``{prefix}.log``.

    Sources are merged into one pressed database in the order given; profile ids
    are namespaced per source (``NCBIREFSEQ__…``, ``VOGDB__…``), so they never
    collide and every row says where it came from.

    Parameters
    ----------
    outdir : str
        Output directory (created if missing).
    prefix : str
        Basename of the outputs.
    sources : list[str]
        Any of :data:`VIRAL_SOURCES`, in merge order.
    taxon : str
        NCBI taxon to build from: downloaded for ``ncbi-refseq``, matched against
        the LCA lineage for ``vogdb``.
    proteins, metadata : str, optional
        An existing protein FASTA and metadata CSV (as written by
        ``get-databases virus``) for the ``ncbi-refseq`` source; when omitted they
        are downloaded for ``taxon``.
    cluster_level : str
        ``"family"`` or ``"genus"`` -- the taxon **within which** proteins are
        clustered, and therefore the level at which the profile taxonomy is exact.
    settings : BuilderSettings
        Clustering/curation/build settings.
    protein_vote : str
        Vote mode for the consensus protein name.
    protein_vote_min_fraction : float
        Minimum canonical-name share for ``canonical-first``.
    keep_unclassified : bool
        Keep proteins whose taxon is missing at ``cluster_level``.
    vogdb_filters : VogdbFilters, optional
        Filters applied to a VOGDB release (eukaryotic scope, virus-only
        stringency, minimum LCA rank, ...).
    release : str
        VOGDB release to pin (``"latest"`` by default).
    metadata_only : bool
        For prebuilt sources: download only the metadata tables and report how
        many profiles would survive, without fetching the profile archive. No
        ``.hmm`` is written.
    threads : int
        Threads for ``cd-hit``/``mafft``/``hmmbuild``.
    """

    source = "ncbi-refseq"

    def __init__(
        self,
        outdir: str,
        prefix: str = "virus_hmm",
        sources: "list[str] | None" = None,
        taxon: str = "10239",
        proteins: "str | None" = None,
        metadata: "str | None" = None,
        cluster_level: str = "family",
        settings: "BuilderSettings | None" = None,
        protein_vote: str = "canonical-first",
        protein_vote_min_fraction: float = 0.3,
        keep_unclassified: bool = False,
        vogdb_filters: "VogdbFilters | None" = None,
        release: str = "latest",
        metadata_only: bool = False,
        threads: int = 1,
    ) -> None:
        if cluster_level not in CLUSTER_LEVELS:
            raise ValueError(f"Unknown cluster level: {cluster_level!r}")
        self.sources = list(sources or [self.source])
        unknown = [s for s in self.sources if s not in VIRAL_SOURCES]
        if unknown:
            raise ValueError(f"Unknown source(s): {', '.join(unknown)}")
        if metadata_only and "ncbi-refseq" in self.sources:
            raise ValueError(
                "--metadata-only sizes a prebuilt source; it does not apply to "
                "ncbi-refseq, which builds its profiles locally."
            )
        settings = settings or BuilderSettings()
        settings.hmmer_version = settings.hmmer_version or hmmer_version()
        super().__init__(outdir, prefix, settings, threads)
        self.taxon = taxon
        self.proteins = proteins
        self.metadata = metadata
        self.cluster_level = cluster_level
        self.protein_vote = protein_vote
        self.protein_vote_min_fraction = protein_vote_min_fraction
        self.keep_unclassified = keep_unclassified
        self.vogdb_filters = vogdb_filters or VogdbFilters(taxon=taxon)
        self.release = release
        self.metadata_only = metadata_only
        self.proteins_in = 0
        self.excluded_uninformative = 0

        self.build()

    def build(self) -> None:
        """Run every requested source, then write the merged database."""
        start_time = time.time()
        steps: "list[StepInfo]" = []

        for source in self.sources:
            if source == "ncbi-refseq":
                self._build_refseq(steps)
            elif source == "vogdb":
                self._add_vogdb(steps)
            elif source == "rvdb":
                self._add_rvdb(steps)
            elif source == "neordrp":
                self._add_neordrp(steps)

        step_start = time.time()
        hmm_path = self._write_outputs(
            HMM_METADATA_COLUMNS, metadata_only=self.metadata_only
        )
        steps.append(
            StepInfo.from_times(
                "Write database",
                step_start,
                time.time(),
                f"Wrote {len(self.rows)} metadata row(s)"
                + (
                    " (--metadata-only: no profiles were fetched)."
                    if self.metadata_only
                    else f" and pressed {len(self.profiles)} profile(s) to {hmm_path}."
                ),
            )
        )

        counts = ProfileCounts(
            proteins_in=self.proteins_in,
            excluded_uninformative=self.excluded_uninformative,
            excluded_transposon_like=0,
            clustered_identical=0,
            clusters_formed=self.clusters_formed,
            clusters_skipped=self.clusters_skipped,
            profiles_built=len(self.rows),
        )
        logger.info(
            f"Profiles: {counts.profiles_built} in the database "
            f"(sources: {', '.join(self.sources)})"
        )
        self._write_log(
            {
                "sources": self.sources,
                "taxon": self.taxon,
                "outdir": self.outdir,
                "prefix": self.prefix,
                "cluster_level": self.cluster_level,
                "protein_vote": self.protein_vote,
                "protein_vote_min_fraction": self.protein_vote_min_fraction,
                "keep_unclassified": self.keep_unclassified,
                "vogdb_filters": asdict(self.vogdb_filters),
                "release": self.release,
                "metadata_only": self.metadata_only,
                "threads": self.threads,
            },
            counts,
            start_time,
            steps,
        )

    # -- ncbi-refseq -------------------------------------------------------- #
    def _acquire(self) -> "tuple[str, str]":
        """Return the protein FASTA and metadata CSV, downloading if needed."""
        if self.proteins and self.metadata:
            return self.proteins, self.metadata
        logger.info(f"Downloading RefSeq viral proteins for taxon '{self.taxon}'")
        download_prefix = f"{self.prefix}_proteins"
        GetDatabases(
            dataset="virus",
            taxon=self.taxon,
            outdir=self.outdir,
            prefix=download_prefix,
            exclude_uninformative=True,
            standardize_proteins=True,
            cluster=True,
            threads=self.threads,
        )
        return (
            f"{self.outdir}/{download_prefix}.fa",
            f"{self.outdir}/{download_prefix}.csv",
        )

    def _build_refseq(self, steps: "list[StepInfo]") -> None:
        """Cluster a RefSeq download per taxon and build a profile per cluster."""
        step_start = time.time()
        proteins, metadata = self._acquire()
        frame = pd.read_csv(metadata, dtype=str).fillna("")
        sequences = read_sequences(proteins)
        self.proteins_in += len(sequences)
        logger.info(
            f"Read {len(frame)} metadata row(s) and {len(sequences)} protein(s)"
        )
        steps.append(
            StepInfo.from_times(
                "Read source database",
                step_start,
                time.time(),
                f"Read {len(sequences)} protein(s) from {proteins} and "
                f"{len(frame)} metadata row(s) from {metadata}.",
            )
        )

        level_column = "Family" if self.cluster_level == "family" else "Genus"
        step_start = time.time()
        buckets = self._bucket(frame, level_column)
        logger.info(f"Clustering inside {len(buckets)} {self.cluster_level} bucket(s)")

        records = {row["Accession"]: row for row in frame.to_dict("records")}
        for taxon, accessions in sorted(buckets.items()):
            self._build_bucket(taxon, accessions, records, sequences)
        steps.append(
            StepInfo.from_times(
                "Build profiles",
                step_start,
                time.time(),
                f"Built {len(self.profiles)} profile(s) from "
                f"{self.clusters_formed} cluster(s) across {len(buckets)} "
                f"{self.cluster_level} bucket(s); {self.clusters_skipped} "
                "cluster(s) skipped as too small/short.",
            )
        )

    def _bucket(self, frame: pd.DataFrame, level_column: str) -> "dict[str, list[str]]":
        """Group accessions by their taxon at the configured level."""
        buckets: dict[str, list[str]] = {}
        for row in frame.to_dict("records"):
            taxon = str(row.get(level_column, "")).strip()
            if is_unclassified(taxon):
                if not self.keep_unclassified:
                    continue
                mol_type = str(row.get("Molecule_type", "")).strip() or "unknown"
                taxon = f"Unclassified-{mol_type}"
            buckets.setdefault(taxon, []).append(str(row["Accession"]))
        return buckets

    def _build_bucket(
        self,
        taxon: str,
        accessions: "list[str]",
        records: "dict[str, dict]",
        sequences: "dict[str, str]",
    ) -> None:
        """Cluster one taxon bucket and build a profile per cluster."""
        bucket_fasta = os.path.join(self.workdir, f"{taxon.replace('/', '_')}.faa")
        present = [acc for acc in accessions if acc in sequences]
        if not present:
            return
        with open(bucket_fasta, "w") as out:
            for accession in sorted(present):
                out.write(f">{accession}\n{sequences[accession]}\n")

        clusters = self._cluster(bucket_fasta, taxon.replace("/", "_"))
        counters: dict[str, int] = {}
        for representative, members in sorted(clusters.items()):
            member_rows = [records[m] for m in members if m in records]
            consensus = consensus_protein_name(
                [row.get("Protein", "") for row in member_rows],
                mode=self.protein_vote,
                min_fraction=self.protein_vote_min_fraction,
            )
            counters[consensus.name] = counters.get(consensus.name, 0) + 1
            profile_id = make_profile_id(
                SOURCE_NAMESPACES["ncbi-refseq"],
                taxon,
                consensus.name,
                counters[consensus.name],
            )
            profile = self._build_cluster(
                sequences,
                ClusterSpec(
                    profile_id=profile_id,
                    representative=representative,
                    members=members,
                ),
            )
            if profile is None:
                continue
            self._record_votes(profile_id, consensus, taxon)
            self.rows.append(self._metadata_row(profile, taxon, consensus, member_rows))

    def _metadata_row(
        self,
        profile: BuiltProfile,
        taxon: str,
        consensus,
        member_rows: "list[dict]",
    ) -> dict:
        """Aggregate the member metadata into one row for the profile."""
        if self.cluster_level == "family":
            family = taxon
            genus = plurality(
                [row.get("Genus", "") for row in member_rows], fallback="Unclassified"
            )
        else:
            genus = taxon
            family = plurality(
                [row.get("Family", "") for row in member_rows], fallback="Unclassified"
            )
        mol_type = molecule_type_for_family(family) or plurality(
            [row.get("Molecule_type", "") for row in member_rows], minimum=0.5
        )
        return {
            "Accession": profile.profile_id,
            "Species": plurality(
                [row.get("Species", "") for row in member_rows],
                fallback="Unclassified",
            ),
            "Genus": genus,
            "Family": family,
            "Molecule_type": mol_type,
            "Protein": consensus.name,
            "Host": plurality(
                [row.get("Host", "") for row in member_rows],
                minimum=0.5,
                fallback="Unknown",
            ),
            "Protein_votes": format_votes(consensus.votes),
            "Protein_agreement": f"{consensus.agreement:.2f}",
            "Profile_seqs": profile.n_seqs,
            "Profile_length": profile.length,
            "Source": "ncbi-refseq",
            "LCA_rank": self.cluster_level,
        }

    # -- vogdb -------------------------------------------------------------- #
    def _add_vogdb(self, steps: "list[StepInfo]") -> None:
        """Download a VOGDB release, filter it, and add the surviving profiles.

        Every filter runs on the small metadata tables, so ``--metadata-only``
        reports the surviving count without fetching the 554 MB profile archive.
        """
        step_start = time.time()
        workdir = check_outdir(f"{self.outdir}/vogdb")
        paths = download_vogdb(
            workdir, release=self.release, metadata_only=self.metadata_only
        )
        release = vogdb_release(workdir)
        records = load_vogdb(paths)
        hosts = parse_hosts(paths["vogdb.host.txt"])
        species = parse_species(paths["vogdb.species.txt"])
        steps.append(
            StepInfo.from_times(
                "Download VOGDB",
                step_start,
                time.time(),
                f"Downloaded VOGDB release {release} "
                f"({'metadata only' if self.metadata_only else 'with profiles'}): "
                f"{len(records)} group(s).",
            )
        )

        step_start = time.time()
        kept, counts = filter_vogs(records, hosts, self.vogdb_filters)
        self.proteins_in += sum(record.protein_count for record in records.values())
        self.excluded_uninformative += counts.dropped_uninformative
        logger.info(
            f"VOGDB {release}: {counts.kept} of {counts.total} group(s) kept "
            f"({counts.dropped_not_eukaryotic} not eukaryotic, "
            f"{counts.dropped_not_virus_only} not virus-only, "
            f"{counts.dropped_lca_rank} above {self.vogdb_filters.min_lca_rank} "
            f"rank, {counts.dropped_too_small} too small, "
            f"{counts.dropped_uninformative} uninformative, "
            f"{counts.dropped_taxon} outside the taxon)"
        )
        steps.append(
            StepInfo.from_times(
                "Filter VOGDB",
                step_start,
                time.time(),
                f"Kept {counts.kept} of {counts.total} group(s): {asdict(counts)}.",
            )
        )

        extracted: dict[str, str] = {}
        if not self.metadata_only:
            step_start = time.time()
            extracted = extract_profiles(
                paths[VOGDB_PROFILE_ARCHIVE],
                {record.group for record in kept},
                os.path.join(self.workdir, "vogdb"),
                namespace=SOURCE_NAMESPACES["vogdb"],
            )
            steps.append(
                StepInfo.from_times(
                    "Extract VOGDB profiles",
                    step_start,
                    time.time(),
                    f"Extracted {len(extracted)} profile(s) from "
                    f"{VOGDB_PROFILE_ARCHIVE}.",
                )
            )

        for record in sorted(kept, key=lambda item: item.group):
            hmm_path = extracted.get(record.group)
            if hmm_path is None and not self.metadata_only:
                logger.debug(f"{record.group} was not in the profile archive")
                continue
            profile_id = f"{SOURCE_NAMESPACES['vogdb']}__{record.group}"
            if hmm_path is not None:
                self.profiles.append(
                    BuiltProfile(
                        profile_id=profile_id,
                        hmm_path=hmm_path,
                        representative=record.group,
                        members=record.member_taxids,
                        n_seqs=record.protein_count,
                        length=profile_length(hmm_path),
                    )
                )
            member_species = [
                species[taxid] for taxid in record.member_taxids if taxid in species
            ]
            member_hosts = [
                hosts[taxid][1] for taxid in record.member_taxids if taxid in hosts
            ]
            self.rows.append(
                self._prebuilt_row(
                    record,
                    profile_id,
                    "vogdb",
                    length=profile_length(hmm_path) if hmm_path else 0,
                    species=plurality(member_species, fallback="Unclassified"),
                    host=plurality(member_hosts, minimum=0.5, fallback="Unknown"),
                )
            )

    # -- rvdb --------------------------------------------------------------- #
    def _add_rvdb(self, steps: "list[StepInfo]") -> None:
        """Download an RVDB-prot release, filter it, and add its profiles.

        RVDB records each family's LCA as a bare taxid, so the taxon names the
        mandatory taxonomy field needs are resolved through the ``datasets`` CLI
        and cached; the protein name is composed from RVDB's scored keyword bag.
        """
        step_start = time.time()
        workdir = check_outdir(f"{self.outdir}/rvdb")
        release = self.release if self.release != "latest" else RVDB_DEFAULT_RELEASE
        paths = download_rvdb(
            workdir, release=release, metadata_only=self.metadata_only
        )
        records = load_rvdb(
            paths["rvdb.sqlite.xz"], cache_path=os.path.join(workdir, "taxonomy.json")
        )
        steps.append(
            StepInfo.from_times(
                "Download RVDB-prot",
                step_start,
                time.time(),
                f"Downloaded RVDB-prot v{release} "
                f"({'metadata only' if self.metadata_only else 'with profiles'}): "
                f"{len(records)} family(ies).",
            )
        )
        self._add_prebuilt(
            "rvdb",
            records,
            steps,
            archive=paths.get("rvdb.hmm.xz"),
            release=release,
            id_from_group=lambda group: f"{SOURCE_NAMESPACES['rvdb']}__{group}",
        )

    # -- neordrp ------------------------------------------------------------ #
    def _add_neordrp(self, steps: "list[StepInfo]") -> None:
        """Download NeoRdRp, infer what taxonomy it allows, and add its profiles.

        NeoRdRp's annotation carries no taxonomy field; the only signal is in the
        names of the CDD/RdRp-scan signatures its seeds matched, so most profiles
        do not reach genus and the default ``--min-lca-rank genus`` drops them.
        It is an RdRp-sensitivity booster to merge with another source, not a
        standalone database.
        """
        step_start = time.time()
        workdir = check_outdir(f"{self.outdir}/neordrp")
        release = self.release if self.release != "latest" else NEORDRP_DEFAULT_RELEASE
        paths = download_neordrp(
            workdir, release=release, metadata_only=self.metadata_only
        )
        records = load_neordrp(paths["neordrp.annotation.tsv.xz"])
        steps.append(
            StepInfo.from_times(
                "Download NeoRdRp",
                step_start,
                time.time(),
                f"Downloaded NeoRdRp {release} "
                f"({'metadata only' if self.metadata_only else 'with profiles'}): "
                f"{len(records)} profile(s) with an annotation.",
            )
        )
        counter = {"n": 0}

        def _identifier(group: str) -> str:
            counter["n"] += 1
            return (
                f"{SOURCE_NAMESPACES['neordrp']}__{slugify(group)}"
                f"__{counter['n']:05d}"
            )

        self._add_prebuilt(
            "neordrp",
            records,
            steps,
            archive=paths.get("neordrp.hmm.xz"),
            release=release,
            id_from_group=_identifier,
        )

    # -- shared prebuilt path ------------------------------------------------ #
    def _add_prebuilt(
        self,
        source: str,
        records: "dict[str, ProfileRecord]",
        steps: "list[StepInfo]",
        archive: "str | None",
        release: str,
        id_from_group,
    ) -> None:
        """Filter a prebuilt source, extract its profiles and add the rows.

        Shared by every source that ships one concatenated ``.hmm`` file. The
        filters run on the metadata alone, which is what makes
        ``--metadata-only`` able to size a release without the profile archive.
        """
        step_start = time.time()
        kept, counts = filter_vogs(records, {}, self.vogdb_filters)
        self.proteins_in += sum(record.protein_count for record in records.values())
        self.excluded_uninformative += counts.dropped_uninformative
        logger.info(
            f"{source} {release}: {counts.kept} of {counts.total} profile(s) kept "
            f"({counts.dropped_lca_rank} above "
            f"{self.vogdb_filters.min_lca_rank} rank, "
            f"{counts.dropped_too_small} too small, "
            f"{counts.dropped_uninformative} uninformative, "
            f"{counts.dropped_taxon} outside the taxon)"
        )
        steps.append(
            StepInfo.from_times(
                f"Filter {source}",
                step_start,
                time.time(),
                f"Kept {counts.kept} of {counts.total} profile(s): {asdict(counts)}.",
            )
        )

        identifiers = {record.group: id_from_group(record.group) for record in kept}
        extracted: dict[str, str] = {}
        if not self.metadata_only and archive:
            step_start = time.time()
            extracted = extract_profiles_from_hmm(
                archive, identifiers, os.path.join(self.workdir, source)
            )
            steps.append(
                StepInfo.from_times(
                    f"Extract {source} profiles",
                    step_start,
                    time.time(),
                    f"Extracted {len(extracted)} profile(s) from {archive}.",
                )
            )

        for record in sorted(kept, key=lambda item: item.group):
            profile_id = identifiers[record.group]
            hmm_path = extracted.get(record.group)
            if hmm_path is None and not self.metadata_only:
                logger.debug(f"{record.group} was not in {archive}")
                continue
            if hmm_path is not None:
                self.profiles.append(
                    BuiltProfile(
                        profile_id=profile_id,
                        hmm_path=hmm_path,
                        representative=record.group,
                        members=[record.group],
                        n_seqs=record.protein_count,
                        length=profile_length(hmm_path),
                    )
                )
            self.rows.append(
                self._prebuilt_row(
                    record,
                    profile_id,
                    source,
                    length=profile_length(hmm_path) if hmm_path else 0,
                )
            )

    def _prebuilt_row(
        self,
        record: "ProfileRecord",
        profile_id: str,
        source: str,
        length: int = 0,
        species: str = "",
        host: str = "",
        description: "str | None" = None,
    ) -> dict:
        """Map one prebuilt profile onto the metadata columns ``screening`` consumes.

        Shared by every prebuilt source: the taxonomy comes from the record's
        lineage (whatever the source could supply), ``Molecule_type`` from the
        family via the bundled ICTV table with the source's own genome type as a
        fallback, and ``Host`` is optional (``Unknown``).
        """
        genus, family = record.genus_family
        family = family or "Unclassified"
        genus = genus or "Unclassified"
        mol_type = molecule_type_for_family(family) or record.mol_type
        return {
            "Accession": profile_id,
            "Species": species or record.species or "Unclassified",
            "Genus": genus,
            "Family": family,
            "Molecule_type": mol_type,
            "Protein": standardize_protein(
                clean_vogdb_description(
                    record.description if description is None else description
                ),
                mol_type,
                "virus",
            ),
            "Host": host or "Unknown",
            "Protein_votes": f"category:{record.category}" if record.category else "",
            "Protein_agreement": "1.00",
            "Profile_seqs": record.protein_count,
            "Profile_length": length,
            "Source": source,
            "LCA_rank": record.lca_rank,
        }


class GetHostHmmDatabase(HmmDatabaseBuilder):
    """Build the host profile database compared against by the EE filter.

    Runs on instantiation. Writes ``{outdir}/{prefix}.hmm`` (+ index),
    ``{prefix}.csv``, ``{prefix}.clusters.tsv``, ``{prefix}.protein_votes.tsv``
    and ``{prefix}.log``.

    Parameters
    ----------
    outdir : str
        Output directory.
    prefix : str
        Basename of the outputs.
    taxa : list[str]
        Host taxa to download (repeatable): several related proteomes make host
        families as deep as the viral ones, which is the most effective fix for
        the profile-depth asymmetry of the comparison.
    proteome : str, optional
        A local protein FASTA to build from instead of downloading.
    exclude_uninformative : bool
        Drop hypothetical/uncharacterised products (default, see the module
        docstring).
    keep_conserved_hypothetical : bool
        Keep uninformative products that end up in clusters of at least
        ``settings.min_profile_seqs`` members -- the conserved ones, which are
        the risky false-positive source.
    exclude_transposon_like : bool
        Drop transposon/retroelement-like products (default).
    settings : BuilderSettings
        Clustering/curation/build settings; must match the viral database's.
    dedup : bool
        Collapse 100%-identical duplicates before clustering.
    threads : int
        Threads for the external tools.
    """

    source = "host"

    def __init__(
        self,
        outdir: str,
        prefix: str = "host_hmm",
        taxa: "list[str] | None" = None,
        proteome: "str | None" = None,
        exclude_uninformative: bool = True,
        keep_conserved_hypothetical: bool = False,
        exclude_transposon_like: bool = True,
        settings: "BuilderSettings | None" = None,
        dedup: bool = True,
        threads: int = 1,
    ) -> None:
        settings = settings or BuilderSettings()
        settings.hmmer_version = settings.hmmer_version or hmmer_version()
        super().__init__(outdir, prefix, settings, threads)
        self.taxa = list(taxa or [])
        self.proteome = proteome
        self.exclude_uninformative = exclude_uninformative
        self.keep_conserved_hypothetical = keep_conserved_hypothetical
        self.exclude_transposon_like = exclude_transposon_like
        self.dedup = dedup
        if not self.taxa and not self.proteome:
            raise ValueError("host-hmm needs at least one --taxon or a --proteome")

        self.build()

    def _acquire(self) -> "list[str]":
        """Return the host protein FASTA(s), downloading one per taxon if needed."""
        if self.proteome:
            return [self.proteome]
        paths = []
        for index, taxon in enumerate(self.taxa):
            download_prefix = f"{self.prefix}_proteome{index}"
            logger.info(f"Downloading host proteins for taxon '{taxon}'")
            GetDatabases(
                dataset="host",
                taxon=taxon,
                outdir=self.outdir,
                prefix=download_prefix,
                exclude_uninformative=self.exclude_uninformative,
                standardize_proteins=False,
                cluster=False,
                threads=self.threads,
            )
            paths.append(f"{self.outdir}/{download_prefix}.fa")
        return paths

    def build(self) -> None:
        """Filter, cluster, curate, build and write every output."""
        start_time = time.time()
        steps: list[StepInfo] = []

        step_start = time.time()
        sources = self._acquire()
        products: dict[str, str] = {}
        sequences: dict[str, str] = {}
        proteins_in = 0
        excluded_uninformative = excluded_transposon = 0
        filtered_fasta = os.path.join(self.workdir, f"{self.prefix}.filtered.faa")
        with open(filtered_fasta, "w") as out:
            for path in sources:
                for record in SeqIO.parse(path, "fasta"):
                    proteins_in += 1
                    product = parse_protein_header(record.description).product
                    if self.exclude_uninformative and is_uninformative(product):
                        # Conserved uninformative families can be recovered later
                        # from the cluster sizes; see keep_conserved_hypothetical.
                        if not self.keep_conserved_hypothetical:
                            excluded_uninformative += 1
                            continue
                    if self.exclude_transposon_like and is_transposon_like(product):
                        excluded_transposon += 1
                        continue
                    products[record.id] = product
                    sequences[record.id] = str(record.seq)
                    out.write(f">{record.id}\n{record.seq}\n")
        logger.info(
            f"Host proteins: {len(sequences)} kept, {excluded_uninformative} "
            f"uninformative and {excluded_transposon} transposon-like dropped "
            f"(of {proteins_in})"
        )
        steps.append(
            StepInfo.from_times(
                "Filter host proteome",
                step_start,
                time.time(),
                f"{proteins_in} protein(s) read from {len(sources)} proteome(s); "
                f"{excluded_uninformative} uninformative and "
                f"{excluded_transposon} transposon-like record(s) dropped.",
            )
        )

        clustered_identical = 0
        if self.dedup and sequences:
            step_start = time.time()
            clustered_identical = cluster_identical_proteins(
                filtered_fasta, self.threads
            )
            sequences = read_sequences(filtered_fasta)
            logger.info(
                f"Removed {clustered_identical} identical duplicate(s); "
                f"{len(sequences)} representative(s) kept"
            )
            steps.append(
                StepInfo.from_times(
                    "Deduplicate",
                    step_start,
                    time.time(),
                    f"cd-hit (100%/100%) removed {clustered_identical} duplicate(s).",
                )
            )

        step_start = time.time()
        clusters = self._cluster(filtered_fasta, self.prefix) if sequences else {}
        counters: dict[str, int] = {}
        for representative, members in sorted(clusters.items()):
            names = [products.get(member, "") for member in members]
            consensus = consensus_protein_name(names, mode="plurality")
            if (
                self.exclude_uninformative
                and self.keep_conserved_hypothetical
                and is_uninformative(consensus.name)
                and len(members) < self.settings.min_profile_seqs
            ):
                continue
            counters[consensus.name] = counters.get(consensus.name, 0) + 1
            profile_id = make_profile_id(
                SOURCE_NAMESPACES[self.source],
                "",
                consensus.name,
                counters[consensus.name],
            )
            profile = self._build_cluster(
                sequences,
                ClusterSpec(
                    profile_id=profile_id,
                    representative=representative,
                    members=members,
                ),
            )
            if profile is None:
                continue
            self._record_votes(profile_id, consensus)
            self.rows.append(
                {
                    "Accession": profile.profile_id,
                    "Protein": consensus.name,
                    "Protein_votes": format_votes(consensus.votes),
                    "Protein_agreement": f"{consensus.agreement:.2f}",
                    "Profile_seqs": profile.n_seqs,
                    "Profile_length": profile.length,
                    "Source": ",".join(self.taxa) or str(self.proteome),
                }
            )
        steps.append(
            StepInfo.from_times(
                "Build host profiles",
                step_start,
                time.time(),
                f"Built {len(self.profiles)} profile(s) from "
                f"{self.clusters_formed} cluster(s); {self.clusters_skipped} "
                "skipped as too small/short.",
            )
        )

        step_start = time.time()
        hmm_path = self._write_outputs(HOST_METADATA_COLUMNS)
        steps.append(
            StepInfo.from_times(
                "Write database",
                step_start,
                time.time(),
                f"Wrote and pressed {len(self.profiles)} profile(s) to {hmm_path}.",
            )
        )

        counts = ProfileCounts(
            proteins_in=proteins_in,
            excluded_uninformative=excluded_uninformative,
            excluded_transposon_like=excluded_transposon,
            clustered_identical=clustered_identical,
            clusters_formed=self.clusters_formed,
            clusters_skipped=self.clusters_skipped,
            profiles_built=len(self.profiles),
        )
        self._write_log(
            {
                "source": self.source,
                "taxa": self.taxa,
                "proteome": self.proteome,
                "outdir": self.outdir,
                "prefix": self.prefix,
                "exclude_uninformative": self.exclude_uninformative,
                "keep_conserved_hypothetical": self.keep_conserved_hypothetical,
                "exclude_transposon_like": self.exclude_transposon_like,
                "dedup": self.dedup,
                "threads": self.threads,
            },
            counts,
            start_time,
            steps,
        )


def read_builder_settings(hmm_path: str) -> "BuilderSettings | None":
    """Load the :class:`BuilderSettings` recorded next to a profile database.

    The settings live in the ``{prefix}.log`` written beside ``{prefix}.hmm``.
    Returns ``None`` when there is no log (e.g. a database built by hand), which
    callers treat as "cannot verify" rather than as a mismatch.
    """
    log_path = f"{os.path.splitext(hmm_path)[0]}.log"
    if not Path(log_path).is_file():
        return None
    try:
        with open(log_path) as handle:
            payload = json.load(handle)
        recorded = payload.get("builder_settings")
    except (json.JSONDecodeError, OSError):
        return None
    if not isinstance(recorded, dict):
        return None
    known = {
        key: value
        for key, value in recorded.items()
        if key in BuilderSettings().to_dict()
    }
    return BuilderSettings(**known)
