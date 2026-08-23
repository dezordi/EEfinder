"""End-to-end integration scenarios for the profile-HMM engine (``-md hmmer``).

These build real profile databases from ``test_files/`` with cd-hit, mafft,
CIAlign and hmmbuild, then run the pipeline against them. They are skipped when
the HMM toolchain is absent, so the rest of the suite still runs on a bare
install.

The scenarios cover what a maintainer would not want to regress:

* a viral profile database is built with exact taxonomy and a voted protein name;
* the host database excludes hypothetical and transposon-like products;
* ``-md hmmer`` runs end to end and produces the documented outputs;
* the symmetric host filter refuses databases built with different settings;
* ``-md hmmer`` without ``-bth`` is an error, not a silent downgrade.
"""

from __future__ import annotations

import json
import shutil
import subprocess

import pandas as pd
import pytest

from conftest import binaries_available

REQUIRED_BINARIES = (
    "eefinder",
    "bedtools",
    "cd-hit",
    "mafft",
    "CIAlign",
    "hmmbuild",
    "hmmpress",
    "hmmsearch",
)

PREFIX = "Ae_aeg_Aag2_ctg_1913"

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not binaries_available(*REQUIRED_BINARIES),
        reason=f"requires {', '.join(REQUIRED_BINARIES)} on PATH",
    ),
]

#: Builder options shared by the viral and host databases. They MUST match for
#: the symmetric host filter to run (see proposal.md §5.4).
SHARED_BUILDER_OPTIONS = (
    "--min-profile-seqs",
    "1",
    "--keep-singletons",
    "-p",
    "2",
)


def _run(*args):
    """Run the CLI, failing the test with its output when it exits non-zero."""
    result = subprocess.run(
        ["eefinder", *[str(arg) for arg in args]],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def _build_viral_db(outdir, virus_db, virus_metadata, *extra):
    _run(
        "get-databases",
        "virus-hmm",
        "--proteins",
        virus_db,
        "--metadata",
        virus_metadata,
        "-od",
        outdir,
        "-pr",
        "virus_hmm",
        *SHARED_BUILDER_OPTIONS,
        *extra,
    )
    return outdir / "virus_hmm.hmm", outdir / "virus_hmm.csv"


def _build_host_db(outdir, proteome, *extra):
    _run(
        "get-databases",
        "host-hmm",
        "--proteome",
        proteome,
        "-od",
        outdir,
        "-pr",
        "host_hmm",
        *SHARED_BUILDER_OPTIONS,
        *extra,
    )
    return outdir / "host_hmm.hmm"


@pytest.fixture(scope="module")
def hmm_databases(tmp_path_factory, virus_db, virus_metadata, filter_db):
    """A viral and a host profile database built with identical settings."""
    outdir = tmp_path_factory.mktemp("hmm_db")
    viral_hmm, viral_csv = _build_viral_db(outdir, virus_db, virus_metadata)
    host_hmm = _build_host_db(outdir, filter_db)
    return viral_hmm, viral_csv, host_hmm


def test_viral_database_is_built_with_exact_taxonomy(hmm_databases):
    """Clustering inside one family makes Family exact, not inferred."""
    viral_hmm, viral_csv, _ = hmm_databases
    assert viral_hmm.is_file()
    for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
        assert viral_hmm.with_suffix(f".hmm{suffix}").is_file()

    metadata = pd.read_csv(viral_csv)
    assert not metadata.empty
    # Every profile id is namespaced, and the family in the id matches the column.
    for row in metadata.itertuples():
        assert row.Accession.startswith("NCBIREFSEQ__")
        assert row.Accession.split("__")[1] == row.Family.replace(" ", "_")
        assert "|" not in row.Accession
    assert set(metadata.LCA_rank) == {"family"}
    assert (metadata.Profile_seqs >= 1).all()
    # No profile id may be a substring of another (GetFinalTaxonomy matches with
    # a substring test).
    ids = list(metadata.Accession)
    for candidate in ids:
        assert not any(candidate in other for other in ids if other != candidate)


def test_viral_database_records_its_builder_settings(hmm_databases):
    viral_hmm, _, _ = hmm_databases
    log = json.loads(viral_hmm.with_suffix(".log").read_text())
    settings = log["builder_settings"]
    assert settings["cluster_identity"] == 0.5
    assert settings["msa_curation"] is True
    assert settings["curation_remove_divergent"] is False
    assert log["profile_counts"]["profiles_built"] > 0


def test_protein_votes_report_lists_close_calls(hmm_databases):
    viral_hmm, _, _ = hmm_databases
    votes = pd.read_csv(viral_hmm.parent / "virus_hmm.protein_votes.tsv", sep="\t")
    assert list(votes.columns) == [
        "Profile",
        "Taxon",
        "Protein",
        "Agreement",
        "Votes",
    ]


def test_host_database_excludes_uninformative_and_transposon_products(
    tmp_path, filter_db
):
    """Hypothetical and retroelement-like host products must not become profiles."""
    proteome = tmp_path / "host_proteome.faa"
    proteome.write_text(
        filter_db.read_text()
        + ">HYPO_1 hypothetical protein LOC123 [Aedes aegypti]\n"
        + "MKAILVGTSGAGKSTLLQALNRLYELDSGSIRIDGVDIRDLDPVELRRHIGYVPQDPFLFS\n"
        + ">RT_1 endonuclease-reverse transcriptase [Aedes aegypti]\n"
        + "MTAVTVAQAFVSSWIARFGVPVKLTTDLGRQFESELFRELTRILGITHLKTTPYHPQANG\n"
    )
    outdir = tmp_path / "host_db"
    host_hmm = _build_host_db(outdir, proteome)

    log = json.loads(host_hmm.with_suffix(".log").read_text())
    counts = log["profile_counts"]
    assert counts["excluded_uninformative"] == 1
    assert counts["excluded_transposon_like"] == 1
    assert counts["profiles_built"] == 5  # the five original baits

    metadata = pd.read_csv(outdir / "host_hmm.csv")
    products = " ".join(metadata.Protein.astype(str)).lower()
    assert "hypothetical" not in products
    assert "reverse transcriptase" not in products


def test_hmm_screening_runs_end_to_end(tmp_path, genome_file, hmm_databases):
    viral_hmm, viral_csv, host_hmm = hmm_databases
    outdir = tmp_path / "run"
    _run(
        "screening",
        "-in",
        genome_file,
        "-od",
        outdir,
        "-md",
        "hmmer",
        "-dbh",
        viral_hmm,
        "-mth",
        viral_csv,
        "-bth",
        host_hmm,
        "-ln",
        "1000",
        "-p",
        "2",
    )

    for output in (
        f"{PREFIX}.EEs.fa",
        f"{PREFIX}.EEs.tax.tsv",
        f"{PREFIX}.EEs.gff3",
        f"{PREFIX}.EEs.flanks.fa",
    ):
        assert (outdir / output).is_file(), output

    taxonomy = pd.read_csv(outdir / f"{PREFIX}.EEs.tax.tsv", sep="\t")
    assert not taxonomy.empty
    # Taxonomy resolved through the profile metadata, not left unjoined.
    assert (taxonomy["Family"] != "Unclassified").any()
    assert taxonomy["Protein-IDs"].str.contains("NCBIREFSEQ__").all()

    run_log = json.loads((outdir / "eefinder.log").read_text())
    assert run_log["arguments"]["mode"] == "hmmer"
    assert run_log["arguments"]["hmm_host_filter"] == "hmm"
    assert run_log["arguments"]["translation_method"] == "sixframe"


def test_hmm_screening_requires_a_host_database(tmp_path, genome_file, hmm_databases):
    """No -bth is an error, not a silent downgrade to a weaker filter."""
    viral_hmm, viral_csv, _ = hmm_databases
    result = subprocess.run(
        [
            "eefinder",
            "screening",
            "-in",
            str(genome_file),
            "-od",
            str(tmp_path / "out"),
            "-md",
            "hmmer",
            "-dbh",
            str(viral_hmm),
            "-mth",
            str(viral_csv),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "host-hmm" in result.stdout + result.stderr


def test_mismatched_builder_settings_stop_the_symmetric_filter(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db, hmm_databases
):
    """A host database built differently is not comparable, so the run stops."""
    viral_hmm, viral_csv, _ = hmm_databases
    other = tmp_path / "other_host"
    _run(
        "get-databases",
        "host-hmm",
        "--proteome",
        filter_db,
        "-od",
        other,
        "-pr",
        "host_hmm",
        "--min-profile-seqs",
        "1",
        "--keep-singletons",
        "--cluster-identity",
        "0.7",
        "-p",
        "2",
    )

    command = [
        "eefinder",
        "screening",
        "-in",
        str(genome_file),
        "-od",
        str(tmp_path / "out"),
        "-md",
        "hmmer",
        "-dbh",
        str(viral_hmm),
        "-mth",
        str(viral_csv),
        "-bth",
        str(other / "host_hmm.hmm"),
        "-ln",
        "1000",
        "-p",
        "2",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 1
    assert "cluster_identity" in result.stdout + result.stderr

    # ... unless the user explicitly accepts the mismatch.
    accepted = subprocess.run(
        command + ["--allow_builder_mismatch"], capture_output=True, text=True
    )
    assert accepted.returncode == 0, accepted.stdout + accepted.stderr


def test_host_filter_modes_remove_candidates_at_the_comparison_step(
    tmp_path, genome_file, filter_db, hmm_databases
):
    """Each host filter can only drop candidates, never invent them.

    The assertion is made on the table CompareResults writes, not on the final
    elements: merging joins neighbouring hits of the same family, so removing one
    hit changes the *boundaries* of the surviving elements and their ids are not
    comparable across filter modes.
    """
    viral_hmm, viral_csv, host_hmm = hmm_databases

    def _candidates(outdir):
        table = next(outdir.glob("tmp_files/*.concat.nr"))
        return set(pd.read_csv(table, sep="\t")["qseqid"])

    common = [
        "screening",
        "-in",
        genome_file,
        "-md",
        "hmmer",
        "-dbh",
        viral_hmm,
        "-mth",
        viral_csv,
        "-ln",
        "1000",
        "-p",
        "2",
    ]

    unfiltered = tmp_path / "none"
    _run(*common, "-od", unfiltered, "--hmm_host_filter", "none")

    symmetric = tmp_path / "hmm"
    _run(*common, "-od", symmetric, "--hmm_host_filter", "hmm", "-bth", host_hmm)

    fallback = tmp_path / "blastp"
    _run(*common, "-od", fallback, "--hmm_host_filter", "blastp", "-bt", filter_db)

    assert _candidates(symmetric) <= _candidates(unfiltered)
    assert _candidates(fallback) <= _candidates(unfiltered)


def test_calibration_writes_thresholds_and_tightens_the_search(
    tmp_path, genome_file, filter_db, hmm_databases
):
    """`calibrate-hmm` + `--hmm_use_ga` is stricter than the plain E-value cutoff.

    The point of a calibrated threshold is that a profile with cellular homologs
    has to beat the score it achieves on the host's own genes, so the calibrated
    run can only ever keep fewer candidates than the uncalibrated one.
    """
    viral_hmm, viral_csv, host_hmm = hmm_databases

    calibrated_dir = tmp_path / "calibrated"
    calibrated_dir.mkdir()
    calibrated = calibrated_dir / "virus_hmm.hmm"
    shutil.copyfile(viral_hmm, calibrated)
    shutil.copyfile(viral_hmm.with_suffix(".log"), calibrated.with_suffix(".log"))

    _run("calibrate-hmm", "-dbh", calibrated, "-bt", filter_db, "-p", "2")

    text = calibrated.read_text()
    assert text.count("\nGA ") == text.count("\nNAME ")
    thresholds = [
        float(line.split()[1]) for line in text.splitlines() if line.startswith("GA ")
    ]
    # Every profile carries a threshold at least as strict as the E-value floor.
    assert all(threshold > 0 for threshold in thresholds)

    def _candidate_count(outdir):
        table = next(outdir.glob("tmp_files/*.concat.nr"))
        return len(pd.read_csv(table, sep="\t"))

    common = [
        "screening",
        "-in",
        genome_file,
        "-md",
        "hmmer",
        "-mth",
        viral_csv,
        "-bth",
        host_hmm,
        "-ln",
        "1000",
        "-p",
        "2",
    ]
    plain = tmp_path / "plain"
    _run(*common, "-dbh", viral_hmm, "-od", plain)
    with_ga = tmp_path / "with_ga"
    _run(*common, "-dbh", calibrated, "-od", with_ga, "--hmm_use_ga")

    # Raising a profile's threshold can only remove hits, so the calibrated run
    # cannot yield more candidates. Their coordinates are NOT comparable: which
    # hit represents a window changes when a stronger one is filtered out, so the
    # assertion is on counts rather than on region ids.
    assert _candidate_count(with_ga) <= _candidate_count(plain)
    # At least one profile picked up host signal and is now stricter than the floor.
    assert max(thresholds) > min(thresholds) * 1.5


def test_hybrid_mode_takes_taxonomy_from_the_sequence_database(
    tmp_path, genome_file, virus_db, virus_metadata, hmm_databases
):
    """Profiles discover the elements; BLAST names them.

    A profile's taxonomy is only as precise as the group it was built from, so
    hybrid mode re-assigns each element from a sequence-level search of the
    candidate regions — and reports that hit's real percent identity instead of
    the profile's posterior accuracy.
    """
    viral_hmm, viral_csv, host_hmm = hmm_databases
    outdir = tmp_path / "hybrid"
    _run(
        "screening",
        "-in",
        genome_file,
        "-od",
        outdir,
        "-md",
        "hmmer",
        "-dbh",
        viral_hmm,
        "-mth",
        viral_csv,
        "-bth",
        host_hmm,
        "--taxonomy_refine",
        "blast",
        "-db",
        virus_db,
        "-mt",
        virus_metadata,
        "-ln",
        "1000",
        "-p",
        "2",
    )

    taxonomy = pd.read_csv(outdir / f"{PREFIX}.EEs.tax.tsv", sep="\t")
    assert not taxonomy.empty
    # Elements now cite RefSeq accessions rather than profile ids.
    assert taxonomy["Protein-IDs"].str.contains("YP_|NP_").any()
    # Identities are sequence identities again, well below the posterior
    # accuracies (80-100) a pure profile run reports.
    assert taxonomy["Average_pident"].min() < 70

    run_log = json.loads((outdir / "eefinder.log").read_text())
    assert run_log["arguments"]["taxonomy_refine"] == "blast"

    # The pre-refinement assignments are preserved for audit.
    assert list(outdir.glob("tmp_files/*.concat.nr.unrefined"))


def test_hybrid_mode_requires_the_sequence_database(
    tmp_path, genome_file, hmm_databases
):
    viral_hmm, viral_csv, host_hmm = hmm_databases
    result = subprocess.run(
        [
            "eefinder",
            "screening",
            "-in",
            str(genome_file),
            "-od",
            str(tmp_path / "out"),
            "-md",
            "hmmer",
            "-dbh",
            str(viral_hmm),
            "-mth",
            str(viral_csv),
            "-bth",
            str(host_hmm),
            "--taxonomy_refine",
            "blast",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "-db" in result.stdout + result.stderr
