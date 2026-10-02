"""Unit tests for eefinder.utils."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pytest

from eefinder.utils import (
    EXPECTED_METADATA_COLUMNS,
    RunArguments,
    RunInfo,
    StepInfo,
    check_metadata_columns,
    check_metadata_file,
    check_outdir,
)
from eefinder.versions import DependencyVersion, SystemInfo


def test_check_outdir_creates_dir_and_strips_slash(tmp_path):
    target = f"{tmp_path}/results/"

    returned = check_outdir(target)

    assert returned == f"{tmp_path}/results"  # trailing slash removed
    assert Path(returned).is_dir()


def test_check_outdir_is_idempotent(tmp_path):
    target = f"{tmp_path}/results"
    check_outdir(target)
    # Calling again on an existing directory must not raise.
    assert check_outdir(target) == target


def _write_metadata(tmp_path: Path, header: list) -> str:
    """Write a one-row metadata CSV with the given header and return its path."""
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(",".join(header) + "\n" + ",".join("x" for _ in header) + "\n")
    return str(metadata)


def test_check_metadata_columns_accepts_the_expected_header():
    assert check_metadata_columns(list(EXPECTED_METADATA_COLUMNS)) == (
        EXPECTED_METADATA_COLUMNS
    )


def test_check_metadata_columns_reorders_and_drops_extras(caplog):
    shuffled = ["Host", "Extra", "Accession", "Protein", "Molecule_type"]
    shuffled += ["Family", "Genus", "Species"]

    with caplog.at_level("WARNING", logger="eefinder"):
        returned = check_metadata_columns(shuffled)

    # The caller gets the canonical order back, whatever the file order was.
    assert returned == EXPECTED_METADATA_COLUMNS
    assert "extra column(s): Extra" in caplog.text
    assert "not in the expected order" in caplog.text


def test_check_metadata_columns_raises_on_missing_column():
    incomplete = [c for c in EXPECTED_METADATA_COLUMNS if c != "Molecule_type"]

    with pytest.raises(ValueError, match="Molecule_type"):
        check_metadata_columns(incomplete)


def test_check_metadata_file_reads_only_the_header(tmp_path):
    path = _write_metadata(tmp_path, list(EXPECTED_METADATA_COLUMNS))

    assert check_metadata_file(path) == EXPECTED_METADATA_COLUMNS


def test_check_metadata_file_raises_on_bad_header(tmp_path):
    path = _write_metadata(tmp_path, ["Accession", "Species"])

    with pytest.raises(ValueError, match="Genus"):
        check_metadata_file(path)


def test_step_info_computes_minutes():
    info = StepInfo.from_times(
        step="demo",
        start_time=1_000_000_000,
        end_time=1_000_000_060,
        message="done",
    )
    assert info.step == "demo"
    assert info.message == "done"
    assert info.total_time_minutes == "1.0000"


def _run_arguments() -> RunArguments:
    return RunArguments(
        genome_file="genome.fa",
        prefix="PFX",
        outdir="out",
        database="db.fa",
        dbmetadata="meta.csv",
        baits="baits.fa",
        mode="blastx",
        length=10000,
        flank=10000,
        limit=1,
        range_junction=100,
        mask_per=50,
        clean_masked=False,
        threads=1,
        removetmp=False,
        index_databases=True,
        merge_level="genus",
        analysis="virus",
        overlap="keep",
        target_families=[],
        non_target_families=[],
        translation_method="default",
    )


def test_run_info_maps_arguments_and_is_json_serialisable():
    arguments = _run_arguments()
    dependency = DependencyVersion(
        name="bedtools", detected="2.27.1", expected="2.27.1", status="ok"
    )
    system = SystemInfo(
        operating_system="Linux-x", machine="x86_64", hostname="host", user="tester"
    )
    info = RunInfo.from_run(
        eefinder_version="1.1.1",
        system=system,
        arguments=arguments,
        dependencies=[dependency],
        start_time=1_000_000_000,
        end_time=1_000_000_120,
        steps_information=[],
    )

    assert info.eefinder_version == "1.1.1"
    assert info.total_time_minutes == "2.0000"
    assert info.arguments.merge_level == "genus"  # regression: was mis-mapped
    assert info.steps_information == []

    # asdict() yields the nested dict structure written to eefinder.log.
    dumped = asdict(info)
    assert dumped["eefinder_version"] == "1.1.1"
    assert dumped["system"]["operating_system"] == "Linux-x"
    assert dumped["arguments"]["merge_level"] == "genus"
    assert dumped["arguments"]["overlap"] == "keep"
    assert dumped["arguments"]["target_families"] == []
    assert dumped["arguments"]["non_target_families"] == []
    assert dumped["dependencies"] == [
        {
            "name": "bedtools",
            "detected": "2.27.1",
            "expected": "2.27.1",
            "status": "ok",
        }
    ]
