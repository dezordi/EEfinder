"""End-to-end integration scenarios for the ``eefinder`` CLI on ``test_files/``.

These tests shell out to the installed ``eefinder`` console script and require
the external binaries (``blastx``, ``makeblastdb``, ``bedtools``) on ``PATH``
(i.e. an activated EEfinder environment); they are skipped otherwise so the
unit-test suite still runs on a bare Python install.

The set is intentionally small and scenario-driven rather than an exhaustive
parameter sweep. Each test exercises one behaviour a maintainer would not want
to regress:

* reproducibility of the documented run (golden-file comparison);
* the ``--clean_masked`` output being a coherent subset of the full run;
* ``--limit`` (merge distance) actually controlling element merging;
* the family ``--merge_level`` branch running end-to-end;
* ``--removetmp`` housekeeping.
"""

from __future__ import annotations

import filecmp
import shutil
import subprocess

import pandas as pd
import pytest

from conftest import binaries_available

REQUIRED_BINARIES = ("eefinder", "blastx", "makeblastdb", "bedtools")

PREFIX = "Ae_aeg_Aag2_ctg_1913"

#: The four user-facing outputs compared against ``test_files/expected_results``
#: (temporary files and the timestamped run log are intentionally ignored).
MAIN_OUTPUTS = (
    f"{PREFIX}.EEs.fa",
    f"{PREFIX}.EEs.tax.tsv",
    f"{PREFIX}.EEs.flanks.fa",
    f"{PREFIX}.EEs.gff3",
)

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not binaries_available(*REQUIRED_BINARIES),
        reason=f"requires {', '.join(REQUIRED_BINARIES)} on PATH",
    ),
]

EXPECTED_TAX_COLUMNS = {
    "Element-ID",
    "Sense",
    "Protein-IDs",
    "Protein-Products",
    "Molecule_type",
    "Family",
    "Genus",
    "Species",
    "Host",
    "Overlaped_Element_ID",
    "tag",
    "Average_pident",
}


def _run_eefinder(outdir, genome, db, meta, baits, *extra, merge_limit=100):
    """Run the CLI into ``outdir`` with the README example parameters.

    ``-id`` (makeblastdb) writes index files next to the database/baits FASTAs,
    so those are staged into a scratch dir to keep the committed ``test_files/``
    clean. Extra CLI flags may be appended via ``*extra``.
    """
    staged = outdir.parent / "inputs"
    staged.mkdir(parents=True, exist_ok=True)
    db = shutil.copy(db, staged)
    baits = shutil.copy(baits, staged)

    cmd = [
        shutil.which("eefinder"),
        "screening",
        "-in",
        str(genome),
        "-od",
        str(outdir),
        "-db",
        str(db),
        "-mt",
        str(meta),
        "-bt",
        str(baits),
        "-ln",
        "1000",
        "-p",
        "2",
        "-id",
        "-lm",
        str(merge_limit),
        *extra,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"eefinder exited with {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    return result


def _element_ids(outdir, name=f"{PREFIX}.EEs.tax.tsv"):
    """Return the set of ``Element-ID`` values from a taxonomy table."""
    return set(pd.read_csv(outdir / name, sep="\t")["Element-ID"])


def test_matches_expected_results(
    tmp_path,
    genome_file,
    virus_db,
    virus_metadata,
    filter_db,
    expected_results,
    update_expected,
):
    """The documented run reproduces the golden outputs byte-for-byte.

    Run with ``pytest --update-test`` to instead refresh
    ``test_files/expected_results/`` with the freshly produced outputs (e.g.
    after an intentional dependency-version change), then commit them.
    """
    outdir = tmp_path / "out"
    _run_eefinder(outdir, genome_file, virus_db, virus_metadata, filter_db)

    if update_expected:
        expected_results.mkdir(parents=True, exist_ok=True)
        for name in MAIN_OUTPUTS:
            shutil.copy(outdir / name, expected_results / name)
        pytest.skip(f"--update-test: refreshed {expected_results}")

    for name in MAIN_OUTPUTS:
        produced = outdir / name
        assert produced.is_file(), f"missing output: {name}"
        assert filecmp.cmp(
            produced, expected_results / name, shallow=False
        ), f"{name} differs from test_files/expected_results/{name}"

    # Guard against a corrupted golden set.
    tax = pd.read_csv(outdir / f"{PREFIX}.EEs.tax.tsv", sep="\t")
    assert EXPECTED_TAX_COLUMNS.issubset(tax.columns)

    # Without --removetmp the intermediates are archived and the run log written.
    assert (outdir / "tmp_files").is_dir()
    assert (outdir / "eefinder.log").is_file()


def test_clean_masked_is_subset_of_full_run(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """``--clean_masked`` adds a coherent, populated subset of the full run."""
    outdir = tmp_path / "out"
    _run_eefinder(outdir, genome_file, virus_db, virus_metadata, filter_db, "-cm")

    cleaned_fa = outdir / f"{PREFIX}.EEs.cleaned.fa"
    cleaned_tax = outdir / f"{PREFIX}.EEs.cleaned.tax.tsv"
    assert cleaned_fa.is_file()
    assert cleaned_tax.is_file()

    # One taxonomy row per cleaned record, matching the documented schema
    # (regression guard: cleaned IDs once mismatched the table and it was empty).
    tax = pd.read_csv(cleaned_tax, sep="\t")
    assert EXPECTED_TAX_COLUMNS.issubset(tax.columns)
    n_records = sum(
        1 for line in cleaned_fa.read_text().splitlines() if line.startswith(">")
    )
    assert n_records > 0
    assert len(tax) == n_records

    # The cleaned elements must be a subset of the full run's elements.
    assert _element_ids(outdir, cleaned_tax.name) <= _element_ids(outdir)


def test_merge_limit_controls_merging(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """A larger ``--limit`` merges neighbouring elements, reducing the count."""
    strict = tmp_path / "lm1"
    loose = tmp_path / "lm100"
    _run_eefinder(
        strict, genome_file, virus_db, virus_metadata, filter_db, merge_limit=1
    )
    _run_eefinder(
        loose, genome_file, virus_db, virus_metadata, filter_db, merge_limit=100
    )

    ids_strict = _element_ids(strict)
    ids_loose = _element_ids(loose)

    # Merging only ever joins elements, so a looser limit yields fewer of them.
    assert len(ids_loose) < len(ids_strict)
    # The adjacent same-taxon pair (~20 nt apart) merges only at the looser limit.
    merged = "ctg_1913:102863-104096"
    assert merged in ids_loose
    assert merged not in ids_strict


def test_family_merge_level_runs_end_to_end(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """The ``--merge_level family`` branch produces a valid taxonomy table."""
    outdir = tmp_path / "out"
    _run_eefinder(
        outdir, genome_file, virus_db, virus_metadata, filter_db, "-ml", "family"
    )

    tax = pd.read_csv(outdir / f"{PREFIX}.EEs.tax.tsv", sep="\t")
    assert EXPECTED_TAX_COLUMNS.issubset(tax.columns)
    assert len(tax) > 0
    assert (outdir / f"{PREFIX}.EEs.gff3").is_file()


def test_overlap_longest_filters_and_preserves_removed(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """``--overlap longest`` drops overlaping elements but keeps them on disk."""
    keep = tmp_path / "keep"
    longest = tmp_path / "longest"
    _run_eefinder(keep, genome_file, virus_db, virus_metadata, filter_db)
    _run_eefinder(
        longest, genome_file, virus_db, virus_metadata, filter_db, "-ov", "longest"
    )

    full_ids = _element_ids(keep)
    kept_ids = _element_ids(longest)
    assert kept_ids <= full_ids

    # Filtered-out elements are preserved under tmp_outputs/ rather than deleted.
    removed_tax = longest / "tmp_outputs" / f"{PREFIX}.EEs.removed.tax.tsv"
    removed_fa = longest / "tmp_outputs" / f"{PREFIX}.EEs.removed.fa"
    assert removed_tax.is_file()
    assert removed_fa.is_file()

    removed_ids = set(pd.read_csv(removed_tax, sep="\t")["Element-ID"])
    # Kept and removed partition the unfiltered run, with nothing lost.
    assert kept_ids.isdisjoint(removed_ids)
    assert kept_ids | removed_ids == full_ids
    # The test data contains at least one resolvable overlap.
    assert removed_ids


def _targets_cmd(outdir, genome_file, virus_db, virus_metadata, filter_db, *extra):
    staged = outdir.parent / "inputs"
    staged.mkdir(parents=True, exist_ok=True)
    db = shutil.copy(virus_db, staged)
    baits = shutil.copy(filter_db, staged)
    return [
        shutil.which("eefinder"),
        "screening",
        "-in",
        str(genome_file),
        "-od",
        str(outdir),
        "-db",
        str(db),
        "-mt",
        str(virus_metadata),
        "-bt",
        str(baits),
        "-ln",
        "1000",
        "-p",
        "2",
        "-id",
        "-ov",
        "targets",
        *extra,
    ]


def test_overlap_targets_requires_exactly_one_family_list(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """``--overlap targets`` needs exactly one of the two family lists."""
    # Neither list -> error.
    neither = subprocess.run(
        _targets_cmd(
            tmp_path / "neither", genome_file, virus_db, virus_metadata, filter_db
        ),
        capture_output=True,
        text=True,
    )
    assert neither.returncode != 0
    assert "target_families" in (neither.stdout + neither.stderr)

    # Both lists -> error.
    both = subprocess.run(
        _targets_cmd(
            tmp_path / "both",
            genome_file,
            virus_db,
            virus_metadata,
            filter_db,
            "-tf",
            "Flaviviridae",
            "-ntf",
            "Retroviridae",
        ),
        capture_output=True,
        text=True,
    )
    assert both.returncode != 0
    assert "target_families" in (both.stdout + both.stderr)


def test_overlap_targets_non_target_families_drops_that_family(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """``--overlap targets --non_target_families`` drops the listed family."""
    keep = tmp_path / "keep"
    _run_eefinder(keep, genome_file, virus_db, virus_metadata, filter_db)

    keep_tax = pd.read_csv(keep / f"{PREFIX}.EEs.tax.tsv", sep="\t")
    overlaped = keep_tax[keep_tax["tag"] == "overlaped"]
    if overlaped.empty:
        pytest.skip("no overlaping elements in the example run")
    dropped_family = overlaped["Family"].iloc[0]
    # Element-IDs of the overlaped members of that family (a cluster always mixes
    # families, so none are shielded by the "never wipe" rule).
    dropped_ids = set(
        overlaped.loc[overlaped["Family"] == dropped_family, "Element-ID"]
    )

    ntf = tmp_path / "ntf"
    _run_eefinder(
        ntf,
        genome_file,
        virus_db,
        virus_metadata,
        filter_db,
        "-ov",
        "targets",
        "-ntf",
        dropped_family,
    )

    kept_ids = _element_ids(ntf)
    removed_tax = ntf / "tmp_outputs" / f"{PREFIX}.EEs.removed.tax.tsv"
    assert removed_tax.is_file()
    removed_ids = set(pd.read_csv(removed_tax, sep="\t")["Element-ID"])

    # Every overlaped member of the dropped family is filtered out and preserved.
    assert dropped_ids <= removed_ids
    assert kept_ids.isdisjoint(dropped_ids)
    # Kept and removed still partition the unfiltered run.
    assert kept_ids.isdisjoint(removed_ids)
    assert kept_ids | removed_ids == _element_ids(keep)


def test_removetmp_removes_intermediate_files(
    tmp_path, genome_file, virus_db, virus_metadata, filter_db
):
    """``--removetmp`` deletes intermediates instead of archiving them."""
    outdir = tmp_path / "out"
    _run_eefinder(outdir, genome_file, virus_db, virus_metadata, filter_db, "-rm")

    assert not (outdir / "tmp_files").exists()
    # The main outputs are still produced.
    assert (outdir / f"{PREFIX}.EEs.fa").is_file()
