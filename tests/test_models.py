"""Unit tests for eefinder.models."""

from __future__ import annotations
import pytest
from eefinder.models import (
    CANDIDATE_BED,
    ELEMENT_TABLE,
    HIT_TABLE,
    FastaModel,
    ModelError,
    TableModel,
    validate_all,
)


def test_table_model_accepts_a_conforming_header(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text("Element-ID\tSense\tother\nx\tpos\t1\n")

    model = TableModel(name="demo", columns=("Element-ID", "Sense"))

    assert model.validate(str(path)) is None  # no exception


def test_table_model_names_the_missing_columns(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text("Element-ID\n x\n")

    model = TableModel(name="demo", columns=("Element-ID", "Sense", "Host"))

    with pytest.raises(ModelError, match="Sense, Host"):
        model.validate(str(path))


def test_table_model_rejects_extra_columns_when_asked(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text("a\tb\tc\n1\t2\t3\n")

    model = TableModel(name="demo", columns=("a", "b"), extra_columns=False)

    with pytest.raises(ModelError, match="expected exactly"):
        model.validate(str(path))


def test_headerless_model_counts_fields(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("q\ts\t30.0\n")  # only 3 of the 12 outfmt6 fields

    with pytest.raises(ModelError, match="expected at least 12"):
        HIT_TABLE.validate(str(path))


def test_headerless_model_accepts_outfmt6(tmp_path, blast_outfmt6):
    # The shared fixture writes a real outfmt6 table.
    assert HIT_TABLE.validate(str(blast_outfmt6)) is None


def test_empty_file_is_allowed_by_default(tmp_path):
    path = tmp_path / "empty.tsv"
    path.touch()

    assert HIT_TABLE.validate(str(path)) is None


def test_empty_file_can_be_rejected(tmp_path):
    path = tmp_path / "empty.tsv"
    path.touch()

    model = TableModel(name="demo", columns=("a",), allow_empty=False)

    with pytest.raises(ModelError, match="is empty"):
        model.validate(str(path))


def test_missing_file_is_a_violation(tmp_path):
    with pytest.raises(ModelError, match="was not produced"):
        CANDIDATE_BED.validate(str(tmp_path / "absent.bed"))


def test_fasta_model_accepts_a_fasta(tmp_path):
    path = tmp_path / "seqs.fa"
    path.write_text(">a\nACGT\n")

    assert FastaModel(name="demo").validate(str(path)) is None


def test_fasta_model_rejects_a_non_fasta(tmp_path):
    path = tmp_path / "seqs.fa"
    path.write_text("ACGT\n")

    with pytest.raises(ModelError, match="does not start with a FASTA header"):
        FastaModel(name="demo").validate(str(path))


def test_fasta_model_can_require_records(tmp_path):
    path = tmp_path / "seqs.fa"
    path.touch()

    with pytest.raises(ModelError, match="has no records"):
        FastaModel(name="demo", allow_empty=False).validate(str(path))


def test_validate_all_skips_none_paths(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text("\t".join(ELEMENT_TABLE.columns) + "\n")

    # None marks an output the run's options switched off.
    assert validate_all([(ELEMENT_TABLE, str(path)), (ELEMENT_TABLE, None)]) is None


def test_validate_all_reports_the_first_violation(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("nope\n")

    with pytest.raises(ModelError):
        validate_all([(ELEMENT_TABLE, str(bad))])
