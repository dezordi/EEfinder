"""Unit tests for eefinder.prepare_data."""

from __future__ import annotations

from Bio import SeqIO

from eefinder.prepare_data import InsertPrefix


def test_insert_prefix_rewrites_all_headers(fasta_factory, tmp_path):
    fasta = fasta_factory("genome.fasta", {"ctg1": "ACGTACGT", "ctg2": "TTTTGGGG"})

    InsertPrefix(str(fasta), "PFX", str(tmp_path))

    out = tmp_path / "PFX.rn"
    headers = [rec.id for rec in SeqIO.parse(str(out), "fasta")]
    assert headers == ["PFX/ctg1", "PFX/ctg2"]


def test_insert_prefix_preserves_sequences(fasta_factory, tmp_path):
    fasta = fasta_factory("genome.fasta", {"ctg1": "ACGTACGT"})

    InsertPrefix(str(fasta), "PFX", str(tmp_path))

    out = tmp_path / "PFX.rn"
    record = next(SeqIO.parse(str(out), "fasta"))
    assert str(record.seq) == "ACGTACGT"
