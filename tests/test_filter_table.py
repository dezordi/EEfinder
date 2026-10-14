"""Unit tests for eefinder.filter_table.FilterTable."""

from __future__ import annotations

import pandas as pd

from eefinder.filter_table import FilterTable


def test_filter_table_ee_outputs(blast_outfmt6, tmp_path):
    FilterTable(str(blast_outfmt6), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")

    assert set(filtred["qseqid"]) == {"ctg1", "ctg2"}
    assert "PROT_B" not in set(filtred["sseqid"])
    assert "PROT_D" not in set(filtred["sseqid"])


def test_filter_table_negative_strand_is_swapped(blast_outfmt6, tmp_path):
    FilterTable(str(blast_outfmt6), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")
    ctg2 = filtred[filtred["qseqid"] == "ctg2"].iloc[0]

    assert ctg2["sense"] == "neg"
    assert ctg2["qstart"] == 260
    assert ctg2["qend"] == 500


def test_filter_table_bed_name_and_bed_file(blast_outfmt6, tmp_path):
    FilterTable(str(blast_outfmt6), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")
    ctg1 = filtred[filtred["qseqid"] == "ctg1"].iloc[0]
    assert ctg1["bed_name"] == "ctg1:9-310"  # BED start: qstart - 1
    assert (filtred["tag"] == "EE").all()

    bed = pd.read_csv(
        f"{blast_outfmt6}.filtred.bed",
        sep="\t",
        header=None,
        names=["qseqid", "qstart", "qend"],
    )
    assert len(bed) == 2
    assert list(bed.columns) == ["qseqid", "qstart", "qend"]
    assert bed[bed["qseqid"] == "ctg1"].iloc[0].tolist() == ["ctg1", 9, 310]


def test_filter_table_host_tag_uses_qseqid_as_bed_name(blast_outfmt6, tmp_path):
    FilterTable(
        str(blast_outfmt6), rangejunction=100, tag="HOST", out_dir=str(tmp_path)
    )

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")
    assert (filtred["tag"] == "HOST").all()
    assert (filtred["bed_name"] == filtred["qseqid"]).all()


def test_filter_table_without_negative_hits_keeps_integer_coordinates(tmp_path):
    # with no negative-sense hit to swap, the coordinates used to become floats
    # (6.0), which bedtools rejects downstream in `merge`
    raw = tmp_path / "pos_only.tsv"
    raw.write_text(
        "ctg1\tPROT_A\t51.6\t91\t44\t0\t6\t278\t1739\t1829\t1.39e-23\t92.5\n"
        "ctg2\tPROT_B\t36.9\t92\t58\t0\t3\t278\t253\t344\t6.39e-10\t53.4\n"
    )
    FilterTable(str(raw), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{raw}.filtred", sep="\t")
    assert filtred["qstart"].dtype == "int64"
    assert set(filtred["bed_name"]) == {"ctg1:5-278", "ctg2:2-278"}
    bed = (tmp_path / "pos_only.tsv.filtred.bed").read_text()
    assert ".0" not in bed
