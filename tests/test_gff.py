"""Unit tests for eefinder.gff.WriteGFF3."""

from __future__ import annotations

import pandas as pd

from eefinder.gff import WriteGFF3


def _write_tax(path, rows):
    columns = [
        "element_id",
        "sense",
        "protein_ids",
        "protein_products",
        "molecule_type",
        "taxonomy",
        "host",
        "overlaped_element_id",
        "tag",
        "average_pident",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _read_features(path):
    lines = path.read_text().splitlines()
    header, *features = lines
    return header, [line.split("\t") for line in features]


def test_write_gff3_columns_and_coordinates(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:100-200",
                "pos",
                "P1|80.0",
                "polyprotein",
                "ssRNA(+)",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA",
                "Aedes",
                "",
                "unique",
                80.0,
            ],
            [
                "ctg-x:5-40",
                "neg",
                "P2|50.0 | P3|60.0",
                "glyco",
                "ssRNA(-)",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamB;g__GenB;s__SpB",
                "Culex",
                "ctg1:100-200",
                "overlaped",
                55.0,
            ],
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    header, features = _read_features(out)
    assert header == "##gff-version 3"
    assert len(features) == 2
    by_seqid = {cols[0]: cols for cols in features}

    assert by_seqid["ctg1"][:8] == [
        "ctg1",
        "EEfinder",
        "endogenous_viral_element",
        "101",
        "200",
        "80.0",
        "+",
        ".",
    ]

    assert by_seqid["ctg-x"][3:5] == ["6", "40"]
    assert by_seqid["ctg-x"][6] == "-"


def test_write_gff3_features_sorted_by_start(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:200-300",
                "pos",
                "P|1.0",
                "p",
                "m",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__F;g__G;s__S",
                "H",
                "",
                "u",
                1.0,
            ],
            [
                "ctg1:50-100",
                "pos",
                "P|1.0",
                "p",
                "m",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__F;g__G;s__S",
                "H",
                "",
                "u",
                1.0,
            ],
            [
                "ctg1:100-150",
                "pos",
                "P|1.0",
                "p",
                "m",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__F;g__G;s__S",
                "H",
                "",
                "u",
                1.0,
            ],
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    starts = [int(cols[3]) for cols in features]
    assert starts == sorted(starts) == [51, 101, 201]


def test_write_gff3_id_carries_prefix_to_match_fasta_headers(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:100-200",
                "pos",
                "P|1.0",
                "p",
                "m",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__F;g__G;s__S",
                "H",
                "",
                "u",
                1.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), prefix="Aaeg")

    _, features = _read_features(out)
    attrs = dict(pair.split("=", 1) for pair in features[0][8].split(";"))
    assert attrs["ID"] == "Aaeg/ctg1:100-200"


def test_write_gff3_attributes_and_escaping(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0",
                "polyprotein, partial",
                "ssRNA(+)",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA",
                "Aedes",
                "",
                "unique",
                90.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    attrs = dict(pair.split("=", 1) for pair in features[0][8].split(";"))
    assert attrs["ID"] == "ctg1:0-10"
    assert attrs["Name"] == "SpA"
    assert attrs["family"] == "FamA"
    assert attrs["overlap_status"] == "unique"
    assert attrs["product"] == "polyprotein%2C partial"


def test_write_gff3_missing_score_column(tmp_path):
    tax = tmp_path / "eves.tax"
    pd.DataFrame(
        {
            "element_id": ["ctg1:10-20"],
            "sense": ["pos"],
            "taxonomy": ["r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA"],
        }
    ).to_csv(tax, sep="\t", index=False)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    assert features[0][5] == "."
    assert features[0][6] == "+"


def _single_element_tax(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0",
                "prot",
                "ssRNA",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA",
                "Aedes",
                "",
                "unique",
                90.0,
            ]
        ],
    )
    return tax


def test_write_gff3_analysis_virus_is_endogenous_viral_element(tmp_path):
    tax = _single_element_tax(tmp_path)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), analysis="virus")

    _, features = _read_features(out)
    assert features[0][2] == "endogenous_viral_element"


def test_write_gff3_analysis_bacteria_is_endogenous_bacterial_element(tmp_path):
    tax = _single_element_tax(tmp_path)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), analysis="bacteria")

    _, features = _read_features(out)
    assert features[0][2] == "endogenous_bacterial_element"


def test_write_gff3_explicit_feature_type_overrides_analysis(tmp_path):
    tax = _single_element_tax(tmp_path)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), analysis="bacteria", feature_type="match")

    _, features = _read_features(out)
    assert features[0][2] == "match"


def test_write_gff3_custom_source_and_type(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0",
                "prot",
                "ssRNA",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA",
                "Aedes",
                "",
                "unique",
                90.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), source="MyTool", feature_type="match")

    _, features = _read_features(out)
    assert features[0][1] == "MyTool"
    assert features[0][2] == "match"


def test_write_gff3_protein_ids_are_a_comma_list(tmp_path):
    """The ``;`` of the TSV would have to be escaped; GFF3 lists use commas."""
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0;P2|80.0",
                "prot",
                "ssRNA",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA",
                "Aedes",
                "",
                "unique",
                85.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    attrs = dict(pair.split("=", 1) for pair in features[0][8].split(";"))
    assert attrs["protein_ids"] == "P1|90.0,P2|80.0"
    assert "%3B" not in features[0][8]


def test_write_gff3_escapes_a_comma_inside_a_list_entry(tmp_path):
    """A literal comma must not be read as a list separator."""
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0;P2, odd|80.0",
                "prot",
                "ssRNA",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__SpA",
                "Aedes",
                "",
                "unique",
                85.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    attrs = dict(pair.split("=", 1) for pair in features[0][8].split(";"))
    assert attrs["protein_ids"] == "P1|90.0,P2%2C odd|80.0"
