"""Unit tests for eefinder.tag_elements."""

from __future__ import annotations

import pandas as pd

from eefinder.tag_elements import TagElements, _list_to_string


def test_list_to_string():
    assert _list_to_string(["a", "b"]) == "a,b"
    assert _list_to_string([]) == ""


def _write_tax(path):
    pd.DataFrame(
        {
            "element_id": ["ctg1:100-200", "ctg1:250-300", "ctg2:100-200"],
            "sense": ["pos", "pos", "pos"],
            "protein_ids": ["P1|30.0", "P2|40.0;P3|50.0", "P4|20.0"],
            "protein_products": ["prot", "prot", "prot"],
            "molecule_type": ["ssRNA", "ssRNA", "ssRNA"],
            "taxonomy": [
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamA;g__GenA;s__spA",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamB;g__GenB;s__spB",
                "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__FamC;g__GenC;s__spC",
            ],
            "host": ["h", "h", "h"],
        }
    ).to_csv(path, sep="\t", index=False)


def test_tag_elements_flags_overlaps(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(tax)

    TagElements(str(tax))

    df = pd.read_csv(tax, sep="\t").set_index("element_id")
    assert df.loc["ctg1:100-200", "tag"] == "overlaped"
    assert df.loc["ctg1:250-300", "tag"] == "overlaped"
    assert df.loc["ctg2:100-200", "tag"] == "unique"
    assert "ctg1:250-300" in df.loc["ctg1:100-200", "overlaped_element_id"]


def test_tag_elements_average_pident(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(tax)

    TagElements(str(tax))

    df = pd.read_csv(tax, sep="\t").set_index("element_id")
    assert df.loc["ctg1:100-200", "average_pident"] == 30.0
    assert df.loc["ctg1:250-300", "average_pident"] == 45.0
    assert df.loc["ctg2:100-200", "average_pident"] == 20.0
