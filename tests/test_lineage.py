"""Unit tests for eefinder.lineage."""

from __future__ import annotations
import pytest
from eefinder.lineage import (
    RANKS,
    UNKNOWN,
    format_lineage,
    is_ranked,
    lca,
    lineage_from_nodes,
    merge_key,
    parse_lineage,
    rank_name,
    trim_lineage,
)

FULL = {
    "realm": "Riboviria",
    "kingdom": "Orthornavirae",
    "phylum": "Pisuviricota",
    "class": "Pisoniviricetes",
    "order": "Nidovirales",
    "family": "Coronaviridae",
    "genus": "Betacoronavirus",
    "species": "Betacoronavirus pandemicum",
}


def test_a_complete_lineage_renders_every_rank():
    assert format_lineage(FULL) == (
        "r__Riboviria;k__Orthornavirae;p__Pisuviricota;c__Pisoniviricetes;"
        "o__Nidovirales;f__Coronaviridae;g__Betacoronavirus;"
        "s__Betacoronavirus pandemicum"
    )


def test_every_rank_is_always_present():
    rendered = format_lineage({"family": "Coronaviridae"})

    assert rendered.count(";") == len(RANKS) - 1
    for _, prefix in RANKS:
        assert f"{prefix}__" in rendered


def test_a_missing_rank_becomes_unk():
    assert format_lineage({"family": "Coronaviridae"}).startswith(f"r__{UNKNOWN};")


def test_a_literal_unknown_becomes_unk():
    assert f"f__{UNKNOWN}" in format_lineage({"family": "Unknown"})


def test_an_empty_name_becomes_unk():
    assert f"g__{UNKNOWN}" in format_lineage({"genus": "   "})


def test_parse_round_trips_the_named_ranks():
    assert parse_lineage(format_lineage(FULL)) == FULL


def test_parse_treats_unk_as_absent():
    parsed = parse_lineage(format_lineage({"family": "Coronaviridae"}))

    assert parsed == {"family": "Coronaviridae"}


def test_parse_ignores_malformed_fields():
    assert parse_lineage("nonsense;f__Coronaviridae;;x__y") == {
        "family": "Coronaviridae"
    }


def test_parse_handles_an_empty_value():
    assert parse_lineage("") == {}


def test_rank_name_returns_the_taxon():
    assert rank_name(format_lineage(FULL), "order") == "Nidovirales"


def test_rank_name_is_empty_for_unk():
    assert rank_name(format_lineage({"family": "X"}), "order") == ""


@pytest.mark.parametrize("level", ["family", "genus"])
def test_merge_key_is_the_taxon_when_ranked(level):
    assert merge_key(format_lineage(FULL), level) == FULL[level]


def test_merge_key_falls_back_to_the_whole_lineage():
    unranked = format_lineage(
        {"realm": "Duplodnaviria", "class": "Caudoviricetes", "species": "phage sp."}
    )

    assert merge_key(unranked, "family") == unranked


def test_unranked_viruses_of_different_orders_do_not_share_a_key():
    """The point of carrying the lineage: unknowns stay distinguishable."""
    one = format_lineage(
        {"realm": "Duplodnaviria", "order": "Crassvirales", "species": "phage a"}
    )
    two = format_lineage(
        {"realm": "Duplodnaviria", "order": "Tubulavirales", "species": "phage b"}
    )

    assert merge_key(one, "family") != merge_key(two, "family")


def test_unranked_viruses_of_the_same_lineage_share_a_key():
    lineage = format_lineage(
        {"realm": "Duplodnaviria", "order": "Crassvirales", "species": "phage a"}
    )

    assert merge_key(lineage, "family") == merge_key(lineage, "family")


def test_is_ranked():
    full = format_lineage(FULL)

    assert is_ranked(full, "family")
    assert not is_ranked(format_lineage({"species": "x"}), "family")


def test_lineage_from_nodes_labels_the_entries():
    nodes = [
        {"name": "Viruses", "taxId": 10239},
        {"name": "Riboviria", "taxId": 2559587},
        {"name": "Coronaviridae", "taxId": 11118},
    ]
    ranks = {10239: "acellular_root", 2559587: "realm", 11118: "family"}

    rendered = lineage_from_nodes(nodes, ranks, species="Some virus")

    assert rank_name(rendered, "realm") == "Riboviria"
    assert rank_name(rendered, "family") == "Coronaviridae"
    assert rank_name(rendered, "species") == "Some virus"
    assert "Viruses" not in rendered


def test_lineage_from_nodes_without_ranks_is_all_unk():
    rendered = lineage_from_nodes([{"name": "Viruses", "taxId": 10239}], {})

    assert set(parse_lineage(rendered)) == set()


def test_trim_lineage_drops_the_unresolved_tail():
    assert trim_lineage(format_lineage({"realm": "Riboviria", "family": "Fam"})) == (
        "r__Riboviria;k__Unk;p__Unk;c__Unk;o__Unk;f__Fam"
    )


def test_trim_lineage_keeps_an_internal_gap():
    trimmed = trim_lineage(format_lineage({"realm": "Riboviria", "species": "Sp"}))

    assert trimmed.endswith("s__Sp")
    assert "f__Unk" in trimmed


def test_trim_lineage_of_an_empty_lineage():
    assert trim_lineage(format_lineage({})) == ""


def test_lca_of_one_lineage_is_itself():
    assert lca([format_lineage(FULL)]) == format_lineage(FULL)


def test_lca_of_two_species_is_their_genus():
    one = format_lineage({**FULL, "species": "Betacoronavirus a"})
    two = format_lineage({**FULL, "species": "Betacoronavirus b"})

    assert rank_name(lca([one, two]), "genus") == FULL["genus"]
    assert rank_name(lca([one, two]), "species") == ""


def test_lca_of_two_genera_is_their_family():
    one = format_lineage({**FULL, "genus": "Alphacoronavirus", "species": "a"})
    two = format_lineage({**FULL, "genus": "Betacoronavirus", "species": "b"})
    resolved = lca([one, two])

    assert rank_name(resolved, "family") == FULL["family"]
    assert rank_name(resolved, "genus") == ""
    assert resolved.endswith(f"f__{FULL['family']}")


def test_lca_stops_where_a_rank_is_unnamed_in_one_lineage():
    """A lineage with no genus cannot vouch for the other's genus."""
    one = format_lineage({"family": "Fam", "genus": "Gen", "species": "a"})
    two = format_lineage({"family": "Fam", "species": "b"})

    assert rank_name(lca([one, two]), "family") == "Fam"
    assert rank_name(lca([one, two]), "genus") == ""


def test_lca_keeps_deeper_ranks_when_no_lineage_names_the_upper_ones():
    """Viral lineages often lack a realm; a shared family must still survive."""
    one = format_lineage({"family": "Fam", "genus": "Gen", "species": "a"})
    two = format_lineage({"family": "Fam", "genus": "Gen", "species": "b"})

    assert lca([one, two]) == "r__Unk;k__Unk;p__Unk;c__Unk;o__Unk;f__Fam;g__Gen"


def test_lca_of_unrelated_lineages_is_empty():
    one = format_lineage({"realm": "Riboviria", "family": "A"})
    two = format_lineage({"realm": "Monodnaviria", "family": "B"})

    assert lca([one, two]) == ""


def test_lca_ignores_empty_values():
    full = format_lineage(FULL)

    assert lca([full, "", None]) == full


def test_lca_of_nothing_is_empty():
    assert lca([]) == ""
