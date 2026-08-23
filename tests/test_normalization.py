"""Unit tests for the consensus protein naming and the TE-like predicate."""

from __future__ import annotations

import pytest

from eefinder.normalization import (
    canonical_protein_names,
    consensus_protein_name,
    format_votes,
    is_transposon_like,
)

# The worked example of proposal.md §4.5.3: ten clustered chuvirus proteins whose
# names survive standardisation as three distinct labels, because the bundled map
# catches "glycoprotein" but neither the bare "gly" nor the bare "G".
CHUVIRUS_CLUSTER = ["Gly"] * 4 + ["Glycoprotein"] * 3 + ["G"] * 3


class TestConsensusProteinName:
    def test_plurality_takes_the_most_frequent_label(self):
        result = consensus_protein_name(CHUVIRUS_CLUSTER, mode="plurality")
        assert result.name == "Gly"
        assert result.agreement == pytest.approx(0.4)

    def test_canonical_first_prefers_a_canonical_label(self):
        result = consensus_protein_name(CHUVIRUS_CLUSTER, mode="canonical-first")
        assert result.name == "Glycoprotein"
        assert result.agreement == pytest.approx(0.3)

    def test_canonical_first_falls_back_to_plurality_below_the_fraction(self):
        # "Glycoprotein" holds 1/11 < 0.3, so plurality decides.
        names = ["Gly"] * 7 + ["Glycoprotein"] + ["G"] * 3
        assert consensus_protein_name(names).name == "Gly"

    def test_votes_are_reported_with_counts(self):
        result = consensus_protein_name(CHUVIRUS_CLUSTER, mode="plurality")
        assert result.votes == {"Gly": 4, "Glycoprotein": 3, "G": 3}
        assert format_votes(result.votes) == "Gly:4|Glycoprotein:3|G:3"

    def test_result_is_independent_of_input_order(self):
        shuffled = [
            "G",
            "Glycoprotein",
            "Gly",
            "G",
            "Gly",
            "Glycoprotein",
            "Gly",
            "G",
            "Gly",
            "Glycoprotein",
        ]
        for mode in ("plurality", "canonical-first"):
            assert (
                consensus_protein_name(shuffled, mode=mode).name
                == consensus_protein_name(CHUVIRUS_CLUSTER, mode=mode).name
            )

    def test_ties_prefer_the_canonical_name(self):
        result = consensus_protein_name(["RdRp", "Zzz"], mode="plurality")
        assert result.name == "RdRp"

    def test_unknown_labels_never_win(self):
        result = consensus_protein_name(["Unknown"] * 5 + ["Capsid Protein"])
        assert result.name == "Capsid Protein"

    def test_all_unknown_resolves_to_unknown(self):
        result = consensus_protein_name(["Unknown", "", None])
        assert result.name == "Unknown"
        assert result.votes == {}
        assert result.agreement == 0.0

    def test_unknown_vote_mode_is_rejected(self):
        with pytest.raises(ValueError):
            consensus_protein_name(["Capsid Protein"], mode="majority")

    def test_canonical_names_come_from_the_bundled_map(self):
        canonical = canonical_protein_names()
        assert "RdRp" in canonical
        assert "Glycoprotein" in canonical
        assert "Gly" not in canonical


class TestIsTransposonLike:
    @pytest.mark.parametrize(
        "product",
        [
            "reverse transcriptase",
            "endonuclease-reverse transcriptase",
            "RNA-directed DNA polymerase",
            "integrase core domain protein",
            "transposase",
            "retrotransposon protein",
            "Gag",
            "gag-pol polyprotein",
        ],
    )
    def test_transposon_like_products_are_flagged(self, product):
        assert is_transposon_like(product)

    @pytest.mark.parametrize(
        "product",
        [
            "Glycoprotein",
            "DNA-directed RNA polymerase subunit beta",
            "cytochrome c oxidase subunit 1",
            "igE-binding protein-like",
            "",
        ],
    )
    def test_ordinary_host_products_are_not_flagged(self, product):
        assert not is_transposon_like(product)
