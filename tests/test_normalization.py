"""Unit tests for the data-driven protein-name rules in eefinder.normalization."""

from __future__ import annotations
import pytest
from eefinder.get_databases import molecule_type_for_family, strip_placeholder_prefix
from eefinder.normalization import (
    _PROTEIN_RULES,
    _REWRITES,
    _TYPO_CORRECTIONS,
    lineage_taxa,
    standardize_protein,
)

CHUVIRIDAE = "r__Riboviria;k__Orthornavirae;p__Negarnaviricota;c__Monjiviricetes;o__Jingchuvirales;f__Chuviridae;g__Mivirus;s__Mivirus sp."


def test_the_bundled_rule_file_loads():
    assert _TYPO_CORRECTIONS and _REWRITES and _PROTEIN_RULES


def test_lineage_taxa_reads_a_lineage_string():
    taxa = lineage_taxa(CHUVIRIDAE)

    assert "chuviridae" in taxa
    assert "jingchuvirales" in taxa


def test_lineage_taxa_accepts_a_bare_taxon_name():
    """A caller holding only a family can pass it directly."""
    assert lineage_taxa("Chuviridae") == frozenset({"chuviridae"})


def test_lineage_taxa_of_nothing():
    assert lineage_taxa("") == frozenset()


@pytest.mark.parametrize(
    "raw", ["G", "G protein", "Envelope protein", "S protein", "Gp4"]
)
def test_chuviridae_glycoprotein_synonyms(raw):
    assert standardize_protein(raw, "ssRNA(-)", taxonomy="Chuviridae") == "Glycoprotein"


@pytest.mark.parametrize(
    "raw",
    ["Coat protein", "Capsid protein", "Nucleocoprotein", "Nucleprotein", "N protein"],
)
def test_chuviridae_nucleocapsid_synonyms(raw):
    assert standardize_protein(raw, "ssRNA(-)", taxonomy="Chuviridae") == (
        "Nucleocapsid Protein"
    )


def test_a_family_rule_wins_over_the_group_wide_one():
    """Chuviridae has no separate coat protein; (+)RNA viruses do."""
    assert standardize_protein("Coat protein", "ssRNA(-)", taxonomy="Chuviridae") == (
        "Nucleocapsid Protein"
    )
    assert standardize_protein("Coat protein", "ssRNA(+)") == "Capsid Protein"


def test_a_family_rule_does_not_leak_into_another_family():
    assert standardize_protein(
        "Membrane protein", "ssRNA(-)", taxonomy="Rhabdoviridae"
    ) == ("Matrix protein")
    assert standardize_protein(
        "Membrane protein", "ssRNA(-)", taxonomy="Paramyxoviridae"
    ) == ("Membrane protein")


def test_a_family_rule_needs_the_taxonomy_to_fire():
    """Without a lineage only the unrestricted rules apply."""
    assert standardize_protein("Gp4", "ssRNA(-)") == "Gp4"
    assert standardize_protein("Gp4", "ssRNA(-)", taxonomy=CHUVIRIDAE) == "Glycoprotein"


def test_a_rule_scoped_to_a_family_matches_any_rank_of_the_lineage():
    assert standardize_protein("S protein", "ssRNA(-)", taxonomy=CHUVIRIDAE) == (
        "Glycoprotein"
    )


def test_the_earliest_match_wins_over_file_order():
    """A product string leads with the protein's name."""
    assert standardize_protein("nucleocapsid phosphoprotein", "ssRNA(-)") == (
        "Nucleocapsid Protein"
    )


def test_exact_beats_contains():
    assert standardize_protein("RNA polymerase sigma factor", "dsDNA") == (
        "RNA polymerase sigma factor"
    )


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("ORF3 protein", "ORF3"),
        ("ORF12 protein", "ORF12"),
        ("14 kDa protein unknown protein", "14 kDa protein"),
        ("52 kDa protein unknown in function protein", "52 kDa protein"),
        ("Matrix putative protein", "Matrix protein"),
        ("U1 putative protein", "U1 protein"),
    ],
)
def test_shape_rewrites(raw, expected):
    assert standardize_protein(raw, "ssRNA(-)", taxonomy="Rhabdoviridae") == expected


def test_a_bare_unknown_protein_is_kept():
    """The tail rule needs something before it, so this is not emptied out."""
    assert standardize_protein("Unknown protein", "dsDNA") == "Unknown protein"


def test_a_reversed_designation_is_turned_round():
    assert standardize_protein("Protein L", "ssRNA(-)") == "RdRp"
    assert standardize_protein("Protein M1", "ssRNA(-)", taxonomy="Rhabdoviridae") == (
        "Matrix protein"
    )


def test_a_long_tail_is_not_mistaken_for_a_designation():
    assert standardize_protein("Protein of unknown function", "dsDNA") == (
        "Protein of unknown function"
    )


def test_a_word_fused_onto_protein_is_split():
    assert standardize_protein(
        "G proteinGFP fusion protein", "ssRNA(-)", taxonomy="Rhabdoviridae"
    ) == ("Glycoprotein")


def test_the_fused_word_split_leaves_a_plural_alone():
    assert standardize_protein("Structural proteins", "dsDNA") == "Structural proteins"


@pytest.mark.parametrize(
    "raw",
    [
        "Glyciprotein",
        "Glycotrotein",
        "Glycpprotein",
        "Glyocprotein",
        "Glyoprotein",
        "Gylcoprotein",
    ],
)
def test_glycoprotein_misspellings(raw):
    assert standardize_protein(raw, "ssRNA(-)") == "Glycoprotein"


@pytest.mark.parametrize(
    "raw",
    [
        "Nucleprotein",
        "Nucleoproteinn",
        "Nucleotprotein",
        "Nucloprotein",
        "Nulcleoprotein",
        "Nucleicapsid",
    ],
)
def test_nucleocapsid_misspellings(raw):
    assert standardize_protein(raw, "ssRNA(-)") == "Nucleocapsid Protein"


def test_a_misspelled_hypothetical_is_still_dropped():
    assert standardize_protein("Hypothtetical protein", "ssRNA(-)") == "Unknown"


def test_a_typo_that_prefixes_the_correct_word_leaves_it_alone():
    assert standardize_protein("Membrane fusion protein", "dsDNA") == (
        "Membrane fusion protein"
    )


def test_strip_placeholder_prefix():
    assert strip_placeholder_prefix("unclassified Chuviridae") == "Chuviridae"
    assert strip_placeholder_prefix("Chuviridae") == "Chuviridae"


def test_a_placeholder_family_keeps_its_molecule_type():
    """Losing it would silently disable every molecule-type-scoped rule."""
    assert molecule_type_for_family("unclassified Chuviridae") == (
        molecule_type_for_family("Chuviridae")
    )


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("RNA-dependant RNA polymelase", "RdRp"),
        ("RNA-depedent RNA polymerase", "RdRp"),
        ("Nuclecapsid protein", "Nucleocapsid Protein"),
        ("Neucleocapsid protein", "Nucleocapsid Protein"),
        ("Phosphorprotein", "Phosphoprotein"),
        ("Glicoprotein precursor", "Glycoprotein"),
    ],
)
def test_misspellings_found_in_a_whole_refseq_build(raw, expected):
    assert standardize_protein(raw, "ssRNA(-)") == expected


@pytest.mark.parametrize("raw", ["Revrse transcriptase", "Reverse transciptase"])
def test_reverse_transcriptase_misspellings(raw):
    assert standardize_protein(raw, "ssRNA-RT") == "Reverse Transcriptase"


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("Base plate hub", "Baseplate Protein"),
        ("Tapemeasure protein", "Tape Measure Protein"),
        ("Tail length tapemeasure protein", "Tape Measure Protein"),
        ("Majorcapsid protein", "Major Capsid Protein"),
        ("Major capid protein", "Major Capsid Protein"),
        ("Terminase large subunuit", "Terminase Large Subunit"),
    ],
)
def test_split_and_fused_word_misspellings(raw, expected):
    assert standardize_protein(raw, "dsDNA") == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("Tail fibers protein", "Tail Fiber Protein"),
        ("Tail fibers", "Tail Fiber Protein"),
    ],
)
def test_plural_variants(raw, expected):
    assert standardize_protein(raw, "dsDNA") == expected


@pytest.mark.parametrize(
    "raw", ["DnaB helicase", "DnaB-like helicase", "DnaG-like primase"]
)
def test_a_named_replication_protein_is_not_read_as_a_dna_misspelling(raw):
    """DnaB/DnaG name proteins; they are not a mangled "DNA"."""
    assert standardize_protein(raw, "dsDNA") == raw
