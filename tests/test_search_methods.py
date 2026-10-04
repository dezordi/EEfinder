"""Unit tests for eefinder.search_methods."""

from __future__ import annotations
import pytest
from eefinder.models import ModelError
from eefinder.search_methods import (
    SEARCH_METHODS,
    BlastxSearch,
    DiamondSearch,
    PredictedProteinSearch,
    SearchMethod,
    register_search_method,
    resolve_search_method,
)


def test_the_builtin_methods_are_registered():
    assert {"blastx", "diamond", "predicted"} <= set(SEARCH_METHODS)


@pytest.mark.parametrize(
    "mode,translation,expected",
    [
        ("blastx", "default", "blastx"),
        ("very-sensitive", "default", "diamond"),
        ("fast", "default", "diamond"),
        ("blastx", "gv", "predicted"),
        ("blastx", "gv-rv", "predicted"),
        ("more-sensitive", "rv", "predicted"),
    ],
)
def test_resolve_search_method(mode, translation, expected):
    assert resolve_search_method(mode, translation).name == expected


def test_diamond_keeps_the_requested_sensitivity():
    assert resolve_search_method("ultra-sensitive", "default").mode == "ultra-sensitive"


def test_predicted_requires_cd_hit_only_for_gv_rv():
    assert "cd-hit" not in PredictedProteinSearch("gv").requires
    assert "cd-hit" in PredictedProteinSearch("gv-rv").requires


def test_predicted_needs_diamond_when_the_aligner_is_diamond():
    assert PredictedProteinSearch("gv", aligner="sensitive").requires == ("diamond",)
    assert PredictedProteinSearch("gv", aligner="blastx").requires == ("blastp",)


def test_predicted_rejects_an_unknown_translation_method():
    with pytest.raises(ValueError, match="unknown translation method"):
        PredictedProteinSearch("nope")


def test_check_available_names_the_missing_binary():
    class Missing(SearchMethod):
        name = "missing-tool"
        requires = ("definitely-not-on-path",)

        def search(self, query, database, threads, out_table):  # pragma: no cover
            raise AssertionError("should not be reached")

    with pytest.raises(FileNotFoundError, match="definitely-not-on-path"):
        Missing().check_available()


def test_run_rejects_a_method_that_violates_the_hit_model(tmp_path):
    class Wrong(SearchMethod):
        name = "wrong-shape"

        def search(self, query, database, threads, out_table):
            # Three fields instead of the twelve outfmt6 requires.
            with open(out_table, "w") as handle:
                handle.write("q\ts\t30.0\n")

    out = tmp_path / "hits.tsv"
    with pytest.raises(ModelError, match="violated its output model"):
        Wrong().run("query.fa", "db.fa", 1, str(out))


def test_run_accepts_a_conforming_method(tmp_path):
    class Fine(SearchMethod):
        name = "fine-shape"

        def search(self, query, database, threads, out_table):
            with open(out_table, "w") as handle:
                handle.write("\t".join(["x"] * 12) + "\n")

    out = tmp_path / "hits.tsv"
    assert Fine().run("query.fa", "db.fa", 1, str(out)) == str(out)


def test_register_search_method_rejects_a_nameless_class():
    with pytest.raises(ValueError, match="non-empty 'name'"):

        @register_search_method
        class Nameless(SearchMethod):
            def search(self, query, database, threads, out_table):  # pragma: no cover
                raise AssertionError("should not be reached")


def test_describe_mentions_the_required_binaries():
    assert "blastx" in BlastxSearch().describe()
    assert "diamond" in DiamondSearch("fast").describe()
