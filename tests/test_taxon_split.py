"""Unit tests for download splitting (no network: the CLI is faked)."""

import json
import subprocess
import pytest
from eefinder import taxon_exclusion, taxon_split
from eefinder.taxon_split import NO_SPLIT, SkippedTaxon, plan_split, taxa_at_rank

TREE = {
    1: {
        "name": "Root",
        "rank": "acellular_root",
        "parents": [],
        "children": [10, 20, 30],
    },
    10: {
        "name": "OrderA",
        "rank": "order",
        "parents": [1],
        "children": [100, 200, 300],
    },
    100: {
        "name": "FamA",
        "rank": "family",
        "parents": [1, 10],
        "children": [1000],
        "asm": 7,
    },
    200: {
        "name": "FamB",
        "rank": "family",
        "parents": [1, 10],
        "children": [],
        "asm": 3,
    },
    1000: {"name": "SpA", "rank": "species", "parents": [1, 10, 100], "children": []},
    20: {"name": "OrderB", "rank": "order", "parents": [1], "children": [201]},
    201: {"name": "UnclassB", "rank": "", "parents": [1, 20], "children": [], "asm": 5},
    30: {"name": "FamC", "rank": "family", "parents": [1], "children": [], "asm": 11},
    300: {"name": "FamEmpty", "rank": "family", "parents": [1, 10], "children": []},
}
BY_NAME = {node["name"].lower(): tax_id for tax_id, node in TREE.items()}


def _record(tax_id):
    node = TREE[tax_id]
    return {
        "taxonomy": {
            "tax_id": tax_id,
            "current_scientific_name": {"name": node["name"]},
            "rank": node["rank"].upper(),
            "parents": node["parents"],
            "children": node["children"],
            "counts": [{"type": "COUNT_TYPE_ASSEMBLY", "count": node.get("asm", 0)}],
        }
    }


def _descendants(tax_id):
    out = []
    for child in TREE[tax_id]["children"]:
        out.append(child)
        out.extend(_descendants(child))
    return out


def _fake_datasets(monkeypatch, calls=None):
    """Replace subprocess.run with a lookup into TREE, honouring --rank."""

    def fake_run(args, **kwargs):
        wanted = args[args.index("taxon") + 1].split(",")
        rank = args[args.index("--rank") + 1] if "--rank" in args else None
        if calls is not None:
            calls.append((tuple(wanted), rank))
        ids = []
        for item in wanted:
            tax_id = BY_NAME.get(item.lower())
            if tax_id is None:
                try:
                    tax_id = int(item)
                except ValueError:
                    tax_id = None
            if tax_id is None or tax_id not in TREE:
                return subprocess.CompletedProcess(
                    args, 1, "", f"Error: '{item}' is not recognized."
                )
            if rank:
                ids.extend(d for d in _descendants(tax_id) if TREE[d]["rank"] == rank)
            else:
                ids.append(tax_id)
        lines = [json.dumps(_record(i)) for i in ids]
        return subprocess.CompletedProcess(args, 0, "\n".join(lines) + "\n", "")

    monkeypatch.setattr(taxon_exclusion.subprocess, "run", fake_run)
    monkeypatch.setattr(taxon_split.subprocess, "run", fake_run)


def test_taxa_at_rank_lists_only_that_rank(monkeypatch):
    _fake_datasets(monkeypatch)

    found = taxa_at_rank("1", "family")

    assert set(found) == {100, 200, 30, 300}
    assert all(node.rank == "family" for node in found.values())


def test_taxa_at_rank_carries_the_assembly_count(monkeypatch):
    _fake_datasets(monkeypatch)

    assert taxa_at_rank("1", "family")[30].assembly_count == 11


def test_plan_splits_a_broad_taxon_into_families(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", "family")

    assert plan.level == "family"
    assert plan.split
    assert [node.tax_id for node in plan.taxa] == [30, 100, 200]


def test_plan_reports_subtrees_with_no_family(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", "family")

    by_reason = {}
    for entry in plan.skipped:
        by_reason.setdefault(entry.reason, []).append(entry.tax_id)
    assert by_reason["no family in its lineage"] == [20]


def test_skipped_assemblies_sums_the_subtrees(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", "family")

    assert plan.skipped_assemblies == sum(e.assembly_count for e in plan.skipped)


def test_a_taxon_already_at_the_level_is_not_split(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("FamA", "family")

    assert not plan.split
    assert [node.tax_id for node in plan.taxa] == [100]
    assert plan.skipped == ()


def test_no_split_requests_the_root_once(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", NO_SPLIT)

    assert plan.level == NO_SPLIT
    assert [node.tax_id for node in plan.taxa] == [1]


def test_a_taxon_with_no_family_below_is_not_split(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("OrderB", "family")

    assert not plan.split
    assert [node.tax_id for node in plan.taxa] == [20]


def test_unknown_level_is_rejected(monkeypatch):
    _fake_datasets(monkeypatch)

    with pytest.raises(ValueError, match="unknown split level"):
        plan_split("1", "kingdom")


def test_unknown_root_is_reported(monkeypatch):
    _fake_datasets(monkeypatch)

    with pytest.raises(RuntimeError, match="9999"):
        plan_split("9999", "family")


def test_the_rank_listing_is_one_call(monkeypatch):
    calls = []
    _fake_datasets(monkeypatch, calls)

    plan_split("1", "family")

    ranked = [c for c in calls if c[1] == "family"]
    assert len(ranked) == 1


def test_skipped_taxon_fields():
    entry = SkippedTaxon(
        tax_id=1, name="n", rank="order", assembly_count=2, reason="why"
    )

    assert entry.tax_id == 1 and entry.reason == "why"


def test_a_family_with_no_records_is_not_requested(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", "family")

    assert 300 not in [node.tax_id for node in plan.taxa]


def test_a_family_with_no_records_is_reported(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", "family")

    empty = [e for e in plan.skipped if e.reason == "no records in NCBI"]
    assert [e.tax_id for e in empty] == [300]
    assert empty[0].assembly_count == 0


def test_skipping_empty_taxa_does_not_inflate_the_skipped_assemblies(monkeypatch):
    _fake_datasets(monkeypatch)

    plan = plan_split("1", "family")

    assert plan.skipped_assemblies == sum(
        e.assembly_count for e in plan.skipped if e.reason != "no records in NCBI"
    )


def test_all_counts_zero_requests_everything(monkeypatch):
    """A missing count field must not look like an empty taxonomy."""
    stripped = {
        tax_id: {k: v for k, v in node.items() if k != "asm"}
        for tax_id, node in TREE.items()
    }
    monkeypatch.setattr(taxon_split, "TREE", stripped, raising=False)
    original = dict(TREE)
    TREE.clear()
    TREE.update(stripped)
    try:
        _fake_datasets(monkeypatch)
        plan = plan_split("1", "family")
        assert [node.tax_id for node in plan.taxa] == [30, 100, 200, 300]
    finally:
        TREE.clear()
        TREE.update(original)
