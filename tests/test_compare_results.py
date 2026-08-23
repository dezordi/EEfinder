"""Unit tests for the EE-vs-host decision rules."""

from __future__ import annotations

import pandas as pd
import pytest

from eefinder.compare_results import (
    ResolveSourcePriority,
    RULE_BITSCORE,
    RULE_COVERAGE,
    RULE_DENSITY,
    CompareResults,
    element_length,
)
from eefinder.filter_table import FILTERED_COLUMNS


def _row(qseqid, sseqid, bitscore, length, tag, evalue=1e-30, qstart=1, qend=300):
    """One filtered-table row (the schema FilterTable emits)."""
    return {
        "qseqid": qseqid,
        "sseqid": sseqid,
        "pident": 90.0,
        "length": length,
        "mismatch": 0,
        "gapopen": 0,
        "qstart": qstart,
        "qend": qend,
        "sstart": 1,
        "send": length,
        "evalue": evalue,
        "bitscore": bitscore,
        "sense": "pos",
        "bed_name": qseqid,
        "tag": tag,
    }


def _write(tmp_path, name, rows):
    path = tmp_path / name
    pd.DataFrame(rows, columns=FILTERED_COLUMNS).to_csv(path, sep="\t", index=False)
    return str(path)


def _kept(host_table):
    return pd.read_csv(f"{host_table}.concat.nr", sep="\t")


def test_element_length_from_the_bed_name():
    assert element_length("ctg1:100-400") == 300
    assert element_length("ctg1") == 0


class TestBitscoreRule:
    """The historical BLAST behaviour must stay byte-identical."""

    def test_element_survives_when_the_viral_hit_scores_higher(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-300", "PROT_A", 200, 100, "EE")]
        )
        host = _write(
            tmp_path, "host.tsv", [_row("ctg1:1-300", "HOST_A", 150, 100, "HOST")]
        )
        CompareResults(vir, host, rule=RULE_BITSCORE)
        assert list(_kept(host).sseqid) == ["PROT_A"]

    def test_element_is_dropped_when_the_host_hit_scores_higher(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-300", "PROT_A", 100, 100, "EE")]
        )
        host = _write(
            tmp_path, "host.tsv", [_row("ctg1:1-300", "HOST_A", 300, 100, "HOST")]
        )
        CompareResults(vir, host, rule=RULE_BITSCORE)
        assert _kept(host).empty

    def test_output_keeps_the_filtered_schema(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-300", "PROT_A", 200, 100, "EE")]
        )
        host = _write(tmp_path, "host.tsv", [])
        CompareResults(vir, host, rule=RULE_BITSCORE)
        assert list(_kept(host).columns) == FILTERED_COLUMNS


class TestDensityRule:
    """The symmetric profile-vs-profile comparison used by ``-md hmmer``."""

    def test_compares_bits_per_residue_not_raw_bitscore(self, tmp_path):
        # The host hit has the higher RAW score but a much longer alignment, so
        # per residue the viral hit wins: 200/50 = 4.0 vs 300/300 = 1.0.
        vir = _write(tmp_path, "vir.tsv", [_row("ctg1:1-900", "PROT_A", 200, 50, "EE")])
        host = _write(
            tmp_path, "host.tsv", [_row("ctg1:1-900", "HOST_A", 300, 300, "HOST")]
        )
        CompareResults(vir, host, rule=RULE_DENSITY)
        assert list(_kept(host).sseqid) == ["PROT_A"]

    def test_drops_the_element_when_the_host_density_wins(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-900", "PROT_A", 100, 100, "EE")]
        )
        host = _write(
            tmp_path, "host.tsv", [_row("ctg1:1-900", "HOST_A", 300, 100, "HOST")]
        )
        CompareResults(vir, host, rule=RULE_DENSITY)
        assert _kept(host).empty

    def test_keeps_the_element_when_there_is_no_host_hit(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-900", "PROT_A", 100, 100, "EE")]
        )
        host = _write(tmp_path, "host.tsv", [])
        CompareResults(vir, host, rule=RULE_DENSITY)
        assert list(_kept(host).sseqid) == ["PROT_A"]

    def test_margin_requires_the_viral_side_to_win_by_more(self, tmp_path):
        # 2.0 vs 1.9 bits/residue: kept without a margin, dropped with 10%.
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-900", "PROT_A", 200, 100, "EE")]
        )
        host = _write(
            tmp_path, "host.tsv", [_row("ctg1:1-900", "HOST_A", 190, 100, "HOST")]
        )
        CompareResults(vir, host, rule=RULE_DENSITY, margin=0.0)
        assert list(_kept(host).sseqid) == ["PROT_A"]
        CompareResults(vir, host, rule=RULE_DENSITY, margin=0.1)
        assert _kept(host).empty

    def test_best_viral_hit_per_element_is_reported(self, tmp_path):
        vir = _write(
            tmp_path,
            "vir.tsv",
            [
                _row("ctg1:1-900", "PROT_WEAK", 100, 100, "EE"),
                _row("ctg1:1-900", "PROT_BEST", 300, 100, "EE"),
            ],
        )
        host = _write(tmp_path, "host.tsv", [])
        CompareResults(vir, host, rule=RULE_DENSITY)
        assert list(_kept(host).sseqid) == ["PROT_BEST"]

    def test_helper_columns_do_not_leak_into_the_output(self, tmp_path):
        """GetFinalTaxonomy reads the taxonomy table by column position."""
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:1-900", "PROT_A", 200, 100, "EE")]
        )
        host = _write(
            tmp_path, "host.tsv", [_row("ctg1:1-900", "HOST_A", 100, 100, "HOST")]
        )
        CompareResults(vir, host, rule=RULE_DENSITY)
        assert list(_kept(host).columns) == FILTERED_COLUMNS


class TestCoverageRule:
    """The sequence-space fallback: significance + coverage, never a score mix."""

    def test_drops_an_element_a_significant_host_hit_covers(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:100-400", "PROT_A", 50, 100, "EE")]
        )
        host = _write(
            tmp_path,
            "host.tsv",
            [_row("ctg1:100-400", "HOST_A", 40, 100, "HOST", qstart=1, qend=280)],
        )
        CompareResults(vir, host, rule=RULE_COVERAGE, host_min_coverage=0.5)
        assert _kept(host).empty

    def test_keeps_an_element_with_only_a_short_host_hit(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:100-400", "PROT_A", 50, 100, "EE")]
        )
        host = _write(
            tmp_path,
            "host.tsv",
            [_row("ctg1:100-400", "HOST_A", 400, 30, "HOST", qstart=1, qend=90)],
        )
        CompareResults(vir, host, rule=RULE_COVERAGE, host_min_coverage=0.5)
        assert list(_kept(host).sseqid) == ["PROT_A"]

    def test_keeps_an_element_whose_host_hit_is_not_significant(self, tmp_path):
        vir = _write(
            tmp_path, "vir.tsv", [_row("ctg1:100-400", "PROT_A", 50, 100, "EE")]
        )
        host = _write(
            tmp_path,
            "host.tsv",
            [
                _row(
                    "ctg1:100-400",
                    "HOST_A",
                    400,
                    100,
                    "HOST",
                    evalue=0.5,
                    qstart=1,
                    qend=290,
                )
            ],
        )
        CompareResults(vir, host, rule=RULE_COVERAGE)
        assert list(_kept(host).sseqid) == ["PROT_A"]


def test_unknown_rule_is_rejected(tmp_path):
    vir = _write(tmp_path, "vir.tsv", [_row("ctg1:1-300", "PROT_A", 200, 100, "EE")])
    host = _write(tmp_path, "host.tsv", [])
    with pytest.raises(ValueError):
        CompareResults(vir, host, rule="hmm-vs-blast")


class TestResolveSourcePriority:
    """Which source names an element when several hit it with similar scores."""

    def _tables(self, tmp_path, all_rows, chosen):
        all_hits = _write(tmp_path, "host.tsv.concat", all_rows)
        ee_table = _write(tmp_path, "host.tsv.concat.nr", chosen)
        return ee_table, all_hits

    def test_a_preferred_source_within_the_margin_wins(self, tmp_path):
        # The deep VOGDB profile scores higher, but only by 2%.
        rows = [
            _row("ctg1:1-900", "VOGDB__VOG00001", 204, 100, "EE"),
            _row(
                "ctg1:1-900",
                "NCBIREFSEQ__Chuviridae__Glycoprotein__001",
                200,
                100,
                "EE",
            ),
        ]
        ee_table, all_hits = self._tables(tmp_path, rows, [rows[0]])
        step = ResolveSourcePriority(
            ee_table, all_hits, ["NCBIREFSEQ", "VOGDB"], margin=0.05
        )
        assert list(pd.read_csv(ee_table, sep="\t").sseqid) == [
            "NCBIREFSEQ__Chuviridae__Glycoprotein__001"
        ]
        assert step.reassigned == 1

    def test_outside_the_margin_the_score_still_decides(self, tmp_path):
        rows = [
            _row("ctg1:1-900", "VOGDB__VOG00001", 400, 100, "EE"),
            _row(
                "ctg1:1-900",
                "NCBIREFSEQ__Chuviridae__Glycoprotein__001",
                200,
                100,
                "EE",
            ),
        ]
        ee_table, all_hits = self._tables(tmp_path, rows, [rows[0]])
        step = ResolveSourcePriority(ee_table, all_hits, ["NCBIREFSEQ", "VOGDB"])
        assert list(pd.read_csv(ee_table, sep="\t").sseqid) == ["VOGDB__VOG00001"]
        assert step.reassigned == 0

    def test_host_hits_are_never_chosen(self, tmp_path):
        rows = [
            _row("ctg1:1-900", "VOGDB__VOG00001", 200, 100, "EE"),
            _row("ctg1:1-900", "HOST__Actin__001", 205, 100, "HOST"),
        ]
        ee_table, all_hits = self._tables(tmp_path, rows, [rows[0]])
        ResolveSourcePriority(ee_table, all_hits, ["HOST", "VOGDB"])
        assert list(pd.read_csv(ee_table, sep="\t").sseqid) == ["VOGDB__VOG00001"]

    def test_unlisted_sources_rank_last(self, tmp_path):
        rows = [
            _row("ctg1:1-900", "RVDB__FAM1", 202, 100, "EE"),
            _row("ctg1:1-900", "VOGDB__VOG00001", 200, 100, "EE"),
        ]
        ee_table, all_hits = self._tables(tmp_path, rows, [rows[0]])
        ResolveSourcePriority(ee_table, all_hits, ["VOGDB"])
        assert list(pd.read_csv(ee_table, sep="\t").sseqid) == ["VOGDB__VOG00001"]

    def test_the_schema_and_a_backup_are_preserved(self, tmp_path):
        rows = [
            _row("ctg1:1-900", "VOGDB__VOG00001", 204, 100, "EE"),
            _row("ctg1:1-900", "NCBIREFSEQ__X__Y__001", 200, 100, "EE"),
        ]
        ee_table, all_hits = self._tables(tmp_path, rows, [rows[0]])
        ResolveSourcePriority(ee_table, all_hits, ["NCBIREFSEQ", "VOGDB"])
        assert list(pd.read_csv(ee_table, sep="\t").columns) == FILTERED_COLUMNS
        backup = pd.read_csv(f"{ee_table}.unprioritised", sep="\t")
        assert list(backup.sseqid) == ["VOGDB__VOG00001"]
