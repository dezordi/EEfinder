"""Unit tests for the hmmsearch table conversion (no binaries required)."""

from __future__ import annotations

import textwrap

import pandas as pd
import pytest

from eefinder.filter_table import FILTERED_COLUMNS, OUTFMT6_COLUMNS, FilterTable
from eefinder.hmm_analysis import (
    domtbl_to_aa_table,
    has_ga_thresholds,
    parse_domtbl,
)
from eefinder.translation import traceback

# Two real hmmsearch --domtblout rows (HMMER 3.4): a plus-strand segment matching
# a profile over its full length, and a minus-strand segment matching part of it.
DOMTBL = """\
#                                                                    --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
ctg1__f1__00001      -            145 FAM__Glycoprotein__001 -           145   3.6e-92  311.7   0.7   1   1     4e-92     4e-92  311.5   0.7     1   145     1   145     1   145 0.98 seg
ctg1__f4__00002      -             95 FAM__Glycoprotein__001 -           145   5.7e-55  191.0   0.0   1   1   6.3e-55   6.3e-55  190.9   0.0     1    30    11    40    11    40 0.90 seg
"""


@pytest.fixture
def domtbl(tmp_path):
    path = tmp_path / "hits.domtbl"
    path.write_text(DOMTBL)
    return path


class TestParseDomtbl:
    def test_fields_are_mapped_from_the_right_columns(self, domtbl):
        hits = list(parse_domtbl(str(domtbl)))
        assert len(hits) == 2
        first = hits[0]
        assert first.target == "ctg1__f1__00001"
        assert first.profile == "FAM__Glycoprotein__001"
        assert first.profile_length == 145
        assert first.evalue == pytest.approx(4e-92)
        assert first.bitscore == pytest.approx(311.5)  # the domain score, not 311.7
        assert (first.hmm_from, first.hmm_to) == (1, 145)
        assert (first.ali_from, first.ali_to) == (1, 145)
        assert first.accuracy == pytest.approx(0.98)

    def test_coverage_is_the_profile_fraction_aligned(self, domtbl):
        hits = list(parse_domtbl(str(domtbl)))
        assert hits[0].coverage == pytest.approx(1.0)
        assert hits[1].coverage == pytest.approx(30 / 145)

    def test_malformed_rows_are_skipped(self, tmp_path):
        path = tmp_path / "broken.domtbl"
        path.write_text(DOMTBL + "not a domtbl row\n")
        assert len(list(parse_domtbl(str(path)))) == 2


class TestDomtblToAaTable:
    def test_emits_the_twelve_outfmt6_columns(self, domtbl, tmp_path):
        out = tmp_path / "hits.aa.tsv"
        assert domtbl_to_aa_table(str(domtbl), str(out)) == 2

        table = pd.read_csv(out, sep="\t", header=None, names=OUTFMT6_COLUMNS)
        assert list(table.qseqid) == ["ctg1__f1__00001", "ctg1__f4__00002"]
        assert set(table.sseqid) == {"FAM__Glycoprotein__001"}
        # pident carries the posterior accuracy (acc x 100), not identity.
        assert list(table.pident) == [98.0, 90.0]
        assert list(table.length) == [145, 30]
        assert list(table.qstart) == [1, 11]
        assert list(table.sstart) == [1, 1]
        assert list(table.bitscore) == [311.5, 190.9]

    def test_min_coverage_drops_partial_profile_hits(self, domtbl, tmp_path):
        out = tmp_path / "covered.tsv"
        assert domtbl_to_aa_table(str(domtbl), str(out), min_coverage=0.5) == 1
        table = pd.read_csv(out, sep="\t", header=None, names=OUTFMT6_COLUMNS)
        assert list(table.qseqid) == ["ctg1__f1__00001"]

    def test_append_merges_several_databases(self, domtbl, tmp_path):
        out = tmp_path / "merged.tsv"
        domtbl_to_aa_table(str(domtbl), str(out))
        domtbl_to_aa_table(str(domtbl), str(out), append=True)
        assert sum(1 for _ in open(out)) == 4

    def test_output_is_consumable_by_filter_table(self, domtbl, tmp_path):
        """The emitted table must flow into the pipeline unchanged."""
        coords = tmp_path / "coords.tsv"
        coords.write_text(
            "protein_id\tcontig\tstart\tend\tstrand\ttool\n"
            "ctg1__f1__00001\tctg1\t1\t435\t+\tsixframe\n"
            "ctg1__f4__00002\tctg1\t500\t784\t-\tsixframe\n"
        )
        aa_table = tmp_path / "hits.aa.tsv"
        domtbl_to_aa_table(str(domtbl), str(aa_table))
        nt_table = tmp_path / "hits.blastx"
        traceback(str(aa_table), str(coords), str(nt_table))

        traced = pd.read_csv(nt_table, sep="\t", header=None, names=OUTFMT6_COLUMNS)
        assert set(traced.qseqid) == {"ctg1"}
        # The plus-strand hit keeps qstart < qend and the minus-strand hit gets
        # qstart > qend, which is how FilterTable infers the sense.
        assert traced.loc[0, "qstart"] < traced.loc[0, "qend"]
        assert traced.loc[1, "qstart"] > traced.loc[1, "qend"]

        FilterTable(str(nt_table), 100, "EE", str(tmp_path))
        filtered = pd.read_csv(f"{nt_table}.filtred", sep="\t")
        assert list(filtered.columns) == FILTERED_COLUMNS
        # Only the 145 aa hit survives: the 30 aa one is below MIN_HIT_LENGTH.
        assert list(filtered.qseqid) == ["ctg1"]
        assert list(filtered.sense) == ["pos"]


class TestHasGaThresholds:
    def _write_hmm(self, path, with_ga):
        body = textwrap.dedent("""\
            HMMER3/f [3.4 | Aug 2023]
            NAME  PROFILE_A
            LENG  10
            {ga}HMM          A        C
            //
            HMMER3/f [3.4 | Aug 2023]
            NAME  PROFILE_B
            LENG  10
            HMM          A        C
            //
            """).format(ga="GA    25.00 25.00;\n" if with_ga else "")
        path.write_text(body)
        return str(path)

    def test_false_when_a_profile_lacks_its_threshold(self, tmp_path):
        assert not has_ga_thresholds(self._write_hmm(tmp_path / "a.hmm", True))

    def test_false_for_a_database_without_any(self, tmp_path):
        assert not has_ga_thresholds(self._write_hmm(tmp_path / "b.hmm", False))


PROFILE_TEMPLATE = """\
HMMER3/f [3.4 | Aug 2023]
NAME  {name}
LENG  76
ALPH  amino
NSEQ  2
EFFN  0.503906
CKSUM 3677150915
{ga}STATS LOCAL MSV       -9.1651  0.71860
STATS LOCAL VITERBI   -9.6832  0.71860
STATS LOCAL FORWARD   -4.2617  0.71860
HMM          A        C
//
"""


def _database(tmp_path, names, ga=""):
    path = tmp_path / "db.hmm"
    path.write_text(
        "".join(PROFILE_TEMPLATE.format(name=name, ga=ga) for name in names)
    )
    return str(path)


class TestScoreForEvalue:
    def test_matches_the_hmmer_forward_tail(self):
        from eefinder.hmm_analysis import score_for_evalue

        # E = Z * exp(-lambda * (S - tau)) with the STATS above.
        assert score_for_evalue(-4.2617, 0.71860, 1e-5, 1e6) == pytest.approx(
            30.98, abs=0.01
        )

    def test_a_stricter_evalue_needs_a_higher_score(self):
        from eefinder.hmm_analysis import score_for_evalue

        assert score_for_evalue(-4.2617, 0.7186, 1e-10, 1e6) > score_for_evalue(
            -4.2617, 0.7186, 1e-5, 1e6
        )

    def test_degenerate_stats_do_not_raise(self):
        from eefinder.hmm_analysis import score_for_evalue

        assert score_for_evalue(0.0, 0.0, 1e-5, 1e6) == 0.0


class TestWriteGaThresholds:
    def test_host_score_plus_margin_becomes_the_threshold(self, tmp_path):
        from eefinder.hmm_analysis import has_ga_thresholds, write_ga_thresholds

        database = _database(tmp_path, ["PROF_A", "PROF_B"])
        out = tmp_path / "calibrated.hmm"
        assert (
            write_ga_thresholds(database, {"PROF_A": 100.0}, str(out), margin=0.1) == 2
        )

        text = out.read_text()
        assert "GA    110.00 110.00;\n" in text
        assert has_ga_thresholds(str(out))

    def test_profiles_without_a_host_hit_fall_back_to_the_evalue_score(self, tmp_path):
        """A calibrated database must never be laxer than an uncalibrated one."""
        from eefinder.hmm_analysis import write_ga_thresholds

        database = _database(tmp_path, ["PROF_A"])
        out = tmp_path / "calibrated.hmm"
        write_ga_thresholds(database, {}, str(out), evalue=1e-5, z=1e6)
        threshold = float(
            next(l for l in out.read_text().splitlines() if l.startswith("GA")).split()[
                1
            ]
        )
        assert threshold == pytest.approx(30.98, abs=0.01)

    def test_a_weak_host_score_does_not_lower_the_floor(self, tmp_path):
        from eefinder.hmm_analysis import write_ga_thresholds

        database = _database(tmp_path, ["PROF_A"])
        out = tmp_path / "calibrated.hmm"
        write_ga_thresholds(database, {"PROF_A": 5.0}, str(out))
        threshold = float(
            next(l for l in out.read_text().splitlines() if l.startswith("GA")).split()[
                1
            ]
        )
        assert threshold > 5.0

    def test_existing_thresholds_are_replaced_not_duplicated(self, tmp_path):
        from eefinder.hmm_analysis import write_ga_thresholds

        database = _database(tmp_path, ["PROF_A"], ga="GA    12.00 12.00;\n")
        out = tmp_path / "calibrated.hmm"
        write_ga_thresholds(database, {"PROF_A": 200.0}, str(out))
        text = out.read_text()
        assert text.count("GA ") == 1
        assert "GA    220.00 220.00;\n" in text

    def test_the_threshold_precedes_the_stats_block(self, tmp_path):
        """HMMER expects GA in the header, before STATS/HMM."""
        from eefinder.hmm_analysis import write_ga_thresholds

        database = _database(tmp_path, ["PROF_A"])
        out = tmp_path / "calibrated.hmm"
        write_ga_thresholds(database, {"PROF_A": 50.0}, str(out))
        lines = out.read_text().splitlines()
        assert lines.index("GA    55.00 55.00;") < next(
            i for i, line in enumerate(lines) if line.startswith("STATS")
        )
