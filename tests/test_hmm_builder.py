"""Unit tests for the profile builder (no external binaries required).

The alignment/curation/build steps shell out to mafft/CIAlign/hmmbuild, so the
tests here cover the parts that decide *what* gets built: the profile-id rules
(which the pipeline's substring-based taxonomy lookup depends on), the curation
guards, and the builder-settings comparison the host filter relies on.
"""

from __future__ import annotations

import pytest

from eefinder.get_taxonomy import GetFinalTaxonomy  # noqa: F401  (documents the rule)
from eefinder.hmm_builder import (
    BuilderSettings,
    ClusterSpec,
    alignment_dimensions,
    build_profile,
    concat_profiles,
    curate_alignment,
    make_profile_id,
    parse_cialign_removed,
    read_hmm_header,
    slugify,
    write_member_fasta,
)


class TestProfileIds:
    def test_slug_is_id_safe(self):
        assert slugify("Glycoprotein (G)") == "Glycoprotein_G"
        assert slugify("RNaseH/reverse transcriptase") == "RNaseH_reverse_transcriptase"
        assert slugify("  ") == "Unknown"

    def test_id_has_no_separator_or_whitespace(self):
        """``|`` would be truncated by bed.RemoveAnnotation, whitespace by domtbl."""
        profile_id = make_profile_id("NCBIREFSEQ", "Chuviridae", "Glyco|protein X", 7)
        assert profile_id == "NCBIREFSEQ__Chuviridae__Glyco_protein_X__007"
        assert "|" not in profile_id
        assert not any(char.isspace() for char in profile_id)

    def test_no_id_is_a_substring_of_another(self):
        """GetFinalTaxonomy resolves accessions with ``in``, so ids must not nest."""
        ids = [
            make_profile_id("NCBIREFSEQ", "Chuviridae", "Glycoprotein", index)
            for index in (1, 10, 100)
        ]
        for candidate in ids:
            others = [other for other in ids if other != candidate]
            assert not any(candidate in other for other in others)

    def test_host_ids_omit_the_taxon(self):
        assert make_profile_id("HOST", "", "Actin", 3) == "HOST__Actin__003"


class TestBuilderSettings:
    def test_identical_settings_have_no_differences(self):
        assert BuilderSettings().differences(BuilderSettings()) == []

    def test_differences_name_the_offending_settings(self):
        viral = BuilderSettings(cluster_identity=0.5, min_profile_seqs=3)
        host = BuilderSettings(cluster_identity=0.7, min_profile_seqs=5)
        assert set(viral.differences(host)) == {"cluster_identity", "min_profile_seqs"}

    def test_hmmer_version_is_compared_too(self):
        assert BuilderSettings(hmmer_version="3.4").differences(
            BuilderSettings(hmmer_version="3.3")
        ) == ["hmmer_version"]


class TestMemberFasta:
    def test_members_are_written_sorted_for_determinism(self, tmp_path):
        out = tmp_path / "cluster.faa"
        written = write_member_fasta(
            {"B": "MMM", "A": "KKK", "C": "LLL"}, ["C", "A", "B"], str(out)
        )
        assert written == 3
        assert out.read_text() == ">A\nKKK\n>B\nMMM\n>C\nLLL\n"

    def test_missing_members_are_skipped(self, tmp_path):
        out = tmp_path / "cluster.faa"
        assert write_member_fasta({"A": "KKK"}, ["A", "ghost"], str(out)) == 1


class TestCuration:
    def _write_alignment(self, path, records):
        path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        return str(path)

    def test_alignment_dimensions(self, tmp_path):
        path = self._write_alignment(
            tmp_path / "aln.fa", [("A", "MK-LL"), ("B", "MKALL")]
        )
        assert alignment_dimensions(path) == (2, 5)

    def test_curation_is_skipped_below_three_sequences(self, tmp_path):
        """CIAlign's column-wise functions refuse to run on fewer than three."""
        path = self._write_alignment(
            tmp_path / "aln.fa", [("A", "MKALL"), ("B", "MKALL")]
        )
        result, removed_seqs, removed_cols = curate_alignment(
            path, str(tmp_path / "stem"), BuilderSettings()
        )
        assert result == path
        assert (removed_seqs, removed_cols) == (0, 0)

    def test_curation_is_skipped_when_every_function_is_off(self, tmp_path):
        path = self._write_alignment(
            tmp_path / "aln.fa", [("A", "MKALL"), ("B", "MKALL"), ("C", "MKALL")]
        )
        settings = BuilderSettings(
            curation_crop_ends=False,
            curation_remove_insertions=False,
            curation_remove_short=False,
            curation_remove_divergent=False,
        )
        result, _, _ = curate_alignment(path, str(tmp_path / "stem"), settings)
        assert result == path

    def test_removed_report_is_parsed(self, tmp_path):
        removed = tmp_path / "stem_removed.txt"
        removed.write_text("remove_insertions\t60,61,62\nremove_short\tSEQ_A\n")
        parsed = parse_cialign_removed(str(removed))
        assert parsed["remove_insertions"] == ["60", "61", "62"]
        assert parsed["remove_short"] == ["SEQ_A"]

    def test_missing_removed_report_is_not_an_error(self, tmp_path):
        assert parse_cialign_removed(str(tmp_path / "absent.txt")) == {}


class TestBuildProfileGuards:
    def test_cluster_below_the_minimum_is_skipped(self, tmp_path):
        spec = ClusterSpec(profile_id="P__X__001", representative="A", members=["A"])
        assert (
            build_profile({"A": "MKALL"}, spec, str(tmp_path), BuilderSettings())
            is None
        )

    def test_singletons_are_buildable_only_when_requested(self, tmp_path, monkeypatch):
        """``--keep-singletons`` lowers the threshold for one-member clusters."""
        spec = ClusterSpec(profile_id="P__X__001", representative="A", members=["A"])
        settings = BuilderSettings(keep_singletons=True, min_alignment_columns=1000)
        # The column guard rejects it afterwards, which proves the size threshold
        # was passed without needing hmmbuild on PATH.
        assert build_profile({"A": "MKALL"}, spec, str(tmp_path), settings) is None
        assert (tmp_path / "P__X__001.faa").is_file()


class TestConcatAndHeaders:
    def test_profiles_are_concatenated_in_order(self, tmp_path):
        paths = []
        for name in ("A", "B"):
            path = tmp_path / f"{name}.hmm"
            path.write_text(f"HMMER3/f\nNAME  {name}\nLENG  5\nHMM\n//\n")
            paths.append(str(path))
        out = tmp_path / "db.hmm"
        assert concat_profiles(paths, str(out)) == 2
        assert out.read_text().count("HMMER3/f") == 2

    def test_hmm_header_fields_are_read(self, tmp_path):
        path = tmp_path / "profile.hmm"
        path.write_text(
            "HMMER3/f [3.4 | Aug 2023]\nNAME  FAM__Capsid__001\nLENG  145\n"
            "NSEQ  6\nHMM   A C\n"
        )
        header = read_hmm_header(str(path))
        assert header["NAME"] == "FAM__Capsid__001"
        assert header["LENG"] == "145"
        assert header["NSEQ"] == "6"


@pytest.mark.parametrize(
    "identity,word_size",
    [(1.0, 5), (0.7, 5), (0.65, 4), (0.5, 3), (0.45, 2)],
)
def test_cdhit_word_size(identity, word_size):
    from eefinder.translation import cdhit_word_size

    assert cdhit_word_size(identity) == word_size


def test_cdhit_word_size_rejects_unusable_identity():
    from eefinder.translation import cdhit_word_size

    with pytest.raises(ValueError):
        cdhit_word_size(0.3)
